package tdg.trees;

import com.beust.jcommander.internal.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import pal.tree.Node;
import pal.tree.Tree;
import pal.tree.TreeManipulator;
import tdg.utils.PhyloUtils;

import java.util.*;

/**
 * Takes two trees and compares their branch lengths. They can be rooted differently but must have identical taxa,
 * which is identified by its string label.
 */
public class BranchLengthComparer {
    public static void main(String[] args) {
        BranchLengthComparer blc = new BranchLengthComparer();
        blc.run(args[0], args[1]);
    }

    private void run(String tree1Path, String tree2Path) {
        Tree tree1 = TreeManipulator.getUnrooted(PhyloUtils.readTree(tree1Path));
        Tree tree2 = TreeManipulator.getUnrooted(PhyloUtils.readTree(tree2Path));

        Map<Set<String>, Node> tree1Nodes = getNodeLeaves(tree1);
        Map<Set<String>, Node> tree2Nodes = getNodeLeaves(tree2);

        // Sanity check - the trees to have the same topology
        System.out.printf("Tree 1: %s nodes.\n", tree1Nodes.size());
        System.out.printf("Tree 2: %s nodes.\n", tree2Nodes.size());
        System.out.printf("Intersection: %s nodes.\n", Sets.intersection(tree1Nodes.keySet(), tree2Nodes.keySet()).size());
        System.out.printf("Difference: %s nodes.\n", Sets.difference(tree1Nodes.keySet(), tree2Nodes.keySet()).size());
        System.out.println();


        List<Set<String>> keys = Lists.newArrayList(tree1Nodes.keySet());
        Collections.sort(keys, new Comparator<Set<String>>() {
            @Override
            public int compare(Set<String> s1, Set<String> s2) {
                return s1.iterator().next().compareTo(s2.iterator().next());
            }
        });


        for (Set<String> group : keys) {
            Node node1 = tree1Nodes.get(group);
            Node node2 = tree2Nodes.get(group);

            System.out.printf("%s\t%s\t%s\n", node1.getBranchLength(), node2.getBranchLength(), group);
        }
    }

    private Map<Set<String>, Node> getNodeLeaves(Tree tree) {
        Map<Set<String>, Node> treeNodes = Maps.newHashMap();

        // Get taxa connect to internal nodes
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            Set<String> leaves = Sets.newTreeSet();
            getLeavesFromNode(tree.getInternalNode(i), leaves);
            treeNodes.put(leaves, tree.getInternalNode(i));
        }

        // Get the external nodes themselves
        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            Set<String> leaves = Sets.newTreeSet();
            getLeavesFromNode(tree.getExternalNode(i), leaves);
            treeNodes.put(leaves, tree.getExternalNode(i));
        }

        return treeNodes;
    }

    private void getLeavesFromNode(Node node, Set<String> leaves) {
        if (node.isLeaf()) {
            leaves.add(node.getIdentifier().getName());
        } else {
            for (int i = 0; i < node.getChildCount(); i++) {
                Node child = node.getChild(i);
                getLeavesFromNode(child, leaves);
            }
        }
    }

}
