package tdg.trees;

import pal.tree.Node;
import pal.tree.SimpleNode;
import pal.tree.SimpleTree;
import pal.tree.TreeUtils;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

/**
 * @author Asif Tamuri (atamuri@nimr.mrc.ac.uk)
 */
public class TreeMaker {
    public static void main(String[] args) {
        String fileIn = "/Users/atamuri/Documents/work/mitochondria/paper/response/consistency/typical.AA.normal.F/tdg.8192/tree.8192.tree";

        SimpleTree st = new SimpleTree();

        Node root = new SimpleNode();
        st.setRoot(root);

        // 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048
        add2Nodes(root);

        st.createNodeList();

        for (int i = 0; i < st.getExternalNodeCount(); i++) {
            st.getExternalNode(i).getIdentifier().setName("s" + i);
        }

        FileWriter fw;
        try {
            fw = new FileWriter(fileIn + ".out");
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            throw new RuntimeException();
        }

        PrintWriter pw = new PrintWriter(fw);
        TreeUtils.printNH(st, pw);
        pw.flush();
    }

    private static void add2Nodes(Node root) {
        System.out.printf("root.getNodeHeight() = %s\n", root.getNodeHeight());

        Node n1 = new SimpleNode("", 0.25);
        n1.setNodeHeight(root.getNodeHeight() + 0.25);
        root.addChild(n1);
        if (n1.getNodeHeight() < 3.25) add2Nodes(n1);

        Node n2 = new SimpleNode("", 0.25);
        n2.setNodeHeight(root.getNodeHeight() + 0.25);
        root.addChild(n2);
        if (n2.getNodeHeight() < 3.25) add2Nodes(n2);
        


    }
}
