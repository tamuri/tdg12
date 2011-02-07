package tdg.models;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.primitives.Doubles;
import pal.tree.Node;
import pal.tree.Tree;
import tdg.models.parameters.Fitness;
import tdg.models.parameters.Parameter;
import tdg.utils.GeneticCode;

import java.util.Arrays;
import java.util.List;
import java.util.Map;

/**
 * @author Asif Tamuri
 * @version $Id: LikelihoodCalculator.java 152 2010-11-08 11:10:01Z tamuri $
 */
public class LikelihoodCalculator {
    private static final double[] CLADE_BRANCH_SPLIT = {0.5, 0.5};
    private static final double SCALING_THRESHOLD = 1e-70;
    private static final int SCALING_NODE_STEP = 10;
    private static final double FITNESS_INITIAL_VALUE = 0.0;

    private final Tree tree;
    private String ROOT_MODEL_NAME; 
    private final Map<String, Integer> states;
    private final Map<String, TDGCodonModel> cladeModels = Maps.newHashMap();
    private final double[] probMatrix = new double[GeneticCode.CODON_STATES * GeneticCode.CODON_STATES];
    private final double[] probMatrix0 = new double[GeneticCode.CODON_STATES * GeneticCode.CODON_STATES];
    private final double[] probMatrix1 = new double[GeneticCode.CODON_STATES * GeneticCode.CODON_STATES];

    private List<Parameter> parameters;
    private final Map<Node, String> nodeLabels = Maps.newHashMap();
    private double logScaling = 0.0;
    private boolean useScaling = true;

    private int[] siteCodons;

    private final double[][] tipConditionals;
    private final double[][] internalConditionals;
    
    public LikelihoodCalculator(Tree tree, Map<String, Integer> states) {
        this.tree = tree;
        this.states = states;

        // set up name lookup
        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            Node n = tree.getExternalNode(i);
            nodeLabels.put(n, n.getIdentifier().getName());
        }
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            Node n = tree.getInternalNode(i);
            nodeLabels.put(n, n.getIdentifier().getName());
        }

        this.tipConditionals = new double[tree.getExternalNodeCount()][GeneticCode.CODON_STATES];
        this.internalConditionals = new double[tree.getInternalNodeCount()][GeneticCode.CODON_STATES];
        fillTipConditionals();
    }

    private void fillTipConditionals() {
        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            String parentName = getNodeLabel(tree.getExternalNode(i));
            int codon = states.get(parentName);

            if (GeneticCode.getInstance().isUnknownCodonState(codon)) {
                Arrays.fill(tipConditionals[i], 1.0);
            } else {
                tipConditionals[i][codon] = 1.0;
            }
        }
    }

    public void applyErrorMatrix(double[] errorMatrix) {
        // TODO: make this an option...better way to handle this!!
        double[] cc = new double[GeneticCode.CODON_STATES];
        for (int i : siteCodons) {
            cc[i] = 1.0;
        }


        for (int k = 0; k < tree.getExternalNodeCount(); k++) {

            double[] c = cc.clone();

            for (int i : siteCodons) {
                double branchProb = 0.0;
                for (int j : siteCodons) {
                    branchProb += tipConditionals[k][j] * errorMatrix[i * GeneticCode.CODON_STATES + j];
                }

                c[i] *= branchProb;
            }

            for (int i : siteCodons) {
                // System.out.printf("%s (%s). %s -> %s\n", k, i, tipConditionals[k][i], c[i]);
                tipConditionals[k][i] = c[i];

            }
        }
    }

    private String getNodeLabel(Node n) {
        return nodeLabels.get(n);
    }

    public double function(double[] parameters) {
        updateParameters(parameters);
        double l = calculateLogLikelihood();

       /* if (!useScaling && l == Double.NEGATIVE_INFINITY) {
            useScaling = true;
            System.out.println("Turning on scaling.");
            l = calculateLogLikelihood();
        }*/

       //System.out.printf("function lnL eval = %s\n", -1.0 * l);
//System.exit(0);
        return -1.0 * l;
    }

    private double calculateLogLikelihood() {
        logScaling = 0.0;
        double[] conditionals = downTree();
  //      System.out.printf("conditionals out = %s\n", Doubles.join(",", conditionals));

        double[] f = cladeModels.get(ROOT_MODEL_NAME).getCodonFrequencies();
    //    System.out.printf("codon freqs = %s\n", Doubles.join(",", f));
      //  double sss = 0.;
        //for (double dd:f) sss+=dd;
        //System.out.printf("codon sum = %s\n", sss);

        double sum = 0.0;

        for (int i = 0; i < GeneticCode.CODON_STATES; i++) {
            sum += conditionals[i] * f[i];
        }

        if (sum < 0) sum = 0;

        //System.out.printf("sum = %s, logScaling = %s\n", Math.log(sum), logScaling);

        return Math.log(sum) + logScaling;
    }

    private double[] downTree() {
        //long start = CodeTimer.start();
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {

            Node node = tree.getInternalNode(i);

            double[] partial = new double[GeneticCode.CODON_STATES];
            for (int j : siteCodons) partial[j] = 1.0;

            for (int j = 0; j < node.getChildCount(); j++) {
                Node child = node.getChild(j);

                double[] lowerConditional;

                if (child.isLeaf()) {
                    lowerConditional = tipConditionals[child.getNumber()];
                    //System.out.printf("child %s (tip = %s, %s, %s)\n", Doubles.join(", ", child.getNumber()), states.get(getNodeLabel(child)),
                    //        GeneticCode.getCodonTLA(states.get(getNodeLabel(child))), GeneticCode.getAminoAcidCharByIndex(GeneticCode.getAminoAcidIndexFromCodonIndex(states.get(getNodeLabel(child)))));
                } else {
                    lowerConditional = internalConditionals[child.getNumber()];
                    //System.out.printf("child %s (not tip)\n", Doubles.join(", ", child.getNumber()));
                }

                if (cladeModels.size() == 1) {

                    cladeModels.get(ROOT_MODEL_NAME).getProbabilityMatrix(probMatrix, child.getBranchLength());
                    updateIntraCladeConditionals(lowerConditional, partial, probMatrix);

                } else {

                    if (getNodeLabel(node).length() == 0 // or the root of the tree is a parent without a label
                            || getNodeLabel(child).substring(0, 2).equals(getNodeLabel(node).substring(0, 2))) { // or we're not switching to a different model

                        cladeModels.get(getNodeLabel(child).substring(0, 2)).getProbabilityMatrix(probMatrix, child.getBranchLength());
                        updateIntraCladeConditionals(lowerConditional, partial, probMatrix);

                    } else {

                        cladeModels.get(getNodeLabel(node).substring(0, 2)).getProbabilityMatrix(probMatrix0, child.getBranchLength() * CLADE_BRANCH_SPLIT[0]);
                        cladeModels.get(getNodeLabel(child).substring(0, 2)).getProbabilityMatrix(probMatrix1, child.getBranchLength() * CLADE_BRANCH_SPLIT[1]);
                        updateInterCladeConditionals(lowerConditional, partial, probMatrix0, probMatrix1);

                    }
                }

            }

            if (useScaling) {
                scaleConditionals(node, partial);
            }

            internalConditionals[node.getNumber()] = partial;
            //System.out.printf("%s = %s\n", node.getNumber(), Doubles.join(", ", partial));
        }
        //CodeTimer.store("downTree", start);
        return internalConditionals[tree.getRoot().getNumber()];
    }

    private double[] downTree(Node parent) {
        double[] conditionals = new double[GeneticCode.CODON_STATES];
        String parentName, childName;
        double[] lowerConditional;

        parentName = getNodeLabel(parent);

        if (parent.isLeaf()) {
            int codon = states.get(parentName);

            if (GeneticCode.getInstance().isUnknownCodonState(codon)) {
                Arrays.fill(conditionals, 1.0);
            } else {
                conditionals[codon] = 1.0;
            }
        } else {
            // this is an internal node - initialise the conditionals array for this node
            Arrays.fill(conditionals, 1.0);

            // for each child from this node
            for (int i = 0; i < parent.getChildCount(); i++) {
                Node child = parent.getChild(i);

                // recurse down the tree to get the conditionals for this child
                lowerConditional = downTree(child);

                // if we have a single clade model
                if (cladeModels.size() == 1) {
                    // homogeneous model
                    cladeModels.get(ROOT_MODEL_NAME).getProbabilityMatrix(probMatrix, child.getBranchLength());
                    updateIntraCladeConditionals(lowerConditional, conditionals, probMatrix);
                } else {
                    // clade-specific model
                    childName = getNodeLabel(child);

                    if (parentName.length() == 0 // the root of the tree is a parent without a label
                            || childName.substring(0, 2).equals(parentName.substring(0, 2))) {
                        // we're on a branch in a common clade
                        cladeModels.get(childName.substring(0, 2)).getProbabilityMatrix(probMatrix, child.getBranchLength());
                        updateIntraCladeConditionals(lowerConditional, conditionals, probMatrix);

                    } else {
                        // we're on a branch switching clades
                        cladeModels.get(parentName.substring(0, 2)).getProbabilityMatrix(probMatrix0, child.getBranchLength() * CLADE_BRANCH_SPLIT[0]);
                        cladeModels.get(childName.substring(0, 2)).getProbabilityMatrix(probMatrix1, child.getBranchLength() * CLADE_BRANCH_SPLIT[1]);
                        updateInterCladeConditionals(lowerConditional, conditionals, probMatrix0, probMatrix1);
                    }
                }
            }
        }

        if (useScaling) {
            scaleConditionals(parent, conditionals);
        }

        return conditionals;
    }

    private void scaleConditionals(Node node, double[] conditionals) {

        if (node.getNumber() % SCALING_NODE_STEP == 0) {
            double scalingFactor = 0;
            for (int i = 0; i < conditionals.length; i++) {
                if (conditionals[i] > 0 && conditionals[i] > scalingFactor) {
                    scalingFactor = conditionals[i];
                }
            }

            if (scalingFactor < SCALING_THRESHOLD) {
            //    System.out.printf("need to scale, scalingFactor = %s\n", scalingFactor);
             //   System.out.printf("conditionals before scaling: %s\n", Doubles.join(",", conditionals));

                for (int i = 0; i < conditionals.length; i++) {
                    conditionals[i] = conditionals[i] / scalingFactor;
                }
               // System.out.printf("conditionals after scaling: %s\n", Doubles.join(",", conditionals));
                logScaling += Math.log(scalingFactor);
               // System.out.printf("log scaling = %s (total = %s)\n", Math.log(scalingFactor), logScaling);

            }
        }
    }

    private void updateInterCladeConditionals(double[] lowerConditional, double[] conditionals, double[] probMatrix0, double[] probMatrix1) {
        double[] probMatrix = new double[GeneticCode.CODON_STATES * GeneticCode.CODON_STATES];
        for (int i : siteCodons) {
            double branchProb = 0.0;
            for (int j : siteCodons) {
                for (int k : siteCodons) {
                    probMatrix[i * GeneticCode.CODON_STATES + j] += probMatrix0[i * GeneticCode.CODON_STATES + k] * probMatrix1[k * GeneticCode.CODON_STATES + j];
                }
                branchProb += lowerConditional[j] * probMatrix[i * GeneticCode.CODON_STATES + j];
            }
            conditionals[i] *= branchProb;
        }
    }

    private void updateIntraCladeConditionals(double[] lowerConditional, double[] conditionals, double[] probMatrix) {
        //System.out.printf("(%s\n%s)\n", Doubles.join(", ", lowerConditional), Doubles.join(", ", conditionals));
        for (int i : siteCodons) {
            double branchProb = 0.0;
            for (int j : siteCodons) {
                branchProb += lowerConditional[j] * probMatrix[i * GeneticCode.CODON_STATES + j];
            }
            conditionals[i] *= branchProb;
        }

       /* if (CoreUtils.sum(conditionals) == 0) {
            System.exit(0);
        }*/
    }

    private void updateParameters(double[] params) {
        int offset = 0;

        for (Parameter p : parameters) {
            if (p.getClass() == Fitness.class) {
                int len = ((double[]) p.get()).length; // number of fitness coefficient parameters
                double[] tmp = Doubles.concat(new double[]{FITNESS_INITIAL_VALUE}, Arrays.copyOfRange(params, offset, offset + len - 1));
                p.set(tmp);
                offset = offset + len - 1;
            }
        }

        // We've updated the parameters. Notify every clade model to create new models.
        for (String key : cladeModels.keySet()) {
            cladeModels.get(key).updateModel();
        }
    }

    public void setParameters(Parameter... parameters) {
        this.parameters = Arrays.asList(parameters);
    }

    public void addCladeModel(String cladeName, TDGCodonModel cladeModel) {
        // If this is the first clade we're adding, it becomes the default, the "root" model
        if (cladeModels.isEmpty()) {
            ROOT_MODEL_NAME = cladeName;
        }
        cladeModels.put(cladeName, cladeModel);
        siteCodons = cladeModels.get(ROOT_MODEL_NAME).getSiteCodons();
    }

    public MinimisationParameters getMinimisationParameters() {
        List<Double> minParams = Lists.newArrayList();
        for (Parameter p : parameters) {
            // TODO: currently only optimising fitness parameters
            if (p.getClass() == Fitness.class) {
                double[] fitness = (double[]) p.get();
                // first fitness is always 1.0
                minParams.addAll(Doubles.asList(fitness).subList(1, fitness.length));
            } else {
                throw new RuntimeException("NOT IMPLEMENTED!");
            }
        }

        return new MinimisationParameters(Doubles.toArray(minParams),
                new double[]{0.0},
                new double[]{0.0},
                new double[]{0.0});
    }


}
