package tdg.model;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.primitives.Doubles;
import pal.tree.Node;
import pal.tree.Tree;
import tdg.Constants;
import tdg.MatrixArrayPool;
import tdg.utils.GeneticCode;

import java.util.Arrays;
import java.util.List;
import java.util.Map;

/**
 * Felsenstein's pruning algorithm to calculate the likelihood for codon based models. Can deal with heterogenous models.
 *
 * @author Asif Tamuri (atamuri@nimr.mrc.ac.uk)
 */
public class LikelihoodCalculator {
    private final Tree tree;
    private String ROOT_MODEL_NAME;
    private final Map<String, Integer> states;
    private final Map<String, TDGCodonModel> cladeModels = Maps.newHashMap();
    // private final double[] probMatrix = new double[GeneticCode.CODON_STATES * GeneticCode.CODON_STATES];
    private double[] probMatrix;

    private Map<Node, double[]> rootPartials;

    private List<Parameter> parameters;
    private final Map<Node, String> nodeLabels = Maps.newHashMap();
    private double logScaling = 0.0;

    private int[] siteCodons;

    // TODO: Improve memory usage of these conditionals...do we need to store everything??!
    // TODO: For example, strictly speaking, we only need as many arrays as we have concurrent threads running, and each can be reused
    // Could we have a reusable pool of 64*64 arrays? We would then only create 16 (or whatever), instead of 3598!
   //  private final double[][] tipConditionals;
    // private final double[][] internalConditionals;
    private double[][] internalConditionals;

    private final double[] gapPartial = new double[GeneticCode.CODON_STATES];
    {
        Arrays.fill(gapPartial, 1.0);
    }

    private final double[] tipPartial = new double[GeneticCode.CODON_STATES];

    private Prior prior;

    public LikelihoodCalculator(Tree tree, Map<String, Integer> states, Prior prior) {
        this.tree = tree;
        this.states = states;

        if (prior != null) {
            this.prior = prior;
        }

        // set up name lookup
        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            Node n = tree.getExternalNode(i);
            nodeLabels.put(n, n.getIdentifier().getName());
        }
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            Node n = tree.getInternalNode(i);
            nodeLabels.put(n, n.getIdentifier().getName());
        }

        //this.tipConditionals = new double[tree.getExternalNodeCount()][GeneticCode.CODON_STATES];
        // this.internalConditionals = new double[tree.getInternalNodeCount()][GeneticCode.CODON_STATES];
        //fillTipConditionals();
    }

    private void fillTipConditionals() {
        // TODO: do we need to save this each time? we just need codon index or gap!
        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            String parentName = getNodeLabel(tree.getExternalNode(i));
            int codon = states.get(parentName);

            if (GeneticCode.getInstance().isUnknownCodonState(codon)) {
               // Arrays.fill(tipConditionals[i], 1.0);
            } else {
               // tipConditionals[i][codon] = 1.0;

/*                int aai = GeneticCode.getInstance().getAminoAcidIndexFromCodonIndex(codon);
                int[] codoni = GeneticCode.getInstance().getCodonIndexFromAminoAcidIndex(aai);

                for (int j : codoni) {
                    tipConditionals[i][j] = 1.0;
                }*/
            }
        }
    }

    private String getNodeLabel(Node n) {
        return nodeLabels.get(n);
    }

    public double function(double[] parameters) {



        updateParameters(parameters);
        double l = calculateLogLikelihood();

        double p = 0.0;

        // If prior has been specified
        if (prior != null) {
            p = prior.calculate(parameters);
        }

        // System.out.printf("%s\t %s\n", l+p, Doubles.join(" ", parameters));


        return l + p;
    }

    public double calculateLogLikelihood() {

        if (probMatrix == null) getStorage();

        rootPartials = Maps.newHashMap();
        rootpartials = Lists.newArrayList();
        logScaling = 0.0;
        double[] conditionals = downTree();
        double[] f = cladeModels.get(ROOT_MODEL_NAME).getCodonFrequencies();

        double sum = 0.0;
        for (int i = 0; i < GeneticCode.CODON_STATES; i++)
            sum += conditionals[i] * f[i];

        if (sum < 0) sum = 0;

        // if (logScaling != 0) System.out.printf("log-scaling = %s\n", logScaling);


        return Math.log(sum) + logScaling;
    }

    public void getStorage() {
        probMatrix = MatrixArrayPool.pop();
        internalConditionals = MatrixArrayPool.pop2();

    }

    public void releaseStorage() {
        MatrixArrayPool.push(probMatrix);
        MatrixArrayPool.push2(internalConditionals);
        probMatrix = null;
        internalConditionals = null;
    }

    public Map<Node, double[]> getRootPartials() {
        return rootPartials;
    }

    public Map<String, TDGCodonModel> getCladeModels() {
        return cladeModels;
    }

    public double getLogScaling() {
        return logScaling;
    }

    class Partial {
        public int number;
        public double[] partial;
        public double branchlength;

        Partial(int number, double[] partial, double branchlength) {
            this.number = number;
            this.partial = partial;
            this.branchlength = branchlength;
        }
    }

    private List<Partial> rootpartials = Lists.newArrayList();

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
                    if (GeneticCode.getInstance().isUnknownCodonState(states.get(getNodeLabel(child)))) {
                        lowerConditional = gapPartial;
                    } else {
                        Arrays.fill(tipPartial, 0.0);
                        tipPartial[states.get(getNodeLabel(child))] = 1.0;
                        lowerConditional = tipPartial;
                    }

                } else {
                    lowerConditional = internalConditionals[child.getNumber()];
                }

                if (child.getParent().isRoot()) {
                    rootPartials.put(child, Arrays.copyOf(lowerConditional, GeneticCode.CODON_STATES));
                    rootpartials.add(new Partial(child.getNumber(), Arrays.copyOf(lowerConditional, GeneticCode.CODON_STATES), child.getBranchLength()));
                    // System.out.printf("%s\t%s\n", child.getNumber(), Doubles.join(",", lowerConditional));
                }

                if (cladeModels.size() == 1) { // homogeneous model

                    // TODO: A probability matrix (of sorts) already exists in TdGCodonModel...do we need to create it again?
                    cladeModels.get(ROOT_MODEL_NAME).getProbabilityMatrix(probMatrix, child.getBranchLength());
                    updateIntraCladeConditionals(lowerConditional, partial, probMatrix);

                } else { // non-homogeneous model

                    if (getNodeLabel(node).length() == 0 // the root of the tree is a parent without a label
                            || getNodeLabel(child).substring(0, 2).equals(getNodeLabel(node).substring(0, 2))) { // or we're not switching to a different model

                        cladeModels.get(getNodeLabel(child).substring(0, 2)).getProbabilityMatrix(probMatrix, child.getBranchLength());
                        updateIntraCladeConditionals(lowerConditional, partial, probMatrix);

                    } else { // this is a hostshift!
                        final double[] probMatrix1 = new double[GeneticCode.CODON_STATES * GeneticCode.CODON_STATES];
                        cladeModels.get(getNodeLabel(node).substring(0, 2)).getProbabilityMatrix(probMatrix, child.getBranchLength() * Constants.CLADE_BRANCH_SPLIT);
                        cladeModels.get(getNodeLabel(child).substring(0, 2)).getProbabilityMatrix(probMatrix1, child.getBranchLength() * (1 - Constants.CLADE_BRANCH_SPLIT));
                        updateInterCladeConditionals(lowerConditional, partial, probMatrix, probMatrix1);
                    }
                }
            }

            if (!(node.getParent() == null) && !node.getParent().isRoot() ) {
                scaleConditionals(node, partial);
            }

            internalConditionals[node.getNumber()] = partial;

        }
        //CodeTimer.store("downTree", start);
        return internalConditionals[tree.getRoot().getNumber()];
    }

    public double getNodeLikelihood(int node, double branchLength) {
        // TODO: store the true branchlengths here, so they can be updated with the distribtuedrunner
        probMatrix = MatrixArrayPool.pop();

        double[] sumPartial = new double[GeneticCode.CODON_STATES];
        Arrays.fill(sumPartial, 1.0);
/*

        for (final Map.Entry<Node, double[]> e : getRootPartials().entrySet()) {
            if (e.getKey().getNumber() == node.getNumber()) {
                getCladeModels().get("ALL").getProbabilityMatrix(probMatrix, branchLength);
                // System.out.printf("found!");
            } else
                getCladeModels().get("ALL").getProbabilityMatrix(probMatrix, e.getKey().getBranchLength());

            updateIntraCladeConditionals(e.getValue(), sumPartial, probMatrix);
        }
*/


        for (Partial p : rootpartials) {
            if (p.number == node) {
                getCladeModels().get("ALL").getProbabilityMatrix(probMatrix, branchLength);
            } else {
                getCladeModels().get("ALL").getProbabilityMatrix(probMatrix, p.branchlength);
            }
            updateIntraCladeConditionals(p.partial, sumPartial, probMatrix);

        }

        double lnL = 0;
        final double[] frequencies = getCladeModels().get("ALL").getCodonFrequencies();
        for (int i = 0; i < GeneticCode.CODON_STATES; i++) lnL += sumPartial[i] * frequencies[i];

        if (lnL < 0) lnL = 0;

        lnL = Math.log(lnL) + getLogScaling();

        MatrixArrayPool.push(probMatrix);

        return lnL + prior.calculate(getMinimisationParameters().getParameters());
    }

    public void setBranch(int node, double bl) {
        for (Partial p : rootpartials) {
            if (p.number == node) {
                p.branchlength = bl;
            }
        }
    }

    private void scaleConditionals(Node node, double[] conditionals) {
        if (node.getNumber() % Constants.SCALING_NODE_STEP == 0) {
            double scalingFactor = 0;
            for (double conditional : conditionals) {
                if (conditional > 0 && conditional > scalingFactor) {
                    scalingFactor = conditional;
                }
            }

            if (scalingFactor < Constants.SCALING_THRESHOLD) {
                for (int i = 0; i < conditionals.length; i++) {
                    conditionals[i] = conditionals[i] / scalingFactor;
                }
                logScaling += Math.log(scalingFactor);
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

    public void updateIntraCladeConditionals(double[] lowerConditional, double[] conditionals, double[] probMatrix) {
        for (int i : siteCodons) {
            double branchProb = 0.0;
            for (int j : siteCodons) {
                branchProb += lowerConditional[j] * probMatrix[i * GeneticCode.CODON_STATES + j];
            }
            conditionals[i] *= branchProb;
        }
    }

    private void updateParameters(double[] params) {
        int offset = 0;

        for (Parameter p : parameters) {
            if (p.getClass() == Fitness.class) {
                int len = ((double[]) p.get()).length; // number of fitness coefficient parameters
                double[] tmp = Doubles.concat(new double[]{Constants.FITNESS_FIXED_FOR_RELATIVE}, Arrays.copyOfRange(params, offset, offset + len - 1));
                p.set(tmp);
                offset = offset + len - 1;
            }
        }

        // We've updated the parameters. Notify every clade model to create new model.
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
                // First fitness is fixed to FITNESS_FIXED_FOR_RELATIVE ( = 0.0)
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
