package tdg.model;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.primitives.Doubles;
import org.apache.commons.math.util.MathUtils;
import pal.tree.Node;
import pal.tree.Tree;
import tdg.Constants;
import tdg.cli.AnalyseOptions;
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
    private final double[] probMatrix = new double[GeneticCode.CODON_STATES * GeneticCode.CODON_STATES];
    private final double[] probMatrix0 = new double[GeneticCode.CODON_STATES * GeneticCode.CODON_STATES];
    private final double[] probMatrix1 = new double[GeneticCode.CODON_STATES * GeneticCode.CODON_STATES];

    private List<Parameter> parameters;
    private final Map<Node, String> nodeLabels = Maps.newHashMap();
    private double logScaling = 0.0;

    private int[] siteCodons;

    private final double[][] tipConditionals;
    private final double[][] internalConditionals;

    private Prior prior;

    public LikelihoodCalculator(Tree tree, Map<String, Integer> states, List<String> prior) {
        this.tree = tree;
        this.states = states;

        // Instatiate the prior calculator if necessary
        if (prior != null) {
            if (prior.get(0).equals("normal")) {
                this.prior = new NormalPrior(Double.parseDouble(prior.get(1)));
            } else if (prior.get(0).equals("dirichlet")) {
                this.prior = new DirichletPrior(Double.parseDouble(prior.get(1)));
            }
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

        return l + p;
    }

    private double calculateLogLikelihood() {
        logScaling = 0.0;
        double[] conditionals = downTree();
        double[] f = cladeModels.get(ROOT_MODEL_NAME).getCodonFrequencies();

        double sum = 0.0;
        for (int i = 0; i < GeneticCode.CODON_STATES; i++)
            sum += conditionals[i] * f[i];

        if (sum < 0) sum = 0;

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
                } else {
                    lowerConditional = internalConditionals[child.getNumber()];
                }

                if (cladeModels.size() == 1) { // homogeneous model

                    cladeModels.get(ROOT_MODEL_NAME).getProbabilityMatrix(probMatrix, child.getBranchLength());
                    updateIntraCladeConditionals(lowerConditional, partial, probMatrix);

                } else { // non-homogeneous model

                    if (getNodeLabel(node).length() == 0 // the root of the tree is a parent without a label
                            || getNodeLabel(child).substring(0, 2).equals(getNodeLabel(node).substring(0, 2))) { // or we're not switching to a different model

                        cladeModels.get(getNodeLabel(child).substring(0, 2)).getProbabilityMatrix(probMatrix, child.getBranchLength());
                        updateIntraCladeConditionals(lowerConditional, partial, probMatrix);

                    } else { // this is a hostshift!

                        cladeModels.get(getNodeLabel(node).substring(0, 2)).getProbabilityMatrix(probMatrix0, child.getBranchLength() * Constants.CLADE_BRANCH_SPLIT);
                        cladeModels.get(getNodeLabel(child).substring(0, 2)).getProbabilityMatrix(probMatrix1, child.getBranchLength() * (1 - Constants.CLADE_BRANCH_SPLIT));
                        updateInterCladeConditionals(lowerConditional, partial, probMatrix0, probMatrix1);

                    }
                }
            }

            if (Constants.USE_SCALING) {
                scaleConditionals(node, partial);
            }

            internalConditionals[node.getNumber()] = partial;
        }
        //CodeTimer.store("downTree", start);
        return internalConditionals[tree.getRoot().getNumber()];
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

    private void updateIntraCladeConditionals(double[] lowerConditional, double[] conditionals, double[] probMatrix) {
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
