package tdg.model;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.primitives.Doubles;
import pal.tree.Node;
import pal.tree.Tree;
import tdg.Constants;
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
    private boolean useScaling = true;

    private int[] siteCodons;

    private final double[][] tipConditionals;
    private final double[][] internalConditionals;

    public String hostshiftName;

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
        double lnL = calculateLogLikelihood();
        //double p = calculatePrior();

        /* if (!useScaling && l == Double.NEGATIVE_INFINITY) {
            useScaling = true;
            System.out.println("Turning on scaling.");
            l = calculateLogLikelihood();
        }*/

        return lnL; //(l) or (l - p)
    }

    private double calculatePrior() {
        return calculatePlainFitnessPrior();
        //return calculateTrueFitnessEnt();
    }

    private double calculateTrueFitnessEnt() {
        // MIT
        double[] baseAApi = new double[]{0.018207581862512204, 0.018207581862512204, 0.11236702136442346, 0.014633089085757604,
                0.0056780258819048, 0.07562774261045106, 0.017241116252006836, 0.004150844977304972,
                0.06418769404473233, 0.04360137848137165, 0.07418592860773728, 0.13239397825585356,
                0.05137236783972999, 0.016918488915981182, 0.0798671208133934, 0.06888517921966052,
                0.1398154366551834, 0.0066900094531035185, 0.04360137848137165, 0.012368035335008317};
        double[] pi = cladeModels.get(ROOT_MODEL_NAME).getAminoAcidFrequencies();

        double entropy = 0;

        for (int i = 0; i < GeneticCode.AMINO_ACID_STATES; i++) {
            entropy += pi[i] * Math.log(pi[i] / baseAApi[i]);
        }
        //System.out.printf("%s\n", entropy);
        return entropy;
    }

    private double calculatePlainFitnessPrior() {
        // maximum entropy prior
        //double[] f = cladeModels.get(ROOT_MODEL_NAME).getAminoAcidFrequencies();
        double[] f = cladeModels.get(ROOT_MODEL_NAME).getFitness();
        double[] pi = new double[20];

        double expFSum = 0;
        for (double F : f) {
            expFSum += Math.exp(F);
        }
        for (int i = 0; i < GeneticCode.AMINO_ACID_STATES; i++) {
            pi[i] = Math.exp(f[i]) / expFSum;
        }


        // the distribution should have unit variance...perhaps?
        // from: http://en.wikipedia.org/wiki/Prior_probability#Other_approaches
        /*
        Another idea, championed by Edwin T. Jaynes, is to use the principle of maximum entropy (MAXENT). The motivation
        is that the Shannon entropy of a probability distribution measures the amount of information contained in the
        distribution. The larger the entropy, the less information is provided by the distribution. Thus, by maximizing
        the entropy over a suitable set of probability distributions on X, one finds that distribution that is least
        informative in the sense that it contains the least amount of information consistent with the constraints that
        define the set. For example, the maximum entropy prior on a discrete space, given only that the probability is
        normalized to 1, is the prior that assigns equal probability to each state. And in the continuous case, the
        maximum entropy prior given that the density is normalized with mean zero and variance unity is the standard
        normal distribution. The principle of minimum cross-entropy generalizes MAXENT to the case of "updating" an
        arbitrary prior distribution with suitability constraints in the maximum-entropy sense.
         */

/*
        double mean = 0;
        for (double P : pi) {
            mean += P;
        }
        mean /= pi.length;
        for (int i = 0; i < GeneticCode.AMINO_ACID_STATES; i++) {
            pi[i] = pi[i] - mean;
        }


        DescriptiveStatistics stats = new DescriptiveStatistics();
            for( double item : pi) {
                stats.addValue(item);
            }
        double std = stats.getStandardDeviation();

        for (int i = 0; i < GeneticCode.AMINO_ACID_STATES; i++) {
            pi[i] /= std;
        }
*/

        double sum = 0;
        //System.out.printf("pi = %s\n", Doubles.join(" ", pi));
        for (double P : pi) {
            if (P > 0) sum += P * Math.log(P);
        }


        // normal distribution prior
        /*double sigma = 10;
        double[] f = cladeModels.get(ROOT_MODEL_NAME).getFitness();
        double sum = 0;
        for (double F : f) {
            sum += Math.pow(F, 2) / Math.pow(sigma, 2);
        }*/

        return sum;
    }

    private double calculateLogLikelihood() {
        logScaling = 0.0;
        double[] conditionals = downTree();
        //     System.out.printf("conditionals out = %s\n", Doubles.join(",", conditionals));

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

                if (cladeModels.size() == 1) { // homogeneous model

                    cladeModels.get(ROOT_MODEL_NAME).getProbabilityMatrix(probMatrix, child.getBranchLength());
                    updateIntraCladeConditionals(lowerConditional, partial, probMatrix);

                } else { // non-homogeneous model

                    if (getNodeLabel(node).length() == 0 // the root of the tree is a parent without a label
                            || getNodeLabel(child).substring(0, 2).equals(getNodeLabel(node).substring(0, 2))
                            || getNodeLabel(node).substring(0, 2).equals("HS") && getNodeLabel(child).substring(0, 2).equals(this.hostshiftName)) { // or we're not switching to a different model

                        cladeModels.get(getNodeLabel(child).substring(0, 2)).getProbabilityMatrix(probMatrix, child.getBranchLength());
                        updateIntraCladeConditionals(lowerConditional, partial, probMatrix);

                    } else { // this is a hostshift!

                        cladeModels.get(getNodeLabel(node).substring(0, 2)).getProbabilityMatrix(probMatrix0, child.getBranchLength() * Constants.CLADE_BRANCH_SPLIT[0]);
                        cladeModels.get(getNodeLabel(child).substring(0, 2)).getProbabilityMatrix(probMatrix1, child.getBranchLength() * Constants.CLADE_BRANCH_SPLIT[1]);
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

    private void scaleConditionals(Node node, double[] conditionals) {

        if (node.getNumber() % Constants.SCALING_NODE_STEP == 0) {
            double scalingFactor = 0;
            for (int i = 0; i < conditionals.length; i++) {
                if (conditionals[i] > 0 && conditionals[i] > scalingFactor) {
                    scalingFactor = conditionals[i];
                }
            }

            if (scalingFactor < Constants.SCALING_THRESHOLD) {
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
                double[] tmp = Doubles.concat(new double[]{Constants.FITNESS_INITIAL_VALUE}, Arrays.copyOfRange(params, offset, offset + len - 1));
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
                // First fitness is fixed to FITNESS_INITIAL_VALUE ( = 0.0)
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
