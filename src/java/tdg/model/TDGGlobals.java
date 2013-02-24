package tdg.model;

import com.google.common.primitives.Doubles;
import tdg.utils.GeneticCode;

/**
 * @author Asif Tamuri (atamuri@nimr.mrc.ac.uk)
 */
public class TDGGlobals {
    private final double tau, kappa, mu, nu, gamma;
    private final double[] pi;
    private final double[] codonPiProduct = new double[GeneticCode.CODON_STATES];
    private final double[] neutralMu = new double[GeneticCode.CODON_STATES * GeneticCode.CODON_STATES];
    // private final double[] errorMatrix = new double[GeneticCode.CODON_STATES * GeneticCode.CODON_STATES];

    public TDGGlobals(double tau, double kappa, double[] pi, double mu, double gamma) { 
        this.tau = tau;
        this.kappa = kappa;
        this.pi = pi;
        this.mu = mu;
        this.gamma = gamma;
        calculateCodonPiProduct();
        this.nu = calculateNeutralMutationRate();
        makeErrorMatrix();
    }

    private void makeErrorMatrix() {
        /*
        Arrays.fill(errorMatrix, gamma);
        double diagonal = 1 - (GeneticCode.CODON_STATES * gamma);
        for (int i = 0; i < GeneticCode.CODON_STATES; i++) {
            errorMatrix[i * GeneticCode.CODON_STATES + i] = diagonal;
        }
        */

        // trying a different style of error matrix, with no_of_substitutions^gamma
       /* for (int i = 0; i < GeneticCode.CODON_STATES; i++) {
            double diagonal = 0;
            for (int j = 0; j < GeneticCode.CODON_STATES; j++) {
                if (i == j) continue;
                char[] codonCharsI = GeneticCode.getInstance().getNucleotidesFromCodonIndex(i);
                char[] codonCharsJ = GeneticCode.getInstance().getNucleotidesFromCodonIndex(j);

                int changes = 0;
                for (int k = 0; k < 3; k++) {
                    if (codonCharsI[k] != codonCharsJ[k]) {
                        changes++;
                    }
                }

                errorMatrix[i * GeneticCode.CODON_STATES + j] = Math.pow(gamma, changes);
                diagonal += Math.pow(gamma, changes);
            }
            errorMatrix[i * GeneticCode.CODON_STATES + i] = 1 - diagonal;
        }
*/
        /*
        for (int i = 0; i < GeneticCode.CODON_STATES; i++) {
            for (int j = 0; j < GeneticCode.CODON_STATES; j++) {
                System.out.printf("%s\t", errorMatrix[i * GeneticCode.CODON_STATES + j]);
            }
            System.out.println();
        }
        */

    }

    private void calculateCodonPiProduct() {
        for (int i = 0; i < GeneticCode.CODON_STATES; i++) {
            char[] nuc = GeneticCode.getInstance().getCodonTLA(i).toCharArray();
            codonPiProduct[i] =
                    pi[GeneticCode.getInstance().getNucleotideIndexByChar(nuc[0])] *
                    pi[GeneticCode.getInstance().getNucleotideIndexByChar(nuc[1])] *
                    pi[GeneticCode.getInstance().getNucleotideIndexByChar(nuc[2])];
        }
    }

    private double calculateNeutralMutationRate() {
        double nuDenominator = 0.0;
        for (int codonI = 0; codonI < GeneticCode.CODON_STATES; codonI++) {
            for (int codonJ = 0; codonJ < GeneticCode.CODON_STATES; codonJ++) {
                if (codonI == codonJ) continue;

                char[] codonCharsI = GeneticCode.getInstance().getNucleotidesFromCodonIndex(codonI);
                char[] codonCharsJ = GeneticCode.getInstance().getNucleotidesFromCodonIndex(codonJ);

                double prod_pi = 1.0;
                int changes = 0;
                int transitions = 0;

                // For each nucleotide
                for (int k = 0; k < 3; k++) {
                    // If the nucleotide at this position changes
                    if (codonCharsI[k] != codonCharsJ[k]) {
                        prod_pi *= pi[GeneticCode.getInstance().getNucleotideIndexByChar(codonCharsJ[k])];
                        changes++;
                    }
                    // Total number of transitions
                    if (GeneticCode.getInstance().isTransitionByChar(codonCharsI[k], codonCharsJ[k])) {
                        transitions++;
                    }
                }

                neutralMu[codonI * GeneticCode.CODON_STATES + codonJ] = Math.pow(tau, changes - 1) * Math.pow(kappa, transitions) * prod_pi;
                nuDenominator += neutralMu[codonI * GeneticCode.CODON_STATES + codonJ] * getCodonPiProduct(codonI);
            }
        }
        // System.out.printf("nu = %s\n", 1 / nuDenominator);
        return 1 / nuDenominator;
    }

    public double getNeutralMutationRate(int codonI, int codonJ) {
        return neutralMu[codonI * GeneticCode.CODON_STATES + codonJ];
    }

    public double getTau() {
        return tau;
    }

    public double getKappa() {
        return kappa;
    }

    public double getMu() {
        return mu;
    }

    public double[] getPi() {
        return pi;
    }

    public double getNu() {
        return nu;
    }

    public double getCodonPiProduct(int codon) {
        return codonPiProduct[codon];
    }

    @Override
    public String toString() {
        return "TDGGlobals{" +
                "tau=" + tau +
                ", kappa=" + kappa +
                ", mu=" + mu +
                ", 1/nu=" + nu +
                ", pi={" + Doubles.join(", ", pi) + "}" +
                '}';
    }
}
