package tdg.model;

import com.google.common.primitives.Doubles;
import org.apache.commons.math.util.MathUtils;
import tdg.utils.CoreUtils;
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

        if (pi.length != 4) throw new RuntimeException("-pi should have exactly 4 items.");
        this.pi = MathUtils.normalizeArray(pi, 1);

        this.mu = mu;
        this.gamma = gamma;
        calculateCodonPiProduct();
        this.nu = calculateNeutralMutationRate();
        makeErrorMatrix();
    }

    /**
     * Returns globals initialised with (arbitrary) "sensible" values (used as initial parameters value in MLE)
     */
    public TDGGlobals() {
        this.tau = 0.01;
        this.kappa = 2.0;
        this.pi = CoreUtils.repd(0.25, 4);
        this.mu = 1.0;
        this.gamma = 0;
        calculateCodonPiProduct();
        this.nu = calculateNeutralMutationRate();
        makeErrorMatrix();
    }

    public TDGGlobals(double[] parameters) {
        this.tau = parameters[0];
        this.kappa = parameters[1];

        double[] pi_ = CoreUtils.alr_inv(new double[]{parameters[2], parameters[3], parameters[4]});
        this.pi = new double[]{pi_[0], pi_[1], pi_[2], 1 - (pi_[0] + pi_[1] + pi_[2])};

        this.mu = parameters[5];
        this.gamma = 0;
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
        return "TDGGlobals{ " +
                "-tau " + tau +
                " -kappa " + kappa +
                " -pi " + Doubles.join(",", pi) +
                " -mu " + mu +
                " (1/nu=" + nu + ")" +
                '}';
    }

    public double[] getOptimiserParameters() {
        double[] parameters = new double[6]; // tau, kappa, pi (3), mu

        parameters[0] = this.tau;
        parameters[1] = this.kappa;

        double[] pi_inv = CoreUtils.alr(new double[]{this.pi[0], this.pi[1], this.pi[2]});

        parameters[2] = pi_inv[0];
        parameters[3] = pi_inv[1];
        parameters[4] = pi_inv[2];

        parameters[5] = this.mu;

        return parameters;
    }

    public static void main(String[] args) {
        TDGGlobals g = new TDGGlobals(0.18384062387520766, 4.128219745551803, new double[]{0.28119571640365953, 0.24370918125910898, 0.3741515588418039, 0.10094354349542756}, 2.6629836091051446, 0);
        System.out.printf("in: %s\n", g.toString());
        System.out.printf("out: %s\n", Doubles.join(", ", g.getOptimiserParameters()));
        System.out.printf("in2: %s\n", (new TDGGlobals(g.getOptimiserParameters())).toString());
    }

}
