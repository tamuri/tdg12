package tdg.model;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.EigenvalueDecomposition;
import cern.jet.math.Functions;
import com.google.common.primitives.Ints;
import tdg.Constants;
import tdg.utils.GeneticCode;
import tdg.utils.PhyloUtils;

import java.util.List;

/**
 * TDG10 codon selection model.
 * <p/>
 * The full formulation is:
 * <p/>
 * \f[
 * q_{IJ}  = \left[ {
 * \mu \tau ^{n - 1} \kappa ^{n_t }
 * \frac{1}{{\prod\limits_{k \ne k\prime } {\pi _{j_{k\prime } }^* } }}
 * \times
 * \frac{{F_J  - F_I }}{{{\rm{e}}^{F_J }  - {\rm{e}}^{F_I } }}
 * } \right]
 * \times
 * \left( {\prod\limits_{k = 1}^3 {\pi _{j_k }^* }} \right)
 * {\rm{e}}^{F_J}
 * \f]
 *
 * @author Asif Tamuri (atamuri@nimr.mrc.ac.uk)
 */
public class TDGCodonModel {
    private final TDGGlobals globals;
    private final Fitness fitness;
    private final int[] siteCodons;
    private final int matrixSize;
    private final double[] codonPi;
    private final int[] aminoAcidsToFitness;

    private final double[] Q; // substitution matrix
    private final DoubleMatrix2D B; // PI^0.5 * Q * PI^-0.5

    private DoubleMatrix1D lambda; // eigenvalue of B
    // For calculating the PT matrix
    private final double[] U;
    private final double[] UInv;
    private final double[] PtTemp;

    // Map<Double, double[]> probMatrixStore = Maps.newHashMap();

    public TDGCodonModel(TDGGlobals globals, Fitness fitness, List<Integer> aminoAcids) {
        this.globals = globals;
        this.fitness = fitness;

        this.siteCodons = Ints.toArray(PhyloUtils.getCodonsFromAminoAcids(aminoAcids));

        this.aminoAcidsToFitness = new int[GeneticCode.AMINO_ACID_STATES];

        for (int i = 0; i < GeneticCode.AMINO_ACID_STATES; i++) {
            if (aminoAcids.contains(i)) {
                aminoAcidsToFitness[i] = aminoAcids.indexOf(i);
            } else {
                aminoAcidsToFitness[i] = -1;
            }
        }

        this.matrixSize = siteCodons.length;

        this.codonPi = new double[matrixSize];
        this.Q = new double[matrixSize * matrixSize];
        this.B = DoubleFactory2D.dense.make(matrixSize, matrixSize);
        this.U = new double[matrixSize * matrixSize];
        this.UInv = new double[matrixSize * matrixSize];
        this.PtTemp = new double[matrixSize * matrixSize];
    }

    public void updateModel() {
        makePI();
        makeQ();
        makeB();
        doEigenValueDecomposition();
        // probMatrixStore.clear();
    }

    /**
     * Performs 2 tasks:
     * <p/>
     * 1. Calculates \f[ \pi _J  = \left( {\prod\limits_{k = 1}^3 {\pi _{j_k }^* } } \right){\rm{e}}^{F_J } \f]
     * <p/>
     * 2. Creates diagonal matrix of codon frequencies, which can be used for
     * eigenvalue decomposition of Q via B = PI^0.5 Q PI^-0.5
     */
    private void makePI() {
        double z = 0;
        for (int i = 0; i < matrixSize; i++) {
            codonPi[i] = globals.getCodonPiProduct(siteCodons[i]) *
                    Math.exp(fitness.get()[aminoAcidsToFitness[GeneticCode.getInstance().getAminoAcidIndexFromCodonIndex(siteCodons[i])]]);
            z += codonPi[i];
        }

        for (int i = 0; i < matrixSize; i++) {
            codonPi[i] /= z;
        }
    }

    /**
     * Construct Q matrix directly using formulation:
     * <p/>
     * \f[
     * q_{IJ}  = \nu  \times \left( {\tau ^{n - 1} \kappa ^{n_t } \prod\limits_{k,i_k  \ne j_k } {\pi _{j_k }^* } } \right) \times h(S_{IJ} )
     * \f]
     * <p/>
     * without using separate S and PI matrices. The two methods are equivalent.
     */
    private void makeQ() {
        double[] f = fitness.get();

        for (int i = 0; i < matrixSize; i++) {
            for (int j = 0; j < matrixSize; j++) {

                if (i == j) {
                    Q[i * matrixSize + j] = 0;
                    continue;
                }

                int cI = siteCodons[i];
                int cJ = siteCodons[j];

                double muIJ = globals.getNeutralMutationRate(cI, cJ);

                double fI = f[aminoAcidsToFitness[GeneticCode.getInstance().getAminoAcidIndexFromCodonIndex(cI)]];
                double fJ = f[aminoAcidsToFitness[GeneticCode.getInstance().getAminoAcidIndexFromCodonIndex(cJ)]];

                double hS = getRelativeFixationProbability(fJ - fI);

                Q[i * matrixSize + j] = globals.getNu() * muIJ * hS;
            }
        }

        for (int row = 0; row < matrixSize; row++) {
            double sum = 0;
            for (int j = 0; j < matrixSize; j++) {
                sum += Q[row * matrixSize + j];
            }
            Q[row * matrixSize + row] = -sum;
        }
    }

    private double getRelativeFixationProbability(double Sij) {
        if (Sij == 0) return 1;
        if (Sij < -1e3) return 0;
        if (Sij > 1e3) return Sij;
        return Sij / -Math.expm1(-Sij);
    }

    private void makeB() {
        for (int i = 0; i < matrixSize; i++) {
            for (int j = 0; j < matrixSize; j++) {
                B.setQuick(i, j, Q[i * matrixSize + j] * Math.sqrt(codonPi[i]) / Math.sqrt(codonPi[j]));
            }
        }
    }

    /**
     * Performs eigenvalue decomposition on Q. Two methods:
     * <p/>
     * 1. Diagonalise Q directly.
     * <p/>
     * 2. Gets eigenvalues and eigenvectors from a constructed symmetric matrix,
     * B, based on PI and Q. Diagonalising B is faster than Q because B is
     * symmetric.
     * <p/>
     * Described in Yang 2006, Computational Molecular Evolution, pp 68-9.
     * <p/>
     * To start,
     * \f[ B = \Pi ^{{1 \mathord{\left/ {\vphantom {1 2}} \right. \kern-\nulldelimiterspace} 2}}
     * Q
     * \Pi ^{ - {1 \mathord{\left/ {\vphantom {1 2}} \right. \kern-\nulldelimiterspace} 2}}
     * \f]
     * <p/>
     * then diagonalise B,
     * \f[ B = U\Lambda U^{ - 1} \f]
     * <p/>
     * where U is the U and ∆ is the lambda on the diagonal.
     * <p/>
     * then,
     * <p/>
     * \f[ Q = (\Pi ^{ - {1 \mathord{\left/ {\vphantom {1 2}} \right. \kern-\nulldelimiterspace} 2}} U)
     * \Lambda
     * (\Pi ^{{1 \mathord{\left/ {\vphantom {1 2}} \right. \kern-\nulldelimiterspace} 2}} U^{ - 1} )
     * \f]
     */
    private void doEigenValueDecomposition() {
        EigenvalueDecomposition evdB = new EigenvalueDecomposition(B);

        lambda = DoubleFactory2D.dense.diagonal(evdB.getD());
        // we scale branch length by global parameter mu here so we only have do it once
        lambda.assign(Functions.mult(globals.getMu()));

        DoubleMatrix2D R = evdB.getV();

        for (int i = 0; i < matrixSize; i++) {
            double piSqrt = Math.sqrt(codonPi[i]);
            double piInvSqrt = 1 / piSqrt;
            for (int j = 0; j < matrixSize; j++) {
                U[i * matrixSize + j] = piInvSqrt * R.getQuick(i, j);
                UInv[j * matrixSize + i] = piSqrt * R.getQuick(i, j); // inverse(R) == transpose(R)
            }
        }

        /* // Eigenvalue decomposition of non-symmetric Q
        EigenvalueDecomposition evdQ = new EigenvalueDecomposition(Q);
        lambda = DoubleFactory2D.dense.diagonal(evdQ.getD());
        lambda.assign(Functions.mult(globals.getMu()));
        U = evdQ.getV();
        UInv = algebra.inverse(evdQ.getV());
        */
    }


    public void getProbabilityMatrix(final double[] matrix, final double branchLength) {

        /*if (probMatrixStore.containsKey(branchLength)) {
            double[] ref = probMatrixStore.get(branchLength);
            for (int i = 0; i < matrix.length; i++) {
                matrix[i] = ref[i];
            }
        } else {*/
        //long start = CodeTimer.start();
        // NOTE: If matrixSize were a final static int, then there would be some
        // performance improvement. Java doesn't array bounds check
        // in that case. See http://wikis.sun.com/display/HotSpotInternals/RangeCheckElimination
        // MitoData, site 274, time in this method ~56secs down to ~50secs (54x54 codon matrix)
        // More improvement using jdk1.7.0 ~58secs to ~42secs!!
        for (int i = 0; i < matrixSize; i++) {
            double temp = Math.exp(branchLength * lambda.getQuick(i));
            for (int j = 0; j < matrixSize; j++) {
                PtTemp[j * matrixSize + i] = temp * U[j * matrixSize + i];
            }
        }


        /*
        * instead of multiplying by U[j,i] here and then UInv[k,j] below, we can calculate the product
        * when we do the eigenvalue decomposition above i.e.:
        * for (i <- 1..64; j <- 1..64; k <- 1..64) { UUInv[x++] = U[i][k] * UInv[k][j] }
        */

        //CodeTimer.store("getProbabilityMatrix_1", start);

        //long start2 = CodeTimer.start();x
        for (int j = 0; j < matrixSize; j++) {
            for (int i = 0; i < matrixSize; i++) {
                double temp = 0;
                for (int k = 0; k < matrixSize; k++) {
                    temp += PtTemp[i * matrixSize + k] * UInv[k * matrixSize + j];
                }
                if (temp < 0) temp = 0;
                matrix[siteCodons[i] * GeneticCode.CODON_STATES + siteCodons[j]] = temp;
            }
        }

/*
            probMatrixStore.put(branchLength, Arrays.copyOf(matrix, GeneticCode.CODON_STATES * GeneticCode.CODON_STATES));

            //CodeTimer.store("getProbabilityMatrix_2", start2);
        }*/
// */
        /* This works, but would be horrible to have all these classes though...wouldn't it!?

        Briefly, if we had a class for each size of matrix (from some lower-limit to 64),
        each class having a compile-time static scale for matrix size, then the JVM will
        eliminate range checks. As shown above, this could lead to a "substantial" (hmm..20s)
        performance improvement

        As class member:
            MatrixMult matrixMultiplier;
        In constructor:
            try {
                this.matrixMultiplier = (MatrixMult) Class.forName("tdg.model.MatrixMult" + matrixSize).newInstance();
            } catch (Exception e) { e.printStackTrace(); }
        Then here:
            matrixMultiplier.getProbabilityMatrix(matrix, branchLength, lambda, U, UInv, siteCodons);

        NOTE: this probably only needs to be a pure matrix multiplier by performing
        the quicker lambda multiplication in this TDGCodonModel method
        */
    }

    /**
     * @return array of equilibrium frequencies for all 64 codons
     */
    public double[] getCodonFrequencies() {
        double[] fullF = new double[GeneticCode.CODON_STATES];
        for (int i = 0; i < matrixSize; i++) {
            fullF[siteCodons[i]] = codonPi[i];
        }
        return fullF;
    }

    public double[] getAminoAcidFrequencies() {
        double[] freqs = new double[GeneticCode.AMINO_ACID_STATES];
        double[] codonPis = getCodonFrequencies();
        for (int i = 0; i < GeneticCode.CODON_STATES; i++) {
            if (GeneticCode.getInstance().isUnknownAminoAcidState(GeneticCode.getInstance().getAminoAcidIndexFromCodonIndex(i)))
                continue;
            freqs[GeneticCode.getInstance().getAminoAcidIndexFromCodonIndex(i)] += codonPis[i];
        }
        return freqs;
    }

    public int[] getSiteCodons() {
        return siteCodons;
    }

    public double[] getFullQ() {
        double[] fullQ = new double[GeneticCode.CODON_STATES * GeneticCode.CODON_STATES];

        // This is the entire 60x60 matrix for Q
        for (int i = 0; i < matrixSize; i++) {
            for (int j = 0; j < matrixSize; j++) {
                fullQ[siteCodons[i] * GeneticCode.CODON_STATES + siteCodons[j]] = Q[i * matrixSize + j];
            }
        }

        // add the stop codons explicitly
        double stopS = getRelativeFixationProbability(Constants.FITNESS_BOUND - -Constants.FITNESS_BOUND);

        for (int i = 0; i < GeneticCode.CODON_STATES; i++) {
            for (int j = 0; j < GeneticCode.CODON_STATES; j++) {
                int aa_from = GeneticCode.getInstance().getAminoAcidIndexFromCodonIndex(i);
                int aa_to = GeneticCode.getInstance().getAminoAcidIndexFromCodonIndex(j);
                // if going to a stop codon
                if (aa_to < 0) {
                    // rate is 0
                    fullQ[i * GeneticCode.CODON_STATES + j] = 0;
                    // if coming from a stop codon and going to a non-stop codon
                } else if (aa_from < 0) {
                    fullQ[i * GeneticCode.CODON_STATES + j] = globals.getNeutralMutationRate(i, j) * globals.getNu() * stopS;
                }
            }
        }

        for (int row = 0; row < GeneticCode.CODON_STATES; row++) {
            double sum = 0;
            for (int j = 0; j < GeneticCode.CODON_STATES; j++) {
                sum += fullQ[row * GeneticCode.CODON_STATES + j];
            }
            fullQ[row * GeneticCode.CODON_STATES + row] = -sum;
        }

        return fullQ;
    }

    /**
     * Given the fitness parameters for each amino acid, and hence every codon (STOP codons have a "very bad fitness"),
     * this method constructs a 64*64 matrix of selection coefficients i.e. (Fj - Fi).
     *
     * @return a 64*64 matrix stored as 1D array of doubles
     */
    public double[] getS() {
        double[] fullS = new double[GeneticCode.CODON_STATES * GeneticCode.CODON_STATES];

        for (int i = 0; i < GeneticCode.CODON_STATES; i++) {
            for (int j = 0; j < GeneticCode.CODON_STATES; j++) {
                if (i == j) {
                    fullS[i * GeneticCode.CODON_STATES + j] = 0;
                } else {
                    int aa_from = GeneticCode.getInstance().getAminoAcidIndexFromCodonIndex(i);
                    int aa_to = GeneticCode.getInstance().getAminoAcidIndexFromCodonIndex(j);
                    if (aa_to < 0) {
                        // going to STOP codon
                        fullS[i * GeneticCode.CODON_STATES + j] = -Constants.FITNESS_BOUND;
                    } else if (aa_from < 0) {
                        // going from a STOP codon to any non-STOP codon
                        fullS[i * GeneticCode.CODON_STATES + j] = Constants.FITNESS_BOUND;
                    } else if (aminoAcidsToFitness[aa_from] > -1 && aminoAcidsToFitness[aa_to] > -1) {
                        // both amino acids that occur at this site
                        fullS[i * GeneticCode.CODON_STATES + j] = fitness.get()[aminoAcidsToFitness[aa_to]] - fitness.get()[aminoAcidsToFitness[aa_from]];
                    } else if (aminoAcidsToFitness[aa_from] > -1 && aminoAcidsToFitness[aa_to] == -1) {
                        // from observed to unobserved
                        fullS[i * GeneticCode.CODON_STATES + j] = -Constants.FITNESS_BOUND;
                    } else if (aminoAcidsToFitness[aa_from] == -1 && aminoAcidsToFitness[aa_to] > -1) {
                        // from unobserved to observed
                        fullS[i * GeneticCode.CODON_STATES + j] = Constants.FITNESS_BOUND;
                    } else {
                        // neither of these amino acids are observed at this site
                        fullS[i * GeneticCode.CODON_STATES + j] = -Constants.FITNESS_BOUND;
                    }
                }
            }
        }
        return fullS;
    }

}

