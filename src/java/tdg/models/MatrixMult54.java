package tdg.models;

import cern.colt.matrix.DoubleMatrix1D;
import tdg.utils.CodeTimer;
import tdg.utils.GeneticCode;

/**
 * @author Asif Tamuri
 * @version $Id$
 */
public class MatrixMult54 implements MatrixMult {
    private static final int size = 54;
    
    private final double[] PtTemp = new double[size * size];

    @Override
    public void getProbabilityMatrix(final double[] matrix, final double branchLength, final DoubleMatrix1D lambda, final double[] U, final double[] UInv, final int[] siteCodons) {
        long start = CodeTimer.start();
        for (int i = 0; i < size; i++) {
            double temp = Math.exp(branchLength * lambda.getQuick(i));
            for (int j = 0; j < size; j++) {
                PtTemp[j * size + i] = temp * U[j * size + i];
            }
        }
        CodeTimer.store("getProbabilityMatrix_1", start);

        long start2 = CodeTimer.start();
        for (int j = 0; j < size; j++) {
            for (int i = 0; i < size; i++) {
                double temp = 0;
                for (int k = 0; k < size; k++) {
                    temp += PtTemp[i * size + k] * UInv[k * size + j];
                }
                if (temp < 0) temp = 0;
                matrix[siteCodons[i] * GeneticCode.CODON_STATES + siteCodons[j]] = temp;
            }
        }
        CodeTimer.store("getProbabilityMatrix_2", start2);
    }
}
