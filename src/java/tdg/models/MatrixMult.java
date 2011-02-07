package tdg.models;

import cern.colt.matrix.DoubleMatrix1D;

/**
 * @author Asif Tamuri
 * @version $Id$
 */
public interface MatrixMult {
    public void getProbabilityMatrix(final double[] matrix, final double branchLength, final DoubleMatrix1D lambda, final double[] U, final double[] UInv, final int[] siteCodons);
}
