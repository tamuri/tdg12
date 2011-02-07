package tdg.utils;

import cern.colt.matrix.DoubleMatrix2D;

/**
 * @author Asif Tamuri
 * @version $Id: ColtMatrixUtils.java 149 2010-08-14 18:48:07Z tamuri $
 */
public class ColtMatrixUtils {
    // From c++ source:
    // http://code.google.com/p/numerical-recipes-java/source/browse/trunk/numerical-recipes-j/core/src/main/cpp/eigen_unsym.h
    public static void balance(DoubleMatrix2D a) {
        double RADIX = 2;
        boolean done = false;
        double RADIX_SQ = RADIX * RADIX;

        int n = a.rows();

        while (!done) {
            done = true;
            for (int i = 0; i < n; i++) {
                double r = 0.0, c = 0.0;
                for (int j = 0; j < n; j++) {
                    if (j != i) {
                        c += Math.abs(a.getQuick(j, i));
                        r += Math.abs(a.getQuick(i, j));
                    }
                }
                if (c != 0.0 && r != 0.0) {
                    double g = r / RADIX;
                    double f = 1.0;
                    double s = c + r;
                    while (c < g) {
                        f *= RADIX;
                        c *= RADIX_SQ;
                    }
                    g = r * RADIX;
                    while (c > g) {
                        f /= RADIX;
                        c /= RADIX_SQ;
                    }
                    if ((c + r) / f < 0.95 * s) {
                        done = false;
                        g = 1.0 / f;
                        // scale[i] *= f;
                        for (int j = 0; j < n; j++)
                            a.setQuick(i, j, a.getQuick(i, j) * g);
                        for (int j = 0; j < n; j++)
                            a.setQuick(j, i, a.getQuick(j, i) * f);
                    }
                }
            }
        }

    }
}
