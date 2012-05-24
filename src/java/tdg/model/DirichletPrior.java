package tdg.model;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;

public class DirichletPrior implements Prior {
    private static final int DIM = 19;
    private Algebra algebra = new Algebra();
    private DoubleMatrix2D jacobi = DoubleFactory2D.dense.make(DIM, DIM);
    private final double alpha;

    public DirichletPrior(double alpha) {
        this.alpha = alpha;
    }

    @Override
    public double calculate(double[] fitness) {
        double[] theta = logit(fitness);
        double det = algebra.det(jacobian(theta));
        double density = density(theta);

        return Math.log(det * density);
    }

    private double[] logit(double[] fitness) {
        double[] theta = new double[fitness.length];

        double sumExpF = 0.0;

        for (double f : fitness)
            sumExpF += Math.exp(f);

        for (int i = 0; i < fitness.length; i++)
            theta[i] = Math.exp(fitness[i]) / (1 + sumExpF);

        return theta;
    }

    private DoubleMatrix2D jacobian(double[] theta) {
        for (int i = 0; i < DIM; i++) {
            for (int j = 0; j < DIM; j++) {
                if (i == j) {
                    jacobi.setQuick(i, j, theta[i] * (1 - theta[i]));
                } else {
                    jacobi.setQuick(i, j, -theta[i] * theta[j]);
                }
            }
        }
        return jacobi;
    }

    public double density(double[] theta) {
        double density = 1.0;
        double sum = 0.0;

        for (double x : theta) {
            density = density * Math.pow(x, alpha - 1);
            sum += x;
        }

        density = density * Math.pow(1 - sum, alpha - 1);

        return density;
    }
}
