package tdg.model;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;

public class DirichletPrior implements Prior {
    private Algebra algebra;
    private DoubleMatrix2D jacobi = DoubleFactory2D.dense.make(19, 19);
    private final double alpha;

    public DirichletPrior(double alpha) {
        this.alpha = alpha;
        this.algebra = new Algebra();
    }

    public static double quickDensity(double[] xx, double alpha) {
        double density = 1.0;
        double sum = 0.0;

        for (double x : xx) {
            density = density * Math.pow(x, alpha - 1);
            sum += x;
        }

        density = density * Math.pow(1 - sum, alpha - 1);

        return  density;
    }

    @Override
    public double calculate(double[] parameters) {
        double[] invLogit = getInverseLogit(parameters);
        double det = algebra.det(makeJacobian(invLogit));

        //Arrays.fill(alpha, 1.0);
        //double density = DirichletDist.density(alpha, invLogit);
        double density = quickDensity(invLogit, alpha);
        //System.out.printf("det = %s * den = %s\n", det, density);
        return Math.log(det * density);
    }

    private double[] getInverseLogit(double[] fitnesses) {
        double[] invLogit = new double[fitnesses.length];

        double sumExpF = 0.0;
        for (double f : fitnesses) {
            sumExpF += Math.exp(f);
        }
        //sumExpF += Math.exp(0); // = 1

        for (int i = 0; i < fitnesses.length; i++) {
            invLogit[i] = Math.exp(fitnesses[i]) / (1 + sumExpF);
        }

        //invLogit[invLogit.length - 1] = Math.exp(0) / (1 + sumExpF);

        return invLogit;
    }

    private DoubleMatrix2D makeJacobian(double[] theta) {
        for (int i = 0; i < 19; i++) {
            for (int j = 0; j < 19; j++) {
                if (i == j) {
                    jacobi.setQuick(i, j, theta[i] * (1 - theta[i]));
                } else {
                    jacobi.setQuick(i, j, -theta[i] * theta[j]);
                }
            }
        }

        return jacobi;
    }


}
