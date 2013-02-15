package tdg.model;


import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.jet.math.Functions;

public class NormalPrior implements Prior {

    private final double coefficient;

    public NormalPrior(double sigma) {
        this.coefficient = 1 / (2 * Math.pow(sigma, 2));
    }

    @Override
    public double calculate(double[] parameters) {
        DoubleMatrix1D f = DoubleFactory1D.dense.make(parameters);
        return -coefficient * f.aggregate(Functions.plus, Functions.pow(2.0));
    }

    @Override
    public String toString() {
        return "NormalPrior{coefficient 1/(2 * sigma^2)=" + coefficient + '}';
    }

    public static void main(String[] args) {
        for (double s : new double[]{0.1, 1, 2, 10, 100}) {
            NormalPrior n = new NormalPrior(s);

            double[] f = new double[]{-1.0135940331511981, 0.9254958707503568, -0.29255612654565816, -0.22681279848729385, 0.47387836055700755, -0.15291969247995885, -0.8747788652108226, -1.7177774513544954, 1.3341023647790993, -1.264047307128678, -3.129557842977072, 0.10073966555682189, 1.3617232021873107, 1.1981880441978687, 0.44439692707114636, -0.5022015366261262, 1.676426029781541, 0.020532038066831744, 0.008691116712494517};
            System.out.printf("%s\t%s\n", s, n.calculate(f));
        }
    }
}

