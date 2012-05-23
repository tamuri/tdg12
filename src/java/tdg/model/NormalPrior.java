package tdg.model;


public class NormalPrior implements Prior {

    private final double sigma;

    public NormalPrior(double sigma) {
        this.sigma = sigma;
    }

    @Override
    public double calculate(double[] parameters) {

        double coeff = 1 / (2 * Math.pow(sigma, 2));

        double sum = 0;
        for (double f : parameters) {
            sum += Math.pow(f, 2);
        }

        return -(coeff * sum);
    }
}
