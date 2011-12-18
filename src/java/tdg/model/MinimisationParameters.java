package tdg.model;

import com.google.common.primitives.Doubles;

/**
 * @author Asif Tamuri (atamuri@nimr.mrc.ac.uk)
 */
public class MinimisationParameters {
    private double[] parameters;
    private double[] stepSize;
    private double[] lowerBounds;
    private double[] upperBounds;

    public MinimisationParameters(double[] parameters, double[] stepsize, double[] lowerBounds, double[] upperBounds) {
        this.parameters = parameters;
        this.stepSize = stepsize;
        this.lowerBounds = lowerBounds;
        this.upperBounds = upperBounds;
    }

    public double[] getParameters() {
        return parameters;
    }

    public double[] getStepSize() {
        return stepSize;
    }

    public double getLowerBounds(int i) {
        return lowerBounds[i];
    }

    public double getUpperBounds(int i) {
        return upperBounds[i];
    }

    @Override
    public String toString() {
        return "model.ParametersForMinimisation{" +
                "parameters=" + (parameters == null ? null : Doubles.join(", ", parameters)) +
                ", stepsize=" + (stepSize == null ? null : Doubles.join(", ", stepSize)) +
                ", lowerBounds=" + (lowerBounds == null ? null : Doubles.join(", ", lowerBounds)) +
                ", upperBounds=" + (upperBounds == null ? null : Doubles.join(", ", upperBounds)) +
                '}';
    }
}
