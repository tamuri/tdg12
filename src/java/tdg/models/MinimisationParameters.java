package tdg.models;

import com.google.common.primitives.Doubles;

/**
 * @author Asif Tamuri
 * @version $Id: MinimisationParameters.java 25 2010-04-19 14:39:24Z atamuri $
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
        return "models.ParametersForMinimisation{" +
                "parameters=" + (parameters == null ? null : Doubles.join(", ", parameters)) +
                ", stepsize=" + (stepSize == null ? null : Doubles.join(", ", stepSize)) +
                ", lowerBounds=" + (lowerBounds == null ? null : Doubles.join(", ", lowerBounds)) +
                ", upperBounds=" + (upperBounds == null ? null : Doubles.join(", ", upperBounds)) +
                '}';
    }
}
