package tdg.optim;

import org.apache.commons.math.FunctionEvaluationException;
import org.apache.commons.math.analysis.MultivariateRealFunction;
import tdg.Constants;
import tdg.models.LikelihoodCalculator;

/**
 * Wraps the likelihood function in a MultivariateRealFunction, as needed by the Apache Commons Math optimisation
 * routines.
 *
 * @author Asif Tamuri
 * @version $Id: LikelihoodMaximiser.java 152 2010-11-08 11:10:01Z tamuri $
 * @see LikelihoodCalculator
 */
public class LikelihoodMaximiser implements MultivariateRealFunction{
    private static final double CONSTRAINT = Constants.FITNESS_BOUND + 1;
    private LikelihoodCalculator lc;
    // private Map<DoubleArrayKey, Double> cache = Maps.newHashMap();

    public void setLc(LikelihoodCalculator lc) {
        this.lc = lc;
    }

    @Override
    public double value(double[] point) throws FunctionEvaluationException, IllegalArgumentException {
        // TODO: Theoretically, we should be using something like Powell's COBYLA with proper handling of constraints
        // Keep fitness parameters within our constraints - if any point is outside constraint, return a very bad likelihood value
        for (double d : point) {
            if (d < -CONSTRAINT || d > CONSTRAINT) {
                return Constants.VERY_BAD_LIKELIHOOD;
            }
        }   

        /*
        // If we want to cache results (some optimisers seem to cover the same ground)
        DoubleArrayKey dak = new DoubleArrayKey(point);
        if (cache.containsKey(dak)) {
            return cache.get(dak);
        }
        cache.put(dak, out);
        */

        return lc.function(point);
    }
}
