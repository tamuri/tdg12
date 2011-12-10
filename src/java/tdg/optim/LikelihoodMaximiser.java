package tdg.optim;

import com.google.common.primitives.Doubles;
import org.apache.commons.math.FunctionEvaluationException;
import org.apache.commons.math.analysis.MultivariateRealFunction;
import tdg.models.LikelihoodCalculator;

/**
 * @author Asif Tamuri
 * @version $Id: LikelihoodMaximiser.java 152 2010-11-08 11:10:01Z tamuri $
 */
public class LikelihoodMaximiser implements MultivariateRealFunction{
    private static final int CONSTRAINT = 21;
    private int count = 0;
    private LikelihoodCalculator lc;
    // private Map<DoubleArrayKey, Double> cache = Maps.newHashMap();

    public void setLc(LikelihoodCalculator lc) {
        this.lc = lc;
    }

    @Override
    public double value(double[] point) throws FunctionEvaluationException, IllegalArgumentException {
        // don't give very small or very large fitness coefficients
        for (double d : point) {
            if (d < -CONSTRAINT || d > CONSTRAINT) {
                return 999999;
            }
        }   

        /*DoubleArrayKey dak = new DoubleArrayKey(point);
        if (done.containsKey(dak)) {
            System.out.printf("HIT in site params!\n");
            return done.get(dak);
        }*/

        double out = lc.function(point);
        //System.out.printf("Params= { %s }; \n ", Doubles.join(", ", point));
        //System.out.printf("%s - %s \n", count++, out);
        //count++;
        //if (count % 500 == 0) System.out.printf("%s ", count);
        // done.put(dak, out);

        //System.exit(0);
        return out;
    }
}
