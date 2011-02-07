package tdg.models.parameters;

import com.google.common.primitives.Doubles;

/**
 * @author Asif Tamuri
 * @version $Id: Fitness.java 128 2010-08-10 12:19:12Z tamuri $
 */
public class Fitness extends Parameter {
    public Fitness(double[] fitness, boolean optimise) {
        set(fitness);
        setOptimiseValue(optimise);
    }

    public double[] get() {
        return (double[]) super.get();
    }

    public String toString() {
        return "Fitness = [" + Doubles.join(", ", get()) + "]";
    }

}
