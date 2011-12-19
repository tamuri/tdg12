package tdg.model;

import com.google.common.primitives.Doubles;

/**
 * The LikelihoodCalculator optimises one or more parameters. Fitness is one of the parameters that can be optimised (in
 * fact, the only parameter in the swMutSel0 model). It holds an array of fitness values for each amino acid.
 *
 * @author Asif Tamuri (atamuri@nimr.mrc.ac.uk)
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
