package tdg.optim;

import org.apache.commons.math.optimization.RealConvergenceChecker;
import org.apache.commons.math.optimization.RealPointValuePair;
import org.apache.commons.math.optimization.SimpleRealPointChecker;
import org.apache.commons.math.optimization.SimpleScalarValueChecker;

/**
* User: atamuri
* Date: 08/02/11
* Time: 15:08
*/
public class PointAndValueConvergenceChecker implements RealConvergenceChecker {
    RealConvergenceChecker rcc1;
    RealConvergenceChecker rcc2;
    int sameValueIterations = 0;
    final int maxEquivalentIterations;
    final double equivalenceEpsilon;

    public PointAndValueConvergenceChecker(double epsilon, int maxEquivalentIterations, double equivalenceEpsilon) {
        rcc1 = new SimpleRealPointChecker(epsilon, epsilon);
        rcc2 = new SimpleScalarValueChecker(epsilon, epsilon);
        this.maxEquivalentIterations = maxEquivalentIterations;
        this.equivalenceEpsilon = equivalenceEpsilon;
    }

    @Override
    public boolean converged(int iteration, RealPointValuePair previous, RealPointValuePair current) {

        if (sameValueIterations > maxEquivalentIterations) return true;

        // we've violated a constraint
        if (previous.getValue() > 9990 || current.getValue() > 9990) return false;

        // save the number of iterations where the function evaluation doesn't change
        if (Math.abs(previous.getValue() - current.getValue()) < equivalenceEpsilon) {
            sameValueIterations++;
        } else {
            sameValueIterations = 0;
        }

        return rcc1.converged(iteration, previous, current) && rcc2.converged(iteration, previous, current);
    }
}
