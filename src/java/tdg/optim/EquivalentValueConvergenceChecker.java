package tdg.optim;

import org.apache.commons.math.optimization.RealConvergenceChecker;
import org.apache.commons.math.optimization.RealPointValuePair;
import org.apache.commons.math.optimization.SimpleScalarValueChecker;

/**
 * Same as the standard scalar value convergence checker, but regards function has converged if
 * 50 consecutive evaluations do not change by more than 1E-6.
 *
 * @author Asif Tamuri (atamuri@nimr.mrc.ac.uk)
*/
public class EquivalentValueConvergenceChecker implements RealConvergenceChecker {
    RealConvergenceChecker convergenceChecker = new SimpleScalarValueChecker();
    int sameValueIterations = 0;
    final int maxEquivalentIterations;
    final double epsilon;

    public EquivalentValueConvergenceChecker(double epsilon, int maxEquivalentIterations) {
        this.epsilon = epsilon;
        this.maxEquivalentIterations = maxEquivalentIterations;
    }

    @Override
    public boolean converged(int iteration, RealPointValuePair previous, RealPointValuePair current) {

        if (sameValueIterations > maxEquivalentIterations) return true;

        // we've violated a constraint
        if (previous.getValue() > 9990 || current.getValue() > 9990) return false;

        // save the number of iterations where the function evaluation doesn't change
        if (Math.abs(previous.getValue() - current.getValue()) < epsilon) {
            sameValueIterations++;
        } else {
            sameValueIterations = 0;
        }

        return convergenceChecker.converged(iteration, previous, current);
    }
}
