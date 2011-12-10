package tdg;

/**
 * Various constants used throughout the program but not expected to be changed by end-users.
 * @author Asif Tamuri (atamuri@nimr.mrc.ac.uk)
 */
public class Constants {
    // Random initial value for fitness parameter from range -INITIAL_PARAM_RANGE to INITIAL_PARAM_RANGE
    public static final int INITIAL_PARAM_RANGE = 3;
    public static final double CONVERGENCE_TOL = 1E-6;
    public static final double FITNESS_BOUND = 20; // effectively -20 = -Infinity and 20 = Infinity for fitness parameter
}
