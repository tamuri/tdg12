package tdg;

/**
 * Various constants used throughout the program but not expected to be changed by end-users.
 *
 * @author Asif Tamuri (atamuri@nimr.mrc.ac.uk)
 * @version 1.0
 */
public class Constants {
    // Random initial value for fitness parameter from range -INITIAL_PARAM_RANGE to INITIAL_PARAM_RANGE
    public static final int INITIAL_PARAM_RANGE = 3;
    public static final double CONVERGENCE_TOL = 1E-6;
    public static final double FITNESS_BOUND = 20; // effectively -20 = -Infinity and 20 = Infinity for fitness parameter
    public static final int MAX_EVALUATIONS = 10000;
    public static final int VERY_BAD_LIKELIHOOD = 9999999;
    public static final double[] CLADE_BRANCH_SPLIT = {0.5, 0.5};
    public static final double SCALING_THRESHOLD = 1e-15;
    public static final int SCALING_NODE_STEP = 5;
    public static final double FITNESS_INITIAL_VALUE = 0.0;
    
    public static final String S_FILENAME = "S.txt";
    public static final String Q0_FILENAME = "Q0.txt";
    public static final String QS_FILENAME = "QS.txt";
    public static final String PI_FILENAME = "PiS.txt";
    public static final String PIAA_FILENAME = "PiAA.txt";
    public static final String F_FILENAME = "F.txt";

}
