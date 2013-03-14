package tdg;

/**
 * Various constants used throughout the program but not expected to be changed by end-users.
 *
 * @author Asif Tamuri (atamuri@nimr.mrc.ac.uk)
 */
public class Constants {
    /**
     * One of the residues (usually the most observed) is fixed to this value; other fitnesses are relative to this one.
     */
    public static final double FITNESS_FIXED_FOR_RELATIVE = 0.0;

    /**
     * Random initial value for fitness parameter from range -RANDOM_INITIAL_FITNESS_RANGE to RANDOM_INITIAL_FITNESS_RANGE.
     */
    public static final int RANDOM_INITIAL_FITNESS_RANGE = 3;

    /**
     * When to assume convergence of MLE of fitness parameters (log-likelihood).
     */
    public static final double CONVERGENCE_TOL = 5e-6;

    /**
     * Effectively -20 = -Infinity and 20 = Infinity for fitness parameter.
     */
    public static final double FITNESS_BOUND = 20;

    /**
     * Terminate the optimisation routine after this number of evaluations.
     */
    public static final int MAX_EVALUATIONS = 10000;

    /**
     * To constrain MLE of fitness within FITNESS_BOUND, return this log-likelihood if we step outside of bounds.
     */
    public static final double VERY_BAD_LIKELIHOOD = Double.NEGATIVE_INFINITY;

    /**
     * Yang's method of scaling partials on very large trees to avoid numerical underflow.
     */
    public static final double SCALING_THRESHOLD = 1e-5;
    public static final int SCALING_NODE_STEP = 5;

    /**
     * How to split the branch connecting heterogeneous models (e.g. 0.5 = half-way).
     */
    public static final double CLADE_BRANCH_SPLIT = 0.5;

    /**
     * The filenames for parsed results files.
     */
    public static final String S_FILENAME = "S.txt";
    public static final String Q0_FILENAME = "Q0.txt";
    public static final String QS_FILENAME = "QS.txt";
    public static final String PI_FILENAME = "PiS.txt";
    public static final String PIAA_FILENAME = "PiAA.txt";
    public static final String F_FILENAME = "F.txt";

    /**
     * The name of the program .jar (for command-line option help)
     */
    public static final String PROGRAM_JAR = "tdg12.jar";
    public static final double INITIAL_BRANCH_LENGTH = 0.1;
}
