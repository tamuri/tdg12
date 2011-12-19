package tdg.results;

import com.beust.jcommander.JCommander;
import tdg.Constants;
import tdg.model.TDGGlobals;

/**
 * Runs all the various procedures for creating the distribution of selection coefficients from raw results files
 *
 * @author Asif Tamuri (atamuri@nimr.mrc.ac.uk)
 */
public class All {
    public static void main(String[] args) throws Exception {
        // Parse command-line options
        Options options = new Options();
        JCommander jc = new JCommander(options);
        jc.parse(args);

        // Write fitnesses to F.txt given the raw results file
        FitnessExtractor fe = new FitnessExtractor();
        fe.extract(options.outputFile);

        // Write the model parameter values given the fitnesses and global parameters
        ModelWriter mw = new ModelWriter(options,
                new TDGGlobals(options.globals.tau, options.globals.kappa, options.globals.pi, options.globals.mu, options.globals.gamma),
                Constants.F_FILENAME);
        mw.run();

        // Use the model parameter values to create the distribution of selection coefficients
        DistributionWriter dw = new DistributionWriter();
        dw.run();
    }
}
