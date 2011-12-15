package tdg.results;

import com.beust.jcommander.Parameter;
import tdg.cli.DoubleArrayConverter;
import tdg.cli.DoubleConverter;
import tdg.cli.GeneticCodeConverter;
import tdg.utils.GeneticCode;

/**
 * @author Asif Tamuri (atamuri@nimr.mrc.ac.uk)
 */
class Options {
    @Parameter(names = "-tau", description = "Rate of multiple substitutions.", converter = DoubleConverter.class, required = true)
    public double tau;

    @Parameter(names = "-kappa", description = "Transition/transversion bias.", converter = DoubleConverter.class, required = true)
    public double kappa;

    @Parameter(names = "-pi", description = "Comma-separated base nucleotide frequencies (T,C,A,G).", converter = DoubleArrayConverter.class, required = true)
    public double[] pi;

    @Parameter(names = "-mu", description = "Branch/rate scaling factor.", converter = DoubleConverter.class, required = true)
    public double mu;

    double gamma = 0;

    @Parameter(names = "-gc", description = "The genetic code translation to use (standard or vertebrate_mit).", required = true, converter = GeneticCodeConverter.class)
    public GeneticCode geneticCode;

    @Parameter(names = "-approx", description = "Use the approximate method to optimise the likelihood")
    public boolean approx = false;

}
