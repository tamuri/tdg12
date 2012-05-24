package tdg.cli;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.converters.DoubleConverter;

/**
 * @author Asif Tamuri (atamuri@nimr.mrc.ac.uk)
 */
public class GlobalsOptions {
    @Parameter(names = "-tau", description = "Rate of multiple substitutions.", converter = DoubleConverter.class, required = true)
    public double tau;

    @Parameter(names = "-kappa", description = "Transition/transversion bias.", converter = DoubleConverter.class, required = true)
    public double kappa;

    @Parameter(names = "-pi", description = "Comma-separated base nucleotide frequencies (T,C,A,G).", converter = DoubleArrayConverter.class, required = true)
    public double[] pi;

    @Parameter(names = "-mu", description = "Branch/rate scaling factor.", converter = DoubleConverter.class, required = true)
    public double mu;

    @Parameter(names = "-gamma", description = "Probability of error.", converter = DoubleConverter.class, hidden = true)
    public double gamma = 0;

}
