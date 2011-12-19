package tdg.results;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import tdg.cli.ApproxOption;
import tdg.cli.GeneticCodeOption;
import tdg.cli.GlobalsOptions;

/**
 * @author Asif Tamuri (atamuri@nimr.mrc.ac.uk)
 */
class Options {
    @Parameter(names = "-o", description = "The results file from tdg.Analyse", required = true)
    public String outputFile;

    @ParametersDelegate
    public ApproxOption approx = new ApproxOption();

    @ParametersDelegate
    public GlobalsOptions globals = new GlobalsOptions();

    @ParametersDelegate
    public GeneticCodeOption gc = new GeneticCodeOption();
}
