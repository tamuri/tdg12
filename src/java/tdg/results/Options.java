package tdg.results;

import com.beust.jcommander.ParametersDelegate;
import tdg.cli.ApproxOption;
import tdg.cli.GeneticCodeOption;
import tdg.cli.GlobalsOptions;

/**
 * @author Asif Tamuri (atamuri@nimr.mrc.ac.uk)
 */
class Options {
    @ParametersDelegate
    public ApproxOption approx = new ApproxOption();

    @ParametersDelegate
    public GlobalsOptions globals = new GlobalsOptions();

    @ParametersDelegate
    public GeneticCodeOption gc = new GeneticCodeOption();
}
