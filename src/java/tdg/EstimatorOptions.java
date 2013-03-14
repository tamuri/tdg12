package tdg;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.google.common.collect.Lists;
import tdg.cli.GeneticCodeOption;
import tdg.cli.PriorConverter;
import tdg.model.Prior;

import java.util.List;

public class EstimatorOptions {
    @Parameter(names = "-tree", description = "", required = true)
    public String tree;

    @Parameter(names = "-alignment", description = "", required = true)
    public String alignment;

    @Parameter(names = "-threads", description = "", required = false)
    public int threads = 1;

    @Parameter(names = "-distributed", description = "", required = false)
    public boolean distributed;

    @ParametersDelegate
    public GeneticCodeOption gc = new GeneticCodeOption();

    @Parameter(names = "-hosts", description = "", required = false, variableArity = true)
    public List<String> hosts = Lists.newArrayList();

    @Parameter(names = "-hostsfile", required = false)
    public String hostsFile;

    @Parameter(names = "-checkpoint", required = false)
    public String checkpointFile;

    @Parameter(names = "-prior", description = "Comma-separated prior and prior parameters (e.g. normal,10 or dirichlet,2.0)", required = false, converter = PriorConverter.class)
    public Prior prior;

    /*
    TODO: allow command-line options for:
    1. skip estimation of global parameters
    2. skip estimation of branch length parameters
    4. if you've specified sites, you can't estimate globals or branch lengths (i.e. these are fixed)
    */
}
