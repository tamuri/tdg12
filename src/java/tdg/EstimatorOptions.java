package tdg;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.google.common.collect.Lists;
import tdg.cli.GeneticCodeOption;

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

    /*
    TODO: allow command-line options for:
    1. skip estimation of global parameters
    2. skip estimation of branch length parameters
    3. provide global parameters (as initial parameter value or for fixed globals)
    4. if you've specified sites, you can't estimate globals or branch lengths (i.e. these are fixed)
    */
}
