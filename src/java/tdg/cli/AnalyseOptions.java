package tdg.cli;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import tdg.model.Prior;
import tdg.utils.Pair;

/**
 * AnalyseOptions class used by JCommander to display and load command-line options.
 *
 * @author Asif Tamuri (atamuri@nimr.mrc.ac.uk)
 */
public class AnalyseOptions {
    @Parameter(names = "-t", description = "Tree file in Newick format", required = true)
    public String treeFile;

    @Parameter(names = "-s", description = "Codon alignment file in Phylip sequential format.", required = true)
    public String alignmentFile;

    @Parameter(names = "-site", description = "Site (e.g. 101) or range of sites (e.g. 101-202) to analyse.", required = false, converter = SiteRangeConverter.class)
    public Pair<Integer, Integer> site;

    @Parameter(names = "-heteroclades", description = "Perform analysis using the heterogeneous model with these comma-separated clade prefixes.")
    public String heteroClades = null;

    @Parameter(names = "-optimruns", description = "The number of times to run the optimisation, with random initial parameters.", hidden = true)
    public int optimRuns = 1;

    @Parameter(names = "-threads", description = "The number of threads to use.", required = false)
    public int threads = 1;

    @ParametersDelegate
    public GlobalsOptions globals = new GlobalsOptions();

    @ParametersDelegate
    public GeneticCodeOption gc = new GeneticCodeOption();

    @ParametersDelegate
    public ApproxOption approx = new ApproxOption();

    @Parameter(names = "-prior", description = "Comma-separated prior and prior parameters (e.g. normal,10 or dirichlet,2.0)", required = false, converter = PriorConverter.class)
    public Prior prior;

}
