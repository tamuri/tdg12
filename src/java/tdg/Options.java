package tdg;

import com.beust.jcommander.Parameter;
import tdg.cli.DoubleArrayConverter;
import tdg.cli.DoubleConverter;
import tdg.cli.GeneticCodeConverter;
import tdg.utils.GeneticCode;

/**
 * @author Asif Tamuri
 * @version $Id$
 */
public class Options {
    @Parameter(names = "-t", description = "Tree file in Newick format", required = true)
    public String treeFile;

    @Parameter(names = "-s", description = "Codon alignment file in Phylip sequential format.", required = true)
    public String alignmentFile;

    @Parameter(names = "-site", description = "Location of site to analyse (0 = all sites).", required = false)
    public int site = 0;

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

    @Parameter(names = "-homogonly", description = "Run homogeneous model only.")
    public boolean homogonly = false;

    @Parameter(names = "-approx", description = "Use the approximate method to optimise the likelihood")
    public boolean approx = false;

    @Parameter(names = "-optimruns", description = "The number of times to run the optimisation, with random initial parameters.", hidden = true)
    public int optimRuns = 1;

    @Parameter(names = "-gc", description = "The genetic code translation to use (standard or vertebrate_mit).", required = true, converter = GeneticCodeConverter.class)
    public GeneticCode geneticCode;

    @Parameter(names = "-threads", description = "The number of threads to use.", required = false)
    public int threads = 1;
}
