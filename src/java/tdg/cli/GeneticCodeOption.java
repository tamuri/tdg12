package tdg.cli;

import com.beust.jcommander.Parameter;
import tdg.utils.GeneticCode;

/**
 * @author Asif Tamuri (atamuri@nimr.mrc.ac.uk)
 */
public class GeneticCodeOption {
    @Parameter(names = "-gc", description = "The genetic code translation to use (standard or vertebrate_mit).", required = true, converter = GeneticCodeConverter.class)
    public GeneticCode geneticCode;
}
