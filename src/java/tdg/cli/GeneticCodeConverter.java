package tdg.cli;

import com.beust.jcommander.IStringConverter;
import tdg.utils.GeneticCode;

/**
 * A helper for JCommander: takes the command-line 'gc' option and loads the correct GeneticCode instance.
 *
 * @author Asif Tamuri (atamuri@nimr.mrc.ac.uk)
 */
public class GeneticCodeConverter implements IStringConverter<GeneticCode> {
    @Override
    public GeneticCode convert(String value) {
        if ("vertebrate_mit".equals(value)) {
            GeneticCode.setCode(GeneticCode.VERTEBRATE_MITOCHONDRIAL_CODE);
        } else if ("standard".equals(value)) {
            GeneticCode.setCode(GeneticCode.STANDARD_CODE);
        } else if ("plastid".equals(value)) {
            GeneticCode.setCode(GeneticCode.PLASTID_CODE);
        } else {
            throw new RuntimeException("Unknown genetic code table '" + value + "'. Valid tables are 'standard' and 'vertebrate_mit'.\n");
        }

        System.out.printf("tdg.cli.GeneticCodeConverter - Set code to %s\n", value);

        return GeneticCode.getInstance();
    }
}
