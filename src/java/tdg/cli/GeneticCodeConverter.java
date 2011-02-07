package tdg.cli;

import com.beust.jcommander.IStringConverter;
import tdg.utils.GeneticCode;

/**
 * @author Asif Tamuri
 * @version $Id$
 */
public class GeneticCodeConverter implements IStringConverter<GeneticCode> {
    @Override
    public GeneticCode convert(String value) {
        if (value.equals("vertebrate_mit")) {
            GeneticCode.initialise(GeneticCode.VERTEBRATE_MITOCHONDRIAL_CODE);
        } else if (value.equals("standard")) {
            GeneticCode.initialise(GeneticCode.STANDARD_CODE);
        } else {
            throw new RuntimeException("Unknown genetic code table '"+value+"'. Valid tables are 'standard' and 'vertebrate_mit'.\n");
        }

        return GeneticCode.getInstance();
    }
}
