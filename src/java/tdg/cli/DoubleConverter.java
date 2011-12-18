package tdg.cli;

import com.beust.jcommander.converters.BaseConverter;

/**
 * @author Asif Tamuri (atamuri@nimr.mrc.ac.uk)
 */
public class DoubleConverter extends BaseConverter<Double> {
    public DoubleConverter(String value) { super(value); }

    @Override
    public Double convert(String value) {
        return Double.parseDouble(value);
    }
}
