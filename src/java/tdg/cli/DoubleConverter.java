package tdg.cli;

import com.beust.jcommander.converters.BaseConverter;

/**
 * @author Asif Tamuri
 * @version $Id: DoubleConverter.java 149 2010-08-14 18:48:07Z tamuri $
 */
public class DoubleConverter extends BaseConverter<Double> {
    public DoubleConverter(String value) { super(value); }

    @Override
    public Double convert(String value) {
        return Double.parseDouble(value);
    }
}
