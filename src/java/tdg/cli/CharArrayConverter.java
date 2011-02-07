package tdg.cli;

import com.beust.jcommander.IStringConverter;

/**
 * @author Asif Tamuri
 * @version $Id: CharArrayConverter.java 149 2010-08-14 18:48:07Z tamuri $
 */
public class CharArrayConverter implements IStringConverter<char[]> {

    @Override
    public char[] convert(String value) {
        String[] s = value.split(",");
        char[] c = new char[s.length];
        for (int i = 0; i < s.length; i++) {
            c[i] = s[i].toCharArray()[0];
        }
        return c;
    }

}
