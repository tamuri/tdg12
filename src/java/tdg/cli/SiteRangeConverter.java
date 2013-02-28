package tdg.cli;

import com.beust.jcommander.IStringConverter;
import tdg.utils.Pair;

/**
 * Takes a string and returns a Pair&lt;Integer,Integer&gt; indicating a range
 */
public class SiteRangeConverter implements IStringConverter<Pair<Integer,Integer>> {

    private static final String SEP = "-";
    @Override
    public Pair<Integer,Integer> convert(String s) {
        int pos = s.indexOf(SEP);

        // if this is a range
        if (pos >= 0) {
            String[] split = s.split(SEP);
            return Pair.of(Integer.parseInt(split[0]), Integer.parseInt(split[1]));
        } else {
            // it's a single site
            return Pair.of(Integer.parseInt(s), Integer.parseInt(s));
        }
    }
}
