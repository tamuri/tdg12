package tdg.cli;

import com.beust.jcommander.IStringConverter;

/**
 * @author Asif Tamuri
 * @version $Id: DoubleArrayConverter.java 149 2010-08-14 18:48:07Z tamuri $
 */
public class DoubleArrayConverter implements IStringConverter<double[]> {

  @Override
  public double[] convert(String value) {
      String[] s = value.split(",");
      double[] d = new double[s.length];
      for (int i = 0; i < s.length; i++) {
          d[i] = Double.parseDouble(s[i]);
      }
      return d;
  }

}
