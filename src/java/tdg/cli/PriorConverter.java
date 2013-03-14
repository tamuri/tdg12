package tdg.cli;

import com.beust.jcommander.IStringConverter;
import com.beust.jcommander.ParameterException;
import tdg.model.DirichletPrior;
import tdg.model.NormalPrior;
import tdg.model.Prior;

public class PriorConverter implements IStringConverter<Prior> {
    @Override
    public Prior convert(String s) {
        try {
            Prior p;
            String[] definition = s.split(",");
            if ("normal".equals(definition[0])) {
                p = new NormalPrior(Double.parseDouble(definition[1]));
            } else if ("dirichlet".equals(definition[0])) {
                p = new DirichletPrior(Double.parseDouble(definition[1]));
            }  else {
                throw new ParameterException("Could not create prior '" + s + "'.\n");
            }

            System.out.printf("tdg.cli.PriorConverter - Using %s.\n", p.toString());
            return p;

        } catch (Exception e) {
            throw new ParameterException("Could not create prior '" + s + "'.\n");
        }
    }
}
