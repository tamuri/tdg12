package tdg.cli;

import com.beust.jcommander.IStringConverter;
import com.beust.jcommander.ParameterException;
import tdg.model.DirichletPrior;
import tdg.model.NormalPrior;
import tdg.model.Prior;

public class PriorConverter implements IStringConverter<Prior> {
    @Override
    public Prior convert(String s) {
        System.out.printf("tdg.cli.PriorConverter: Using prior '%s'.\n", s);
        try {
            String[] definition = s.split(",");
            if ("normal".equals(definition[0])) {
                return new NormalPrior(Double.parseDouble(definition[1]));
            } else if ("dirichlet".equals(definition[0])) {
                return new DirichletPrior(Double.parseDouble(definition[1]));
            }  else {
                throw new ParameterException("Could not create prior '" + s + "'.\n");
            }
        } catch (Exception e) {
            throw new ParameterException("Could not create prior '" + s + "'.\n");
        }
    }
}
