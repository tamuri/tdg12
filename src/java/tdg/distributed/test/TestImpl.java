package tdg.distributed.test;

import tdg.model.Fitness;

/**
 * User: atamuri
 * Date: 07/03/2013 00:04
 */
public class TestImpl implements TestAPI {
    @Override
    public void checkFitness(Fitness f) {
        System.out.printf("%s\n", f == Fitness.getMutationOnlyFitness());
    }
}
