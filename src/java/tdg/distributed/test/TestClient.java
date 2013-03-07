package tdg.distributed.test;

import com.caucho.hessian.client.HessianProxyFactory;
import tdg.model.Fitness;

/**
 * User: atamuri
 * Date: 07/03/2013 00:03
 */
public class TestClient {
    public static void main(String[] args) throws Exception {
        String url = "http://localhost:8080/service";

        HessianProxyFactory factory = new HessianProxyFactory();
        final TestAPI service = (TestAPI) factory.create(TestAPI.class, url);

        service.checkFitness(new Fitness(new double[]{0}, true));
        service.checkFitness(Fitness.getMutationOnlyFitness());
    }
}
