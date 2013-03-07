package tdg.distributed.test;

import com.googlecode.hessianserver.HessianServer;
import com.googlecode.hessianserver.HessianServiceDefinition;

/**
 * User: atamuri
 * Date: 07/03/2013 00:03
 */
public class TestServer {
    public static void main(String[] args) {
        HessianServiceDefinition hessianService = new HessianServiceDefinition("/service", new TestImpl(), TestAPI.class);
        HessianServer hessianServer = new HessianServer(8080, hessianService);
        hessianServer.start();
    }
}
