package tdg.distributed;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.googlecode.hessianserver.HessianServer;
import com.googlecode.hessianserver.HessianServiceDefinition;
import pal.alignment.Alignment;
import tdg.cli.GeneticCodeOption;
import tdg.utils.CoreUtils;
import tdg.utils.PhyloUtils;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.net.InetAddress;
import java.util.Random;

/**
 * User: atamuri
 * Date: 05/03/2013 15:55
 */
public class Server {

    public static void main(String[] args) {
        Server s = new Server();
        s.run(args);
    }

    private void run(String[] args) {
        ServerOptions options = new ServerOptions();
        new JCommander(options).parse(args);

        Alignment alignment = PhyloUtils.readAlignment(options.alignment);

        ServiceImpl service = new ServiceImpl(alignment, options.threads);

        HessianServiceDefinition hessianService = new HessianServiceDefinition("/service", service, ServiceAPI.class);

        int port = CoreUtils.findFreePort();
        HessianServer hessianServer = new HessianServer(port, hessianService);

        try {
            writeHostnameFile(port);
        } catch (Exception e) {
            e.printStackTrace();
        }

        hessianServer.start();
    }

    private static void writeHostnameFile(int port) throws Exception {
        InetAddress in = InetAddress.getLocalHost();
        InetAddress[] all = InetAddress.getAllByName(in.getHostName());

        FileWriter writer = new FileWriter("host_" + in.getHostName() + "_" + new Random().nextInt(100000) +".txt");
        BufferedWriter out = new BufferedWriter(writer);

        for (InetAddress anAll : all) {
            out.write(anAll + ":" + port + "\n");
        }

        out.close();
    }


    private class ServerOptions {
        @Parameter(names = "-alignment", description = "", required = true)
        public String alignment;

        @Parameter(names = "-threads", description = "", required = false)
        public int threads = 1;

        @ParametersDelegate
        public GeneticCodeOption gc = new GeneticCodeOption();
    }
}
