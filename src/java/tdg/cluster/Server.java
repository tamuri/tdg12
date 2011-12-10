package tdg.cluster;

import com.beust.jcommander.JCommander;
import com.google.common.base.Joiner;
import com.google.common.collect.Lists;
import org.simpleframework.http.Request;
import org.simpleframework.http.Response;
import org.simpleframework.http.core.Container;
import org.simpleframework.transport.connect.Connection;
import org.simpleframework.transport.connect.SocketConnection;
import org.simpleframework.util.thread.Scheduler;
import pal.alignment.Alignment;
import pal.tree.Tree;
import tdg.Options;
import tdg.analysis.SiteAnalyser;
import tdg.models.TDGGlobals;
import tdg.utils.Functions;
import tdg.utils.PhyloUtils;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.net.InetAddress;
import java.net.InetSocketAddress;
import java.net.SocketAddress;
import java.util.Arrays;
import java.util.List;

/**
 * @author Asif Tamuri
 */
public class Server implements Container {

    static final int SERVER_PORT = 9090; // TODO: Put in config file
    private Scheduler scheduler;
    private Options options;
    private Tree tree;
    private Alignment alignment;
    private TDGGlobals globals;


    private class Task implements Runnable {
        private final Response response;
        private final Request request;

        public Task(Request request, Response response) {
            this.response = response; 
            this.request = request;
        }

        @Override
        public void run() {
            long startTime = System.currentTimeMillis();

            try {
                int site = Integer.parseInt(request.getParameter("site"));

                System.out.printf("Site %s-- Task started.\n", site);

                SiteAnalyser sa = new SiteAnalyser(tree, alignment, globals, site, options);
                sa.run();

                PrintStream body = response.getPrintStream();
                response.set("Content-Type", "text/plain");

                body.println(sa.getHomogeneousLikelihood());
                body.close();

                System.out.printf("Site %s-- Task complete in %s ms.\n", site, System.currentTimeMillis() - startTime);
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }

    public Server(Scheduler scheduler, Options options) {
        this.scheduler = scheduler;
        this.options = options;

        this.tree = PhyloUtils.readTree(options.treeFile);
        this.alignment = PhyloUtils.readAlignment(options.alignmentFile);
        this.globals = new TDGGlobals(options.tau, options.kappa, options.pi, options.mu, options.gamma);
    }

    @Override
    public void handle(Request request, Response response) {
        try {
            if (request.getQuery().containsKey("updateglobals")) {

                List<Double> g = Lists.transform(Arrays.asList(request.getParameter("updateglobals").split(",")), Functions.stringToDouble());

                globals = new TDGGlobals(g.get(0), g.get(1), new double[]{g.get(2), g.get(3), g.get(4), g.get(5)}, g.get(6), g.get(7));
                System.out.printf("Successfully updated TDGGlobals { %s }.", Joiner.on(", ").join(g));

                PrintStream body = response.getPrintStream();
                response.set("Content-Type", "text/plain");
                body.print("OK.");
                body.close();
            } else {
                int site = Integer.parseInt(request.getParameter("site"));
                System.out.printf("Site %s - Received request. Creating task and adding to queue.\n", site);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        Task task = new Task(request, response);
        scheduler.execute(task);
    }

    public static void main(String... args) throws Exception {
        Options options = new Options();
        JCommander jc = new JCommander(options);
        if (args.length == 0) {
            jc.usage();
            System.out.println("Options preceded by an asterisk are required.");
            System.exit(0);
        } else {
            jc.parse(args);
        }

        Scheduler scheduler = new Scheduler(options.threads);
        Container container = new Server(scheduler, options);
        Connection connection = new SocketConnection(container);
        SocketAddress address = new InetSocketAddress(SERVER_PORT);
        connection.connect(address);

        System.out.printf("%s started on port %s with %s worker thread(s).\n", Server.class.getName(), SERVER_PORT, options.threads);
        writeHostnameFile();
    }

    private static void writeHostnameFile() throws Exception {
        InetAddress in  = InetAddress.getLocalHost();
        InetAddress[] all = InetAddress.getAllByName(in.getHostName());

        FileWriter writer = new FileWriter("host" + in.getHostName() + ".txt");
        BufferedWriter out = new BufferedWriter(writer);

        for (InetAddress anAll : all) {
            out.write(anAll + "\n");
        }

        out.close();
    }

}
