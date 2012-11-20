package tdg.distributed;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.ParameterException;
import com.google.common.cache.CacheBuilder;
import com.google.common.cache.CacheLoader;
import com.google.common.cache.LoadingCache;
import com.google.common.collect.Lists;
import org.simpleframework.http.Request;
import org.simpleframework.http.Response;
import org.simpleframework.http.core.Container;
import org.simpleframework.transport.connect.Connection;
import org.simpleframework.transport.connect.SocketConnection;
import org.simpleframework.util.thread.Scheduler;
import pal.alignment.Alignment;
import pal.tree.Tree;
import tdg.SiteAnalyser;
import tdg.cli.AnalyseOptions;
import tdg.model.TDGGlobals;
import tdg.utils.CoreUtils;
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
import java.util.Random;
import java.util.concurrent.ExecutionException;

/**
 * @author Asif Tamuri (atamuri@nimr.mrc.ac.uk)
 */
public class Slave implements Container {
    private Scheduler scheduler;
    private AnalyseOptions options;
    private Tree tree;
    private Alignment alignment;
    private TDGGlobals globals;

    LoadingCache<Integer, SiteAnalyser> siteAnalyserCache = CacheBuilder.newBuilder().maximumSize(1000).build(
            new CacheLoader<Integer, SiteAnalyser>() {
                @Override
                public SiteAnalyser load(Integer site) throws Exception {
                    SiteAnalyser sa = new SiteAnalyser(tree, alignment, globals, site, options);
                    sa.run();
                    return sa;
                }
            }
    );



    private class Task implements Runnable {
        private final Response response;
        private final Request request;

        public Task(Request request, Response response) {
            this.response = response;
            this.request = request;
        }

        @Override
        public void run() {
            // long startTime = System.currentTimeMillis();

            try {
                int site = Integer.parseInt(request.getParameter("site"));

                // System.out.printf("Site %s-- Task started.\n", site);

                SiteAnalyser sa;

                try {
                    sa = siteAnalyserCache.get(site);

                } catch (ExecutionException e) {
                    e.printStackTrace();
                    throw new RuntimeException(e);
                }

                PrintStream body = response.getPrintStream();
                response.set("Content-Type", "text/plain");

                // body.println(sa.getHomogeneousLikelihood());
                double l = sa.evaluate(globals);
                body.println(l);
                body.close();

                // System.out.printf("Site %s-- Task complete in %s ms.\n", site, System.currentTimeMillis() - startTime);
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }

    public Slave(Scheduler scheduler, AnalyseOptions options) {
        this.scheduler = scheduler;
        this.options = options;

        this.tree = PhyloUtils.readTree(options.treeFile);
        this.alignment = PhyloUtils.readAlignment(options.alignmentFile);
        this.globals = new TDGGlobals(options.globals.tau, options.globals.kappa, options.globals.pi, options.globals.mu, options.globals.gamma);
    }

    @Override
    public void handle(Request request, Response response) {
        try {
            if (request.getQuery().containsKey("updateglobals")) {


                List<Double> g = Lists.transform(Arrays.asList(request.getParameter("updateglobals").split(",")), Functions.stringToDouble());

                globals = new TDGGlobals(g.get(0), g.get(1), new double[]{g.get(2), g.get(3), g.get(4), g.get(5)}, g.get(6), g.get(7));
                System.out.printf("Set %s.\n", globals.toString());

                PrintStream body = response.getPrintStream();
                response.set("Content-Type", "text/plain");
                body.print("OK.");
                body.close();

                if (request.getParameter("refresh").equals("true")) {
                    siteAnalyserCache.invalidateAll();
                    System.out.println("Refreshed cache.");
                }

                // siteAnalyserCache.invalidateAll();


            } else {
                int site = Integer.parseInt(request.getParameter("site"));
                // System.out.printf("Site %s - Received request. Creating task and adding to queue.\n", site);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        Task task = new Task(request, response);
        scheduler.execute(task);
    }

    public static void main(String... args) throws Exception {
        AnalyseOptions options = new AnalyseOptions();
        JCommander jc = new JCommander(options);

        try {
            jc.parse(args);
        } catch (ParameterException e) {
            System.out.printf("Error: %s\n", e.getMessage());
            jc.usage();
            System.exit(0);
        }

        Scheduler scheduler = new Scheduler(options.threads);
        Container container = new Slave(scheduler, options);
        Connection connection = new SocketConnection(container);

        int freePort = CoreUtils.findFreePort();

        SocketAddress address = new InetSocketAddress(freePort);
        connection.connect(address);

        System.out.printf("%s started on port %s with %s worker thread(s).\n", Slave.class.getName(), freePort, options.threads);
        writeHostnameFile(freePort);
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

}
