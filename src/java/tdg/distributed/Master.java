package tdg.distributed;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParameterException;
import com.beust.jcommander.ParametersDelegate;
import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.google.common.io.Files;
import com.google.common.primitives.Doubles;
import com.ning.http.client.AsyncCompletionHandler;
import com.ning.http.client.AsyncHttpClient;
import com.ning.http.client.AsyncHttpClientConfig;
import com.ning.http.client.Response;
import com.ning.http.client.extra.ThrottleRequestFilter;
import org.apache.commons.math.FunctionEvaluationException;
import org.apache.commons.math.analysis.MultivariateRealFunction;
import org.apache.commons.math.optimization.GoalType;
import org.apache.commons.math.optimization.RealPointValuePair;
import org.apache.commons.math.optimization.SimpleScalarValueChecker;
import org.apache.commons.math.optimization.direct.DirectSearchOptimizer;
import org.apache.commons.math.optimization.direct.NelderMead;
import pal.alignment.Alignment;
import tdg.Constants;
import tdg.cli.GeneticCodeOption;
import tdg.cli.GlobalsOptions;
import tdg.utils.CoreUtils;
import tdg.utils.PhyloUtils;
import tdg.utils.ValueComparer;

import java.io.File;
import java.io.FilenameFilter;
import java.nio.charset.Charset;
import java.util.*;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

/**
 * @author Asif Tamuri (atamuri@nimr.mrc.ac.uk)
 */
public class Master {

    private Collection<Integer> sites;
    private List<String> servers;
    private boolean refresh = false;
    AsyncHttpClient asyncHttpClient;
    private int totalits = 1;
    private List<Result> results = Lists.newArrayList();

    class Result {
        public int iteration;
        public double lnL;
        public double[] point;

        Result(int iteration, double lnL, double[] point) {
            this.iteration = iteration;
            this.lnL = lnL;
            this.point = point;
        }

        @Override
        public String toString() {
            return "Result{" +
                    iteration + "." +
                    " lnL =\t" + lnL +
                    ", point =\t" + Doubles.join("\t", point) +
                    '}';
        }
    }

    class GlobalParameterOptimiser implements MultivariateRealFunction {
        private int evals = 1;

        @Override
        public double value(double[] point) throws FunctionEvaluationException, IllegalArgumentException {
            System.out.printf("Evaluation %s / %s\n", evals++, totalits++);

            // use very bad log-likelihood to exit quickly if we're near constraint
            if (point[0] <= 0) // tau
                return Constants.VERY_BAD_LIKELIHOOD;

            if (point[1] <= 0) // kappa
                return Constants.VERY_BAD_LIKELIHOOD;

            if (point[5] <= 0) // mu
                return Constants.VERY_BAD_LIKELIHOOD;

            double[] pi = CoreUtils.alr_inv(new double[]{point[2], point[3], point[4]});
            return evaluate(point[0], point[1], pi[0], pi[1], pi[2], point[5]);
        }
    }

    private void optimiseBranchLengths() {
        // We assume that each of the slaves has a cached SiteAnalyser with the correct globals and fitness
        // i.e. we've estimated global parameters already

        // while the sum log-likelihoods of all sites is reducing

            // for every internal node in the tree

                // re-root the tree with this internal node


                // update each slave with the new re-rooted tree
                // (as part of updating tree, each slave should calculate likelihood of site again
                // and save conditionals of child nodes)

                // for each child of this internal node





    }

    private double evaluate(double... point) {
        long startTime = System.currentTimeMillis();
        System.out.printf("Evaluate for: %s\n", Doubles.join(", ", point));

        Set<Future<Response>> results = Sets.newHashSet();
        double sumlnL = 0.0;


        try {
            // Update TdG globals on each server
            updateGlobalsOnServers(servers, asyncHttpClient,
                    point[0], // tau
                    point[1], // kappa
                    point[2], point[3], point[4], 1 - point[2] - point[3] - point[4], // pi_T, pi_C, pi_A, pi_G (= 1 - pi_T - pi_C - pi_A)
                    point[5], // mu
                    0.0); // gamma

            sendRequestToServers(sites.iterator(), servers, asyncHttpClient, results);

            for (Future<Response> r : results) {
                String response = null;
                try {
                    response = r.get().getResponseBody();
                    // System.out.printf("response: %s\n\n", response);
                } catch (Exception e) {
                    System.err.println("ERROR: couldn't read response.");
                    e.printStackTrace();
                }

                sumlnL += Double.parseDouble(response);
                // System.out.printf("%s\n", response);
            }

        } catch (Exception e) {
            e.printStackTrace();
        }


        System.out.printf("Log-likelihood = %s\n", sumlnL);
        System.out.printf("Time: %s ms\n\n", System.currentTimeMillis() - startTime);

        return sumlnL;
    }

    private void updateGlobalsOnServers(List<String> servers, AsyncHttpClient asyncHttpClient, double tau, double kappa, double t, double c, double a, double g, double mu, double gamma) throws Exception {
        Set<Future<Response>> results = Sets.newHashSet();

        if (refresh)
            System.out.println("Refreshing slave cache.");

        for (String server : servers) {
            final String url = "http://" + server + "/?updateglobals=" + Doubles.join(",", tau, kappa, t, c, a, g, mu, gamma) + "&refresh=" + refresh;
            Future<Response> f = asyncHttpClient.prepareGet(url).execute();
            results.add(f);
        }

        refresh = false;

        for (Future<Response> r : results) {
            String response = r.get().getResponseBody();
            // System.out.printf("%s\n", response);
        }

    }

    private void sendRequestToServers(Iterator<Integer> siteIterator, List<String> servers, AsyncHttpClient asyncHttpClient, Collection<Future<Response>> results) throws Exception {
        // Keep cycling over list of servers
        Iterator<String> cycleServers = Iterators.cycle(servers);

        // Loop over every site, in Q matrix size descending order
        while (siteIterator.hasNext()) {
            final int site = siteIterator.next();

            // Pick a server and construct request URL
            final String url = "http://" + cycleServers.next() + "/?site=" + site;

            // Send a request to the selected server
            Future<Response> f = asyncHttpClient.prepareGet(url).execute(
                    new AsyncCompletionHandler<Response>() {
                        @Override
                        public Response onCompleted(Response response) throws Exception {
                            // System.out.printf("Response from %s for site %s. lnL = %s\n", url, site, response.getResponseBody().split(",")[0]);
                            return response;
                        }
                    });

            // Add the Future object to our tdg.results.results collection
            results.add(f);
        }
    }

    /**
     * Configure and intialise the asynchronous http client to use for all
     * client->server requests
     *
     * @return initialised and configured AsyncHttpClient object
     */
    private AsyncHttpClient getAsyncHttpClient() {
        // Configure and startTime the http client, not too many requests, not too many native threads (e.g. OSX barfs on 2560)
        AsyncHttpClientConfig.Builder b = new AsyncHttpClientConfig.Builder()
                //.setMaximumConnectionsPerHost(maxConnectionsPerHost)
                //.setMaximumConnectionsTotal(servers.size() * 20)
                .setRequestTimeoutInMs(requestTimeout)
                .setAllowPoolingConnection(true)
                .addRequestFilter(new ThrottleRequestFilter(servers.size() * 20))
                .setExecutorService(Executors.newFixedThreadPool(100));

        return new AsyncHttpClient(b.build());
    }

    /**
     * @return a list of server addresses running the TdGServer
     * @throws Exception IOException
     */
    private List<String> getServerAddresses() throws Exception {
        // Get all the TdGServer hostname files
        String hostnameFilesPath = "./";
        List<File> allHostFiles = Lists.newArrayList(
                new File(hostnameFilesPath).listFiles(new FilenameFilter() { // could use commons-io: new WildcardFileFilter("host*.txt");

                    @Override
                    public boolean accept(File file, String s) {
                        return s.startsWith("host_") && s.endsWith(".txt");
                    }
                }));

        List<String> hostAddresses = Lists.newArrayList();

        for (File hostFile : allHostFiles) {
            // address written in the form: beo-25/192.168.52.25
            String address = Files.readFirstLine(hostFile, Charset.forName("US-ASCII")).split("/")[1];
            hostAddresses.add(address);
            System.out.printf("Slave server address: %s\n", address);
        }

        return hostAddresses;
    }

    /**
     * Reads the alignment file and returns an iterator of sites
     * ordered by the number of codons that occur at that site, descending
     *
     * @return an iterator of site positions
     */
    private Collection<Integer> getSites() {
        // Load the alignment
        final Alignment alignment = PhyloUtils.readAlignment(alignmentFile);

        // Get the Q matrix size (= number of codons) for each site
        Map<Integer, Integer> siteCodonCount = new HashMap<Integer, Integer>();
        for (int i = 1; i <= (alignment.getSiteCount() / 3); i++) {
            Map<String, Integer> sitePattern = PhyloUtils.getCodonsAtSite(alignment, i);
            List<Integer> aminoAcidsAtSite = PhyloUtils.getDistinctAminoAcids(sitePattern.values());
            List<Integer> allCodons = PhyloUtils.getCodonsFromAminoAcids(aminoAcidsAtSite);

            if (allCodons.size() < codonCutoff)
                siteCodonCount.put(i, allCodons.size());
        }

        // Order sites by number of codons at site, descending
        SortedMap<Integer, Integer> siteCodonCountDesc = new TreeMap<Integer, Integer>(new ValueComparer<Integer,Integer>(siteCodonCount));
        siteCodonCountDesc.putAll(siteCodonCount);

        // DEBUG: Here's what we get:
        /*System.out.printf("Order of sites: ");
        for (Map.Entry<Integer, Integer> e : siteCodonCountDesc.entrySet()) {
            System.out.printf("%s (%s), ", e.getKey(), e.getValue());
        }
        System.out.println();
*/
        // to iterate over sites in their natural order:
        //return siteCodonCount.keySet().iterator();

        return siteCodonCountDesc.keySet();
    }


    public static void main(String[] args) throws Exception {
        // TODO: There should be a 'BaseClient' and then implementations of 'optimise', or task etc.
        // e.g. simple optimise, optim one variable (mu), optim multiple variables (pi + kappa)
        Master c = new Master();
        JCommander jc = new JCommander(c);

        try {
            jc.parse(args);
        } catch (ParameterException e) {
            System.out.printf("Error: %s\n", e.getMessage());
            jc.usage();
            System.exit(0);
        }
        c.run();

    }

    private void run() throws Exception {
        sites = getSites();
        servers = getServerAddresses();

        double[] start = new double[]{globals.tau, globals.kappa, globals.pi[0], globals.pi[1], globals.pi[2], globals.mu};

        asyncHttpClient = getAsyncHttpClient();

        if (optimiseGlobals) {
            int iteration = 0;
            while (iteration++ < 5) {
                System.out.printf("Iteration %s for tolerance %s.\n", iteration, convergence);

                // 6 global parameters to optimise: tau, kappa, pi (3) and mu
                DirectSearchOptimizer dso = new NelderMead();
                dso.setConvergenceChecker(new SimpleScalarValueChecker(-1, convergence));
                dso.setMaxEvaluations(Constants.MAX_EVALUATIONS);

                RealPointValuePair pair;

                // use additive log-transform to ensure that pi parameters sum to 1
                double[] pi = CoreUtils.alr(new double[]{start[2], start[3], start[4]});
                start = new double[]{start[0], start[1], pi[0], pi[1], pi[2], start[5]};

                try {
                    pair = dso.optimize(new GlobalParameterOptimiser(), GoalType.MAXIMIZE, start);
                } catch (FunctionEvaluationException me) {
                    System.out.println("Reached maximum number of evaluations!");
                    throw new RuntimeException(me);
                }

                double[] optima = pair.getPoint();
                System.out.printf("Found optima.\nLog-likelihood = %s\n", pair.getValue());
                pi = CoreUtils.alr_inv(Arrays.copyOfRange(optima, 2, 5));
                System.out.printf("MLE: -tau %s -kappa %s -pi %s,%s,%s,%s -mu %s\n", optima[0], optima[1], pi[0], pi[1], pi[2], (1 - pi[0] - pi[1] - pi[2]), optima[5]);

                System.arraycopy(optima, 0, start, 0, optima.length);
                start = new double[]{start[0], start[1], pi[0], pi[1], pi[2], start[5]};

                results.add(new Result(iteration, pair.getValue(), Arrays.copyOf(start, start.length)));
                refresh = true;
            }

            System.out.println("Results of optimisation iterations:");
            for (Result r : results) {
                System.out.printf("%s\n", r.toString());
            }

        } else {
            evaluate(start);
        }

        System.out.printf("FINISHED!\n");
        asyncHttpClient.close();
    }

    @Parameter(names = "-s", description = "Codon alignment file in Phylip sequential format", required = true)
    public String alignmentFile;

    @ParametersDelegate
    public GeneticCodeOption gc = new GeneticCodeOption();

    @ParametersDelegate
    public GlobalsOptions globals = new GlobalsOptions();

    @Parameter(names = "-timeout", description = "How long to wait before timing out master-to-slave connection", required = false)
    public int requestTimeout = 60000000;

    @Parameter(names = "-max-connections", description = "Maximum number of connections from master to each slave host", required = false)
    public int maxConnectionsPerHost = 10;

    @Parameter(names = "-thread-pool-size", description = "How many threads to have available to connect to each slave", required = false)
    public int threadPoolSize = 100;

    @Parameter(names = "-optimise-globals", description = "Optimise the global parameters", required = false)
    public boolean optimiseGlobals = false;

    @Parameter(names = "-codon-cutoff", description = "Only analyse sites with this maximum number of codons observed (just for testing).", required = false, hidden = true)
    public int codonCutoff = 100;

    @Parameter(names = "-convergence", description = "Convergence tolerance.", required = false, hidden = true)
    public double convergence = 1E-5;
}
