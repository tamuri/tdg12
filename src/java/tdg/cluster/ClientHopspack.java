package tdg.cluster;

import com.google.common.base.Charsets;
import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.google.common.io.Files;
import com.google.common.primitives.Doubles;
import com.ning.http.client.AsyncCompletionHandler;
import com.ning.http.client.AsyncHttpClient;
import com.ning.http.client.AsyncHttpClientConfig;
import com.ning.http.client.Response;
import org.apache.commons.configuration.Configuration;
import org.apache.commons.configuration.PropertiesConfiguration;
import pal.alignment.Alignment;
import tdg.utils.GeneticCode;
import tdg.utils.PhyloUtils;

import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.util.*;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

/**
 * @author Asif Tamuri
 */
public class ClientHopspack {
    private String alignmentFile, hostnameFilesPath;
    private int maxConnectionsPerHost, requestTimeout, threadPoolSize;

    private void run(String[] args) throws Exception {

        Collection<Integer> sites = getSites();
        List<String> servers = getServerAddresses();


        // Read HOPSPACK input file
        List<String> lines = Files.readLines(new File(args[0]), Charsets.UTF_8);

        double kappa = Double.parseDouble(lines.get(2));
        double pi1   = Double.parseDouble(lines.get(3));
        double pi2   = Double.parseDouble(lines.get(4));
        double pi3   = Double.parseDouble(lines.get(5));
        double pi4   = Double.parseDouble(lines.get(6));
        double tau   = Double.parseDouble(lines.get(7));
        double mu    = Double.parseDouble(lines.get(8));
        double gamma = Double.parseDouble(lines.get(9));


        // Update TdG globals on each server
        AsyncHttpClient asyncHttpClient = getAsyncHttpClient();

        try {
            // updateGlobalsOnServers(servers, asyncHttpClient, tau, 8.09325, 0.15828, 0.20551, 0.57731, 0.05890, mu, gamma);
            updateGlobalsOnServers(servers, asyncHttpClient, tau, kappa, pi1, pi2, pi3, pi4, mu, gamma);
        } catch (Exception e) {
            e.printStackTrace();
            throw new RuntimeException(e);
        }

        asyncHttpClient.getConfig().executorService().shutdown();
        asyncHttpClient = getAsyncHttpClient();


        long startTime = System.currentTimeMillis();

        Set<Future<Response>> results = Sets.newHashSet();

        sendRequestToServers(sites.iterator(), servers, asyncHttpClient, results);

        // Output each result
        double sumlnL = 0.0;
        for (Future<Response> r : results) {
            String response = null;
                try {
                    response = r.get().getResponseBody();
                } catch (Exception e) {
                    System.err.println("Error: couldn't read response.");
                    e.printStackTrace();
                }

                sumlnL += Double.parseDouble(response);
                System.out.printf("%s\n", response);
        }

        // We're done
        System.out.printf("Total time: %s ms\n", System.currentTimeMillis() - startTime);

        String out = String.format("1\n%s\n", sumlnL);
        System.out.printf("\n ****\nf(x) = %s \nGLOBALS:\ntau=%s, mu=%s, gamma=%s\n ****\n", sumlnL, tau, mu, gamma);

        // Write HOPSPACK output file
        Files.write(out, new File(args[1]), Charsets.UTF_8);


        // Close the http client thread pool and close all connections
        asyncHttpClient.getConfig().executorService().shutdown();


    }

    private void updateGlobalsOnServers(List<String> servers, AsyncHttpClient asyncHttpClient, double tau, double kappa, double t, double c, double a, double g, double mu, double gamma) throws Exception {
        Set<Future<Response>> results = Sets.newHashSet();

        for (String server : servers) {
            final String url = "http://" + server + ":" + Server.SERVER_PORT + "/?updateglobals=" + Doubles.join(",", tau, kappa, t, c, a, g, mu, gamma);
            Future<Response> f = asyncHttpClient.prepareGet(url).execute();
            results.add(f);
        }

        for (Future<Response> r : results) {
            String response = r.get().getResponseBody();
            System.out.printf("%s\n", response);
        }

    }

    private void sendRequestToServers(Iterator<Integer> siteIterator, List<String> servers, AsyncHttpClient asyncHttpClient, Collection<Future<Response>> results) throws Exception {
        // Keep cycling over list of servers
        Iterator<String> cycleServers = Iterators.cycle(servers);

        // Loop over every site, in Q matrix size descending order
        while (siteIterator.hasNext()) {
            final int site = siteIterator.next();

            // Pick a server and construct request URL
            final String url = "http://" + cycleServers.next() + ":" + Server.SERVER_PORT + "/?site=" + site;

            // Send a request to the selected server
            Future<Response> f = asyncHttpClient.prepareGet(url).execute(
                    new AsyncCompletionHandler<Response>() {
                        @Override
                        public Response onCompleted(Response response) throws Exception {
                            System.out.printf("Response from %s for site %s. lnL = %s\n", url, site, response.getResponseBody().split(",")[0]);
                            return response;
                        }
                    });

            // Add the Future object to our results collection
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
                .setMaximumConnectionsPerHost(maxConnectionsPerHost)
                .setRequestTimeoutInMs(requestTimeout)
                .setExecutorService(Executors.newFixedThreadPool(threadPoolSize));
        return new AsyncHttpClient(b.build());
    }

    /**
     * @return a list of server addresses running the TdGServer
     * @throws Exception IOException
     */
    private List<String> getServerAddresses() throws Exception {
        // Get all the TdGServer hostname files
        List<File> allHostFiles = Lists.newArrayList(
                new File(hostnameFilesPath).listFiles(new FilenameFilter() { // could use commons-io: new WildcardFileFilter("host*.txt");
                    @Override
                    public boolean accept(File file, String s) {
                        return s.startsWith("host") && s.endsWith(".txt");
                    }
                }));

        List<String> hostAddresses = Lists.newArrayList();

        for (File hostFile : allHostFiles) {
            // address written in the form: beo-25/192.168.52.25
            hostAddresses.add(Files.readFirstLine(hostFile, Charset.forName("US-ASCII")).split("/")[1]);
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

        List<String> sitesToDo;
        // Read in the sites not done
        try {
            sitesToDo = Files.readLines(new File("sites.notdone"), Charsets.UTF_8);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }

        // Get the Q matrix size (= number of codons) for each site
        Map<Integer, Integer> siteCodonCount = new HashMap<Integer, Integer>();
        for (int i = 1; i <= (alignment.getSiteCount() / 3); i++) {
            if (!sitesToDo.contains(Integer.toString(i))) continue;
            Map<String, Integer> sitePattern = PhyloUtils.getCodonsAtSite(alignment, i);
            List<Integer> aminoAcidsAtSite = PhyloUtils.getDistinctAminoAcids(sitePattern.values());
            List<Integer> allCodons = PhyloUtils.getCodonsFromAminoAcids(aminoAcidsAtSite);
            siteCodonCount.put(i, allCodons.size());
        }

        // Order sites by number of codons at site, descending
        SortedMap<Integer, Integer> siteCodonCountDesc = new TreeMap<Integer, Integer>(new ValueComparer<Integer, Integer>(siteCodonCount));
        siteCodonCountDesc.putAll(siteCodonCount);

        // DEBUG: Here's what we get:
        System.out.printf("Order of sites: ");
        for (Map.Entry<Integer, Integer> e : siteCodonCountDesc.entrySet()) {
            System.out.printf("%s (%s), ", e.getKey(), e.getValue());
        }
        System.out.println();

        // to iterate over sites in their natural order:
        //return siteCodonCount.keySet().iterator();

        return siteCodonCountDesc.keySet();
    }

    private void loadConfiguration() throws Exception {
        Configuration config = new PropertiesConfiguration("client.properties");
        alignmentFile = config.getString("alignment");
        hostnameFilesPath = config.getString("hostname.file.path");
        maxConnectionsPerHost = config.getInt("http.host.max.connections");
        requestTimeout = config.getInt("http.request.timeout");
        threadPoolSize = config.getInt("http.thread.pool.size");
    }

    public static void main(String[] args) throws Exception {
        // TODO: There should be a 'BaseClient' and then implementations of 'run', or task etc.
        // e.g. simple run, optim one variable (mu), optim multiple variables (pi + kappa)
        // TODO: This needs to be a config item - should use same JCommander option configuration!
        GeneticCode.initialise(GeneticCode.VERTEBRATE_MITOCHONDRIAL_CODE);
        ClientHopspack c = new ClientHopspack();
        c.loadConfiguration();
        c.run(args);
    }

    /**
     * Little utility class that can be used to construct a SortedMap that sorts
     * keys by values, descending. It is type-safe.
     * <p/>
     * e.g.
     * <p/>
     * Map<K, V> unsorted = new HashMap<K,V>();
     * SortedMap<K, V> sorted = new TreeMap<K, V>(new ValueComparer<K, V>(unsorted));
     * sorted.putAll(unsorted);
     *
     * @param <K>
     * @param <V>
     */
    class ValueComparer<K, V extends Comparable<V>> implements Comparator<K> {
        private final Map<K, V> map;

        public ValueComparer(Map<K, V> map) {
            this.map = map;
        }

        // Compare two values in a map (in descending order)
        public int compare(K key1, K key2) {
            V value1 = this.map.get(key1);
            V value2 = this.map.get(key2);
            int c = value2.compareTo(value1);
            if (c != 0) {
                return c;
            }
            Integer hashCode1 = key1.hashCode();
            Integer hashCode2 = key2.hashCode();
            return hashCode1.compareTo(hashCode2);
        }
    }

}