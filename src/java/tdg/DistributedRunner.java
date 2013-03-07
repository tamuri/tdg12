package tdg;

import com.caucho.hessian.client.HessianProxyFactory;
import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;
import com.google.common.primitives.Ints;
import pal.alignment.Alignment;
import pal.tree.Node;
import pal.tree.Tree;
import tdg.distributed.ServiceAPI;
import tdg.model.TDGGlobals;
import tdg.utils.PhyloUtils;
import tdg.utils.ValueComparer;

import java.net.MalformedURLException;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

/**
 * User: atamuri
 * Date: 05/03/2013 19:17
 */
public class DistributedRunner extends AbstractRunner {
    private List<String> slaves;
    private List<ServiceAPI> slaveService = Lists.newArrayList();
    private List<List<Integer>> slaveSites = Lists.newArrayList();
    private final ExecutorService threadPool;

    public DistributedRunner(Alignment alignment, List<String> slaves) {
        setAlignment(alignment);
        this.slaves = slaves;

        // get a list of sites ordered by number of observed amino acids, descending

        List<Integer> slaveThreads = Lists.newArrayList();
        for (String slave : slaves) {
            try {
                HessianProxyFactory factory = new HessianProxyFactory();
                this.slaveService.add((ServiceAPI) factory.create(ServiceAPI.class, slave));
            } catch (MalformedURLException e) {
                throw new RuntimeException(e);
            }

            // TODO: Use threads-on-servers to distribute sites
            slaveThreads.add(this.slaveService.get(this.slaveService.size() - 1).getThreads());
        }

        Collection<Integer> orderedSites = getSites(alignment);

        // assign every site a particular slave for these runs
        for (int i = 0; i < this.slaveService.size(); i++) {
            slaveSites.add(Lists.<Integer>newArrayList());
        }

        Iterator<List<Integer>> cycleServerSites = Iterators.cycle(this.slaveSites);
        for (int site : orderedSites) {
            cycleServerSites.next().add(site);
        }

        for (int i = 0; i < this.slaveService.size(); i++) {
            this.slaveService.get(i).setSites(Ints.toArray(this.slaveSites.get(i)));
        }

        this.threadPool = Executors.newFixedThreadPool(slaveService.size());

        System.out.printf("%s\n%s\n", slaveService, slaveSites);

    }

    @Override
    public void runnerSetTree(final Tree tree) {
        List<Future<Void>> futures = Lists.newArrayList();

        for (final ServiceAPI s : slaveService) {
            Future<Void> future = threadPool.submit(new Callable<Void>() {
                @Override
                public Void call() throws Exception {
                    s.setTree(tree);
                    return null;
                }
            });

            futures.add(future);
        }

        getAllResults(futures);
    }

    @Override
    protected void runnerSetFitnessStore(final FitnessStore fitnessStore) {
        List<Future<Void>> futures = Lists.newArrayList();

        for (int i = 0; i < slaveService.size(); i++) {
            final int service_i = i;


            Future<Void> future = threadPool.submit(new Callable<Void>() {
                @Override
                public Void call() throws Exception {
                    final ServiceAPI s = slaveService.get(service_i);

                    final List<Integer> sites = slaveSites.get(service_i);

                    // each slave should just have the fitness for the sites for which it is responsible
                    FitnessStore slaveFitnessStore = new FitnessStore(sites.size());
                    for (int i = 0; i < sites.size(); i++) {
                        slaveFitnessStore.setFitness(i + 1, fitnessStore.getFitness(sites.get(i)));
                    }
                    s.setFitnessStore(slaveFitnessStore);

                    return null;
                }
            });

            futures.add(future);
        }

        getAllResults(futures);
    }

    @Override
    public double runnerGetLogLikelihood(final Tree tree, final FitnessStore fitnessStore, final TDGGlobals globals) {
        List<Future<Double>> futures = Lists.newArrayList();

        for (final ServiceAPI service : slaveService) {
            Future<Double> future = threadPool.submit(new Callable<Double>() {
                @Override
                public Double call() throws Exception {
                    return service.optimiseMutationModel(globals);
                }
            });

            futures.add(future);
        }

        double total = 0;
        for (double x : getAllResults(futures)) total += x;

        return total;
    }

    /*

     */

    @Override public double updateSiteCalculatorTrees(final Tree tree, final TDGGlobals globals, final FitnessStore fitnessStore) {
        List<Future<Double>> futures = Lists.newArrayList();

        for (int i = 0; i < slaveService.size(); i++) {
            final int service_i = i;

            Future<Double> future = threadPool.submit(new Callable<Double>() {
                @Override
                public Double call() throws Exception {
                    final ServiceAPI s = slaveService.get(service_i);
                    final List<Integer> sites = slaveSites.get(service_i);

                    // each slave should just have the fitness for the sites for which it is responsible
                    FitnessStore slaveFitnessStore = new FitnessStore(sites.size());
                    for (int i = 0; i < sites.size(); i++) {
                        slaveFitnessStore.setFitness(i + 1, fitnessStore.getFitness(sites.get(i)));
                    }
                    s.setFitnessStore(slaveFitnessStore);
                    s.setTree(tree);

                    return s.updateLikelihoodCalculators(globals);
                }
            });

            futures.add(future);
        }

        double total = 0;

        for (Double p : getAllResults(futures)) {
            total += p;
        }

        return total;
    }

    @Override public double getLikelihoodSum(final Node node, final double newBranchLength) {

        final List<Future<Double>> futures = Lists.newArrayList();

        for (int i = 0; i < slaveService.size(); i++) {
            final int service_i = i;

            Future<Double> future = threadPool.submit(new Callable<Double>() {
                @Override
                public Double call() throws Exception {
                    final ServiceAPI s = slaveService.get(service_i);
                    return s.getNodeLikelihood(node, newBranchLength);
                }
            });

            futures.add(future);
        }

        double total = 0;
        for (double x : getAllResults(futures)) total += x;

        return total;
    }

    @Override
    protected void updateBranchLength(final Node child, final double branchlength) {
        List<Future<Void>> futures = Lists.newArrayList();
        for (int j = 0; j < slaveService.size(); j++) {
            final int service_i = j;

            Future<Void> future = threadPool.submit(new Callable<Void>() {
                @Override
                public Void call() throws Exception {
                    final ServiceAPI s = slaveService.get(service_i);
                    s.setBranchLength(child, branchlength);
                    return null;
                }
            });

            futures.add(future);
        }

        getAllResults(futures);
    }

    @Override
    public double optimiseFitness(Tree tree, TDGGlobals globals, FitnessStore fitnessStore) {
        return 0;
    }

    @Override
    public void close() {

    }

    private Collection<Integer> getSites(final Alignment alignment) {
        // Get the Q matrix size (= number of codons) for each site
        Map<Integer, Integer> siteCodonCount = new HashMap<Integer, Integer>();
        for (int i = 1; i <= (alignment.getSiteCount() / 3); i++) {
            Map<String, Integer> sitePattern = PhyloUtils.getCodonsAtSite(alignment, i);
            List<Integer> aminoAcidsAtSite = PhyloUtils.getDistinctAminoAcids(sitePattern.values());
            List<Integer> allCodons = PhyloUtils.getCodonsFromAminoAcids(aminoAcidsAtSite);

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
}
