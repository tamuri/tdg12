package tdg;

import com.google.common.collect.Maps;
import tdg.model.Fitness;

import java.io.Serializable;
import java.util.Iterator;
import java.util.Map;

/**
 * User: atamuri
 * Date: 27/02/2013 15:43
 */
public class FitnessStore implements Serializable, Iterable<Fitness> {
    private static final long serialVersionUID = -2630793049038547404L;

    Map<Integer, Fitness> store;

    public FitnessStore() {
        store = Maps.newTreeMap();
    }

    public void setFitness(int site, Fitness f) {
        store.put(site, f);
    }

    public Fitness getFitness(int site) {
        if (!store.containsKey(site)) throw new RuntimeException("No fitness stored for site " + site);
        return store.get(site);
    }


    @Override
    public Iterator<Fitness> iterator() {
        return store.values().iterator();
    }
}
