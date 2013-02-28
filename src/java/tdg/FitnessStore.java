package tdg;

import com.google.common.collect.Lists;
import tdg.model.Fitness;

import java.util.List;

/**
 * User: atamuri
 * Date: 27/02/2013 15:43
 */
public class FitnessStore {
    List<Fitness> fitnessStore;

    public FitnessStore(int size) {
        fitnessStore = Lists.newArrayList();
        for (int i = 0; i < size; i++) fitnessStore.add(null);
    }

    public void setFitness(int site, Fitness f) {
        fitnessStore.set(site - 1, f);
    }

    public Fitness getFitness(int site) {
        if (fitnessStore.get(site - 1) == null) {
            throw new RuntimeException("No fitness stored for site " + site);
        }

        return fitnessStore.get(site - 1);
    }
}
