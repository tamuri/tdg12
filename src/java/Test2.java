import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.primitives.Doubles;
import tdg.cli.GeneticCodeConverter;
import tdg.utils.GeneticCode;

import javax.annotation.Generated;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: atamuri
 * Date: 09/03/2011
 * Time: 17:17
 * To change this template use File | Settings | File Templates.
 */
public class Test2 {
    public static void main(String[] args) {
        System.out.printf("%s", Double.parseDouble("-Inf".replaceAll("-Inf", "-Infinity")));

        /*
        Test2 t = new Test2();
        //t.test();
        t.test2();
        */
    }

    private void test2() {
        GeneticCode.initialise(GeneticCode.VERTEBRATE_MITOCHONDRIAL_CODE);

        GeneticCode gc = GeneticCode.getInstance();

        // create a neutral codon pis...equal fitnesses & base nucleotide compositions given by:
        double[] pi = new double[]{ 1.85704375000000e-01,2.73384375000000e-01,4.78586875000000e-01,6.23243750000000e-02 };

        double[] Pi = new double[GeneticCode.CODON_STATES];

        for (int i = 0; i < GeneticCode.CODON_STATES; i++) {
            char[] nucs = gc.getNucleotidesFromCodonIndex(i);
            double prod_nuc = pi[gc.getNucleotideIndexByChar(nucs[0])] * pi[gc.getNucleotideIndexByChar(nucs[1])] * pi[gc.getNucleotideIndexByChar(nucs[2])];
            Pi[i] = prod_nuc; // fitness = 0, exp(F) = 1
        }

        double[] refPi = new double[GeneticCode.AMINO_ACID_STATES];
        double sum = 0;
        for (int i = 0; i < GeneticCode.CODON_STATES; i++) {
            if (!gc.isUnknownAminoAcidState(gc.getAminoAcidIndexFromCodonIndex(i))) {
                refPi[gc.getAminoAcidIndexFromCodonIndex(i)] += Pi[i];
                sum += Pi[i];
            }
        }


        // normalise
        for (int i = 0; i < refPi.length; i++) {
            refPi[i] /= sum;
        }


        System.out.printf("%s\n", Doubles.join(", ", refPi));

//        GeneticCode gc = GeneticCode.getInstance();
 //       for (int i = 0; i < GeneticCode.CODON_STATES; i++) {
  //          System.out.printf("%s %s XXX %s\n", i+1, gc.getAminoAcidCharByIndex(gc.getAminoAcidIndexFromCodonIndex(i)), gc.getCodonTLA(i));
   //     }
    }

    private void test() {

/*
        List<Integer> l = Lists.newArrayList(1, 2, 3, 4, 5, 6);
        Queue<Integer> q = new LinkedList<Integer>(l);
        Iterator<Integer> i = q.iterator();

        Iterator j = Iterators.cycle(i);
             int pos = 0;
        while (true) {
            System.out.printf("%s\n", i.next());
            pos++;
            if (pos % 3 == 0) Iterators.addAll(Lists.<Integer>newArrayList(1), i);
            if (pos > 100) {
                break;
            }
        }
*/


        Map<Integer, Integer> m = Maps.newConcurrentMap();
        m.put(1, 0);
        m.put(2, 0);
        m.put(3, 0);
        m.put(4, 0);
        m.put(5, 0);
        m.put(6, 0);

        Iterator<Integer> j = m.keySet().iterator();
        Iterator<Integer> i = Iterators.cycle(m.keySet());
        int pos = 0;

        int offset = 0;

        while (true) {
            System.out.printf("%s\n", i.next());
            pos++;
            if (pos == 10) m.put(7 + offset++, 0);
            if (pos > 700) break;
        }

        System.exit(0);
        MyObject o1 = new MyObject();
        MyObject2 o2 = new MyObject2();
        o2.setName("hello");
        o1.setMyObject2(o2);

        MyObject o3 = null;
        try {
            o3 = (MyObject) o1.clone();
        } catch (CloneNotSupportedException e) {
            e.printStackTrace();
        }

        o1.printHash();
        if (o3 != null) o3.printHash();
    }

    class MyObject implements Cloneable {
        private int test1 = 1;
        private double test2 = 2.3;
        private String test3 = "four.five";
        private MyObject2 test4;

        public void setMyObject2(MyObject2 o) {
            test4 = o;
        }

        public void printHash() {
            System.out.printf("%s = %s\n", test4.hashCode(), test4.name);

        }

        @Override
        public Object clone() throws CloneNotSupportedException {
            MyObject c = (MyObject) super.clone();
            MyObject2 o2 = (MyObject2) test4.clone();
            //o2.setName("world");
            c.setMyObject2(o2);
            return c;
        }
    }

    class MyObject2 implements Cloneable {
        private String name;

        public void setName(String s) {
            name = s;
        }

        @Override
        protected Object clone() throws CloneNotSupportedException {
            MyObject2 o = (MyObject2) super.clone();
            return o;
        }
    }

}
