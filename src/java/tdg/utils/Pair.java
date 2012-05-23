package tdg.utils;

/**
 * A (handy) generic Pair tuple class
 */
public class Pair<F, S> {
    public final F first; //first member of pair
    public final S second; //second member of pair

    public Pair(F first, S second) {
        this.first = first;
        this.second = second;
    }
}

