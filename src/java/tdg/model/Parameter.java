package tdg.model;

/**
 * @author Asif Tamuri (atamuri@nimr.mrc.ac.uk)
 */
public abstract class Parameter {
    private Object value;
    private boolean optimiseValue;

    public Object get() {
        return this.value;
    }

    public void set(Object value) {
        this.value = value;
    }

    public boolean isOptimiseValue() {
        return optimiseValue;
    }

    public void setOptimiseValue(boolean optimiseValue) {
        this.optimiseValue = optimiseValue;
    }
}
