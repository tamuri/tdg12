package tdg.models.parameters;

/**
 * @author Asif Tamuri
 * @version $Id: Parameter.java 128 2010-08-10 12:19:12Z tamuri $
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
