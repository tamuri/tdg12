package tdg.utils;

import com.google.common.collect.Maps;

import java.util.Map;
import java.util.concurrent.atomic.AtomicLong;

/**
 * @author Asif Tamuri
 * @version $Id: CodeTimer.java 148 2010-08-13 16:35:08Z tamuri $
 */
public class CodeTimer {
    private static final Map<String, AtomicLong> times = Maps.newConcurrentMap();

    public final static void store(String key, long startTime) {
        if (!times.containsKey(key)) {
            times.put(key, new AtomicLong());
        }
        times.get(key).getAndAdd(System.currentTimeMillis() - startTime);
    }

    public final static long start() {
        return System.currentTimeMillis();
    }

    public final static void printAll() {
        System.out.println("\n-------------------\ntdg.utils.CodeTimer\n-------------------");
        for (Map.Entry e : times.entrySet()) {
            System.out.printf("%s = %sms\n", e.getKey(), e.getValue());
        }
    }
}
