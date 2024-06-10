package main.java.FDR;

import java.util.HashMap;
import java.util.concurrent.ConcurrentHashMap;

public final class FDPCalcCE implements Runnable {
    public static HashMap<Integer,PMatch> psm_map;
    private int i=0;

    public ConcurrentHashMap<Integer, HashMap<String, Double>> res = new ConcurrentHashMap<>();

    public static double r = 1.0;

    public FDPCalcCE(int i, ConcurrentHashMap<Integer, HashMap<String, Double>> res) {
        this.i = i;
        this.res = res;
    }

    @Override
    public void run() {
        // for method 1
        // entrapment targets
        int n_entrapment_targets = 0;
        // targets from target protein database
        int n_targets = 0;
        for(int k=0;k<=i;k++){
            // Row row = psm_table.row(k);
            PMatch pmatch = psm_map.get(k);
            // for method 1: it doesn't need pair information.
            if (pmatch.is_entrapment) {
                n_entrapment_targets = n_entrapment_targets + 1;
            } else {
                n_targets = n_targets + 1;
            }

        }
        double fdp = n_entrapment_targets * (1+1.0/r) / (n_entrapment_targets + n_targets);
        res.put(i, new HashMap<>());
        res.get(i).put("n_entrapment_targets",n_entrapment_targets+0.0);
        res.get(i).put("n_targets",n_targets+0.0);
        res.get(i).put("fdp",fdp);
        double low_bound_fdp = 1.0*n_entrapment_targets / (n_entrapment_targets + n_targets);
        res.get(i).put("low_bound_fdp",low_bound_fdp);
    }

}
