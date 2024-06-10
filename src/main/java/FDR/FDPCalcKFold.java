package main.java.FDR;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ThreadLocalRandom;

import static main.java.FDR.FDREval.debug;

public final class FDPCalcKFold implements Runnable {
    // public static Table psm_table;

    public static HashMap<Integer,PMatch> psm_map;
    // public static String target_col_name;
    // public static FDRType fdp_level;
    // public static HashMap<String,String> target2label = new HashMap<>(100000);
    public static HashMap<String, PEntrapment> target2entrapment = new HashMap<>();
    // public static PeptidePair peptidePair;

    public static int k = 1;

    private int i=0;

    public ConcurrentHashMap<Integer, HashMap<String, Double>> res = new ConcurrentHashMap<>();

    public FDPCalcKFold(int i, ConcurrentHashMap<Integer, HashMap<String, Double>> res) {
        this.i = i;
        this.res = res;
    }

    @Override
    public void run() {
        if(i % 2000 == 0){
            System.out.println("Processing "+i+" matches");
        }
        // entrapment targets
        int n_entrapment_targets = 0;
        // targets from target protein database
        int n_targets = 0;
        for(int k=0;k<=i;k++){
            // Row row = psm_table.row(k);
            PMatch pmatch = psm_map.get(k);
            if(pmatch.is_entrapment){
                n_entrapment_targets = n_entrapment_targets + 1;
            }else{
                n_targets = n_targets + 1;
            }
        }

        // for paired entrapment method
        PMatch pmatch = psm_map.get(i);
        double si = pmatch.rank_score;
        HashMap<Integer,Integer> i2ni = new HashMap<>();
        for(String target_id: target2entrapment.keySet()){
            PEntrapment p = target2entrapment.get(target_id);
            if(p.entrapment_hits.size()>k){
                // sanity check
                System.err.println("The number of entrapment hits is larger than k: "+p.entrapment_hits.size()+" > k="+k);
                System.err.println("Target: "+target_id + ", score:"+si);
                for(PMatch pm:p.entrapment_hits){
                    System.err.println("Entrapment: "+pm.id+"\t"+pm.rank_score);
                }
                System.exit(1);
            }
            int nl = getNl(p, si, k);
            if(nl >= 1){
                if(i2ni.containsKey(nl)){
                    i2ni.put(nl,i2ni.get(nl)+1);
                }else{
                    i2ni.put(nl,1);
                }
                //System.out.println(i+"\t"+target_id+"\t"+si+"\t"+nl+"\t"+k+"\t"+i2ni.get(nl));
            }
            int nk = getNk(p, si, k);
            if(nk >= 1){
                if(i2ni.containsKey(k+1)){
                    i2ni.put(k+1,i2ni.get(k+1)+1);
                }else{
                    i2ni.put(k+1,1);
                }
                //System.out.println(i+"\t"+target_id+"\t"+si+"\t"+nk+"\tk+1\t"+i2ni.get(k+1));
            }
        }

        //res.put(i, new HashMap<>());
        int vt = 0;
        if(!i2ni.isEmpty()) {
            ArrayList<String> info = new ArrayList<>();
            for (int l : i2ni.keySet()) {
                vt = vt + l * i2ni.get(l);
                info.add(l + ":" + i2ni.get(l));
            }
            //res.get(i).put("info", String.join(",", info));
            if(debug) {
                System.out.println("info:"+String.join(",", info));
            }
        }

        // for combined entrapment method
        double fdp = 1.0* n_entrapment_targets * (1+1.0/k) / (n_entrapment_targets + n_targets);

        res.get(i).put("n_entrapment_targets",n_entrapment_targets+0.0);
        res.get(i).put("n_targets",n_targets+0.0);
        if(k==1) {
            if(i2ni.containsKey(1)) {
                res.get(i).put("n_p_s_t", i2ni.get(1) + 0.0);
            }else{
                res.get(i).put("n_p_s_t", 0.0);
            }

            if(i2ni.containsKey(2)) {
                res.get(i).put("n_p_t_s", i2ni.get(2) + 0.0);
            }else{
                res.get(i).put("n_p_t_s", 0.0);
            }
            res.get(i).put("vt",vt+0.0);
        }else{
            res.get(i).put("vt",vt+0.0);

        }

        res.get(i).put("fdp",fdp);
        // for paired entrapment method
        // double fdp_1b = (n_entrapment_targets + 2.0 * n_p_t_s + n_p_s_t) / (n_entrapment_targets + n_targets);
        double fdp_1b = 1.0* (n_entrapment_targets + vt) / (n_entrapment_targets + n_targets);
        if(debug && k==1){
            System.out.println(n_entrapment_targets+"\t"+n_targets+"\t"+res.get(i).get("n_p_s_t")+"\t"+res.get(i).get("n_p_t_s")+"\t"+vt+"\t"+fdp_1b);
        }
        double low_bound_fdp = 1.0*n_entrapment_targets / (n_entrapment_targets + n_targets);
        res.get(i).put("low_bound_fdp",low_bound_fdp);
        res.get(i).put("paired_fdp",fdp_1b);
    }


    private static int getNl(PEntrapment p, double s, int k) {
        double ti = p.target_hit.rank_score;
        int ni = 0;
        for(PMatch e: p.entrapment_hits){
            double ei = e.rank_score;
            // the entrapment hits with confidence score >= the current confidence score cutoff
            // ei >= s >= ti
            if (ei >= s && s > ti) {
                ni = ni + 1;
            }else{
                ni = 0;
                break;
            }
        }
        if(Double.NEGATIVE_INFINITY == ti){
            // if t is -Inf, we need to consider the case that some of the entrapment sequences have no match.
            int n_entrapment_no_hit = k - p.entrapment_hits.size();
            if(n_entrapment_no_hit>=1) {
                // Random rand = new Random();
                int randomNumber = ThreadLocalRandom.current().nextInt(n_entrapment_no_hit + 1); // This will generate a random number between 0 (inclusive) and 6 (exclusive)
                // target ranks the last
                if (randomNumber == n_entrapment_no_hit) {
                    // the target ranks the last
                } else {
                    // don't count this target
                    ni = 0;
                }
            }
        }
        return ni;
    }

    private static int getNk(PEntrapment p, double s,int k) {
        double ti = p.target_hit.rank_score;
        int ni = 0;
        // all entrapment sequences are present in the list
        if(p.entrapment_hits.size()==k) {
            for (PMatch e : p.entrapment_hits) {
                double ei = e.rank_score;
                // the entrapment hits with confidence score >= the current confidence score cutoff
                // ei >= t >= s
                if (ei >= ti && ti >= s) {
                    ni = ni + 1;
                } else {
                    ni = 0;
                    break;
                }
            }
        }
        return ni;
    }

    private static int get_n_entrapment(PEntrapment p, double si) {
        // The number of entrapment hits confidence score >= the current confidence score cutoff
        int n_entrapment_hits = 0;
        for(PMatch e: p.entrapment_hits){
            double ei = e.rank_score;
            if (ei >= si) {
                n_entrapment_hits = n_entrapment_hits + 1;
            }
        }
        return n_entrapment_hits;
    }

}
