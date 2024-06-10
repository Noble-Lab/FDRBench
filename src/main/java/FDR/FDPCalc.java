package main.java.FDR;

import com.google.common.collect.Sets;
import tech.tablesaw.api.Table;

import java.util.HashMap;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;

public final class FDPCalc implements Runnable {
    public static Table psm_table;

    public static HashMap<Integer,PMatch> psm_map;
    public static String target_col_name;
    public static FDRType fdp_level;
    public static HashMap<String,String> target2label = new HashMap<>(100000);
    public static PeptidePair peptidePair;

    private int i=0;

    public ConcurrentHashMap<Integer, HashMap<String, Double>> res = new ConcurrentHashMap<>();

    public FDPCalc(int i, ConcurrentHashMap<Integer, HashMap<String, Double>> res) {
        this.i = i;
        this.res = res;
    }

    @Override
    public void run() {
        if(i % 2000 == 0){
            System.out.println("Processing "+i+" matches");
        }
        // for method 1
        // entrapment targets
        int n_entrapment_targets = 0;
        // targets from target protein database
        int n_targets = 0;
        // protein or peptide level
        //HashSet<String> all_hits = new HashSet<>();
        Set<String> all_hits = Sets.newHashSetWithExpectedSize(1000000);
        // precursor level
        //HashSet<String> all_precursor_hits = new HashSet<>();
        Set<String> all_precursor_hits = Sets.newHashSetWithExpectedSize(1000000);
        for(int k=0;k<=i;k++){
            // Row row = psm_table.row(k);
            PMatch pmatch = psm_map.get(k);
            // protein or peptide column name
            String target_name;
            if(fdp_level.equals(FDRType.protein)){
                target_name = pmatch.protein;
            }else{
                target_name = pmatch.peptide;
            }
            String label = "target";
            if (target2label.containsKey(target_name) || fdp_level.equals(FDRType.protein)) {
                // for peptide/psm/precursor, we need to retrieve the label from database peptide list file.
                label = target2label.get(target_name);
                if (fdp_level.equals(FDRType.protein)) {
                    // if it's protein level, only the first protein is considered to retrieve the label
                    target_name = target_name.split(";")[0];
                    if (target_name.endsWith("p_target")) {
                        label = "p_target";
                    } else {
                        label = "target";
                    }
                }

            } else {
                System.err.println("Error3: invalid: " + target_name);
            }
            // peptide or protein level
            all_hits.add(target_name);
            // precursor level
            if(fdp_level.equals(FDRType.precursor)){
                String mod_seq = pmatch.mod_peptide;
                int charge = pmatch.charge;
                all_precursor_hits.add(target_name+"|"+mod_seq.length()+"|"+charge);
            }
            // for method 1: it doesn't need pair information.
            if (label.equals("target")) {
                n_targets = n_targets + 1;
            } else if (label.equals("p_target")) {
                n_entrapment_targets = n_entrapment_targets + 1;
            } else {
                System.err.println("Error1: invalid target type:" + label);
                System.exit(1);
            }

        }

        // for method 1B
        //HashSet<String> cur_hits = new HashSet<>();
        Set<String> cur_hits = Sets.newHashSetWithExpectedSize(1000000);
        // precursor level
        // HashSet<String> cur_precursor_hits = new HashSet<>();
        Set<String> cur_precursor_hits = Sets.newHashSetWithExpectedSize(1000000);
        int n_p_t_s = 0; // N_{p>t>s} = all PTs that score > their paired target and both are > s
        int n_p_s_t = 0; // N_{p>s>t} = all PTs that score > s > the paired target
        for(int k=i;k>=0;k--) {
            // from low confident to high confident
            // Row row = psm_table.row(k);
            PMatch pmatch = psm_map.get(k);
            // protein or peptide column name
            String target_name;
            if(fdp_level.equals(FDRType.protein)){
                target_name = pmatch.protein;
            }else{
                target_name = pmatch.peptide;
            }
            if (target2label.containsKey(target_name) || fdp_level.equals(FDRType.protein)) {
                // for peptide/psm/precursor, we need to retrieve the label from database peptide list file.
                String label = target2label.get(target_name);
                if (fdp_level.equals(FDRType.protein)) {
                    // if it's protein level, only the first protein is considered to retrieve the label
                    target_name = target_name.split(";")[0];
                    if (target_name.endsWith("p_target")) {
                        label = "p_target";
                    } else {
                        label = "target";
                    }
                    // label = target2label.get(target_name);

                }
                cur_hits.add(target_name);
                // precursor level
                if(fdp_level.equals(FDRType.precursor)){
                    String mod_seq = pmatch.mod_peptide;
                    int charge = pmatch.charge;
                    cur_precursor_hits.add(target_name+"|"+mod_seq.length()+"|"+charge);
                }

                // for method 1B: it needs pair information.
                if (fdp_level.equals(FDRType.protein)) {
                    // for protein level.
                    /*
                     * N_{p>t>s} = all PTs that score > their paired target and both are > s
                     * N_{p>s>t} = all PTs that score > s > the paired target
                     */
                    if (label.equals("p_target")) {
                        // paired target
                        String paired_target = target_name.replaceAll("_p_target", "");
                        if (cur_hits.contains(paired_target)) {
                            // we found its target before
                            // p > t > s
                            n_p_t_s = n_p_t_s + 1;
                        } else if(!all_hits.contains(paired_target)) {
                            // we don't find its target
                            // p > s > t
                            // N_{p>s>t} = all PTs that score > s > the paired target
                            n_p_s_t = n_p_s_t + 1;
                        }
                    }
                } else if (fdp_level.equals(FDRType.precursor)) {
                    // for method 1b
                    // for precursor level.
                    /*
                     * N_{p>t>s} = all PTs that score > their paired target and both are > s
                     * N_{p>s>t} = all PTs that score > s > the paired target
                     */
                    // need modified peptide sequence
                    if (label.equals("p_target")) {
                        // target_name: peptide
                        String paired_peptide = peptidePair.get_paired_peptide(target_name);
                        String mod_seq = pmatch.mod_peptide;
                        int charge = pmatch.charge;
                        // precursor string: paired_peptide|length(mod_seq)|charge
                        String paired_precursor = paired_peptide+"|"+mod_seq.length()+"|"+charge;
                        if (cur_precursor_hits.contains(paired_precursor)) {
                            // we found its target before
                            // p > t > s
                            n_p_t_s = n_p_t_s + 1;
                        } else if(!all_precursor_hits.contains(paired_precursor)){
                            // we don't find its target before
                            // p > s > t
                            // N_{p>s>t} = all PTs that score > s > the paired target
                            n_p_s_t = n_p_s_t + 1;
                        }
                    }
                } else if (fdp_level.equals(FDRType.peptide)){

                    if (label.equals("p_target")) {
                        // target_name: peptide
                        String paired_peptide = peptidePair.get_paired_peptide(target_name);
                        if (cur_hits.contains(paired_peptide)) {
                            // we found its target before
                            // p > t > s
                            n_p_t_s = n_p_t_s + 1;
                        } else if(!all_hits.contains(paired_peptide)){
                            // we don't find its target before
                            // p > s > t
                            // N_{p>s>t} = all PTs that score > s > the paired target
                            n_p_s_t = n_p_s_t + 1;
                        }
                    }
                }

            } else {
                System.err.println("Error3: invalid: " + target_name);
                //System.exit(1);
            }
        }

        double fdp = 2.0 * n_entrapment_targets / (n_entrapment_targets + n_targets);
        res.put(i, new HashMap<>());
        res.get(i).put("n_entrapment_targets",n_entrapment_targets+0.0);
        res.get(i).put("n_targets",n_targets+0.0);
        res.get(i).put("n_p_t_s",n_p_t_s+0.0);
        res.get(i).put("n_p_s_t",n_p_s_t+0.0);

        res.get(i).put("fdp",fdp);
        double fdp_1b = (n_entrapment_targets + 2.0 * n_p_t_s + n_p_s_t) / (n_entrapment_targets + n_targets);
        res.get(i).put("paired_fdp",fdp_1b);
        double low_bound_fdp = 1.0*n_entrapment_targets / (n_entrapment_targets + n_targets);
        res.get(i).put("low_bound_fdp",low_bound_fdp);
    }

}
