package main.java.FDR;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;

public class PEntrapment {


    // For targets
    // A unique ID for each match
    // Precursor level: peptide + modification +charge
    // Peptide level: peptide;
    // Protein level: protein
    public String id;
    public boolean is_target_present = false;
    public PMatch target_hit = new PMatch();

    // For entrapment
    public ArrayList<PMatch> entrapment_hits = new ArrayList<>();

    public static ArrayList<String> get_entrapment_peptides_by_target_peptide(String target_peptide){
        ArrayList<String> entrapment = new ArrayList<>();
        return entrapment;
    }


    public static ArrayList<String> get_entrapment_proteins_by_target_protein(String target_protein){
        ArrayList<String> entrapment = new ArrayList<>();
        return entrapment;
    }

    public static String get_target_peptide_by_entrapment_peptide(String entrapment_peptide, int k_fold, HashMap<String, String> entrapment2target){
        return entrapment2target.get(entrapment_peptide);
    }

    public static String get_target_protein_by_entrapment_protein(String entrapment_protein, int k_fold, String prefix, String sep) {
        String[] proteins = entrapment_protein.split(sep);
        if (k_fold == 1) {
            if (proteins[0].endsWith(prefix)) {
                // _p_target
                return proteins[0].replaceAll(prefix + "$", "");
            } else {
                System.err.println("Error: " + entrapment_protein);
                //System.exit(1);
                return "";
            }
        } else {
            if (proteins[0].endsWith(prefix)) {
                // _p_target
                return proteins[0].replaceAll("_\\d+" + prefix + "$", "");
            } else {
                System.err.println("Error: " + entrapment_protein);
                //System.exit(1);
                return "";
            }
        }
    }

    public static String format_pg(String protein, String sep, String entrapment_prefix, FDREval.PickOneProtein pick_one_protein_method){
        String [] d = protein.split(sep);
        if(d.length==1){
            return d[0];
        }else{
            if(protein.contains(entrapment_prefix)){
                int n_entrapment = 0;
                for(int i=0;i<d.length;i++){
                    if(d[i].endsWith(entrapment_prefix)){
                        n_entrapment++;
                    }
                }
                if(n_entrapment==d.length) {
                    // all the proteins in the group are entrapment proteins
                    if(pick_one_protein_method == FDREval.PickOneProtein.first){
                        return d[0];
                    }else if(pick_one_protein_method == FDREval.PickOneProtein.last) {
                        return d[d.length - 1];
                    }else if(pick_one_protein_method == FDREval.PickOneProtein.random) {
                        // generate a random number between 0 and d.length-1
                        Random rand = new Random();
                        int randomNumber = rand.nextInt(d.length);
                        return d[randomNumber];
                    }else{
                        return protein;
                    }

                }else{
                    // not all the proteins in the group are entrapment proteins
                    ArrayList<String> target_proteins = new ArrayList<>();
                    for(int i=0;i<d.length;i++){
                        if(!d[i].endsWith(entrapment_prefix)){
                            target_proteins.add(d[i]);
                            //System.out.println(pro+"\t"+protein);
                        }
                    }
                    if(pick_one_protein_method == FDREval.PickOneProtein.first){
                        return target_proteins.get(0);
                    }else if(pick_one_protein_method == FDREval.PickOneProtein.last) {
                        return target_proteins.get(target_proteins.size()-1);
                    }else if(pick_one_protein_method == FDREval.PickOneProtein.random) {
                        // generate a random number between 0 and d.length-1
                        Random rand = new Random();
                        int randomNumber = rand.nextInt(target_proteins.size());
                        return target_proteins.get(randomNumber);
                    }else{
                        return protein;
                    }
                }
            }else{
                // none of the groups in the group is an entrapment protein
                if(pick_one_protein_method == FDREval.PickOneProtein.first){
                    return d[0];
                }else if(pick_one_protein_method == FDREval.PickOneProtein.last) {
                    return d[d.length - 1];
                }else if(pick_one_protein_method == FDREval.PickOneProtein.random) {
                    // generate a random number between 0 and d.length-1
                    Random rand = new Random();
                    int randomNumber = rand.nextInt(d.length);
                    return d[randomNumber];
                }else{
                    return protein;
                }
            }
        }

    }

}
