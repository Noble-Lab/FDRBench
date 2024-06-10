package main.java.FDR;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;


public class CallCrux {

    public static HashMap<String, ArrayList<String>> load_crux_peptides(String crux_peptides_file, int n_random) throws IOException {
        BufferedReader reader = new BufferedReader(new FileReader(crux_peptides_file));
        String line = reader.readLine().trim();
        // target	decoy(s)	mass	proteins
        String [] header = line.split("\t");
        HashMap<String,Integer> col_name2index = new HashMap<>();
        for (int i = 0; i < header.length; i++) {
            col_name2index.put(header[i], i);
        }
        HashMap<String, HashSet<String>> target_peptide2random_peptides = new HashMap<>();
        while ((line = reader.readLine()) != null) {
            String [] fields = line.split("\t");
            String target = fields[col_name2index.get("target")];
            String decoys = fields[col_name2index.get("decoy(s)")];
            // String mass = fields[col_name2index.get("mass")];
            // String proteins = fields[col_name2index.get("proteins")];
            String [] rnd_peptides = decoys.split(",");
            if(rnd_peptides.length< 1 || rnd_peptides[0].length()< target.length()){
                // ,,
                System.out.println("Invalid decoys: " + line);
            }else{
                // GGGGAGK,GGGGAGK,GGAGGGK
                if(!target_peptide2random_peptides.containsKey(target)){
                    target_peptide2random_peptides.put(target, new HashSet<>());
                }
                for(String p: rnd_peptides){
                    if(p.length() < target.length()){
                        System.out.println("Invalid decoy: " + p);
                    }else{
                        target_peptide2random_peptides.get(target).add(p);
                    }
                }
            }
        }
        reader.close();

        HashMap<String, ArrayList<String>> target2decoys = new HashMap<>();
        for(String target: target_peptide2random_peptides.keySet()){
            if(target_peptide2random_peptides.get(target).size() < n_random){
                System.out.println("Not enough random peptides for target: " + target);
            }else{
                if(target_peptide2random_peptides.get(target).size()==n_random){
                    target2decoys.put(target, new ArrayList<>(target_peptide2random_peptides.get(target)));
                }else{
                    for(String p: target_peptide2random_peptides.get(target)){
                        if(target2decoys.containsKey(target)){
                            target2decoys.get(target).add(p);
                            if(target2decoys.get(target).size() == n_random){
                                break;
                            }
                        }else{
                            target2decoys.put(target, new ArrayList<>());
                            target2decoys.get(target).add(p);
                        }
                    }
                }
            }

        }

        return target2decoys;

    }
}
