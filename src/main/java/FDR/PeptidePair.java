package main.java.FDR;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import static main.java.FDR.FDREval.get_column_names;

public class PeptidePair {

    /**
     * peptide -> peptide_pair_index
     */
    HashMap<String,String> pep2peptide_pair_index = new HashMap<>();
    HashMap<String,HashMap<String,String>> peptide_pair_index2peptides = new HashMap<>();

    /**
     * For k-fold
     * peptide_pair_index -> type -> peptides
     */
    HashMap<String,HashMap<String, ArrayList<String>>> peptide_pair_index2type2peptides = new HashMap<>();

    public void load(String pep_file, HashSet<String> peptides){
        HashMap<String,Integer> col2index = get_column_names(pep_file,"\t");
        try {
            BufferedReader br = new BufferedReader(new FileReader(pep_file));
            // sequence decoy proteins peptide_type peptide_pair_index
            // decoy: Yes, No
            // peptide_type: target, p_target, decoy, p_decoy
            String head = br.readLine().trim();
            String line;
            while((line=br.readLine())!=null){
                String []d=line.trim().split("\t");
                String peptide_seq = d[col2index.get("sequence")];
                if(peptides.contains(peptide_seq)){
                    String peptide_type = d[col2index.get("peptide_type")];
                    String peptide_pair_index = d[col2index.get("peptide_pair_index")];
                    pep2peptide_pair_index.put(peptide_seq,peptide_pair_index);
                    if(!peptide_pair_index2peptides.containsKey(peptide_pair_index)){
                        peptide_pair_index2peptides.put(peptide_pair_index,new HashMap<>());
                    }
                    peptide_pair_index2peptides.get(peptide_pair_index).put(peptide_type,peptide_seq);
                }

            }
            br.close();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    public void load(String pep_file, HashSet<String> peptides, int k_fold){
        HashMap<String,Integer> col2index = get_column_names(pep_file,"\t");
        try {
            BufferedReader br = new BufferedReader(new FileReader(pep_file));
            // sequence decoy proteins peptide_type peptide_pair_index
            // decoy: Yes, No
            // peptide_type: target, p_target, decoy, p_decoy
            String head = br.readLine().trim();
            String line;
            HashSet<String> all_peptide_index = new HashSet<>();
            while((line=br.readLine())!=null){
                String []d=line.trim().split("\t");
                String peptide_seq = d[col2index.get("sequence")];
                String decoy_type = d[col2index.get("decoy")];
                if(decoy_type.startsWith("Yes")){
                    continue;
                }
                if(peptides.contains(peptide_seq)){
                    String peptide_pair_index = d[col2index.get("peptide_pair_index")];
                    all_peptide_index.add(peptide_pair_index);
                }

            }
            br.close();

            br = new BufferedReader(new FileReader(pep_file));
            head = br.readLine().trim();
            while((line=br.readLine())!=null){
                String []d=line.trim().split("\t");
                String peptide_pair_index = d[col2index.get("peptide_pair_index")];
                String decoy_type = d[col2index.get("decoy")];
                if(decoy_type.startsWith("Yes")){
                    continue;
                }
                if(all_peptide_index.contains(peptide_pair_index)){
                    String peptide_seq = d[col2index.get("sequence")];
                    String peptide_type = d[col2index.get("peptide_type")];
                    pep2peptide_pair_index.put(peptide_seq,peptide_pair_index);
                    if(!peptide_pair_index2type2peptides.containsKey(peptide_pair_index)){
                        peptide_pair_index2type2peptides.put(peptide_pair_index,new HashMap<>());
                    }
                    if(peptide_pair_index2type2peptides.get(peptide_pair_index).containsKey(peptide_type)){
                        peptide_pair_index2type2peptides.get(peptide_pair_index).get(peptide_type).add(peptide_seq);
                    }else{
                        peptide_pair_index2type2peptides.get(peptide_pair_index).put(peptide_type,new ArrayList<>());
                        peptide_pair_index2type2peptides.get(peptide_pair_index).get(peptide_type).add(peptide_seq);
                    }
                }
            }
            br.close();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    public String get_paired_peptide(String peptide){
        String paired_peptide = "-";
        if(this.pep2peptide_pair_index.containsKey(peptide)){
            String peptide_pair_index = this.pep2peptide_pair_index.get(peptide);
            for(String type: this.peptide_pair_index2peptides.get(peptide_pair_index).keySet()){
                String pep = this.peptide_pair_index2peptides.get(peptide_pair_index).get(type);
                if(!pep.equalsIgnoreCase(peptide)){
                    paired_peptide = pep;
                    break;
                }
            }
        }else{
            // no paired peptide
            paired_peptide = "-";
        }
        return paired_peptide;
    }


    public String get_paired_peptide(String peptide, int k_fold){
        String paired_peptide = "-";
        if(this.pep2peptide_pair_index.containsKey(peptide)){
            String peptide_pair_index = this.pep2peptide_pair_index.get(peptide);
            for(String type: this.peptide_pair_index2type2peptides.get(peptide_pair_index).keySet()){
                String pep = this.peptide_pair_index2peptides.get(peptide_pair_index).get(type);
                if(!pep.equalsIgnoreCase(peptide)){
                    paired_peptide = pep;
                    break;
                }
            }
        }else{
            // no paired peptide
            paired_peptide = "-";
        }
        return paired_peptide;
    }

    public String get_paired_target_peptide(String peptide){
        String paired_peptide = "-";
        if(this.pep2peptide_pair_index.containsKey(peptide)){
            String peptide_pair_index = this.pep2peptide_pair_index.get(peptide);
            paired_peptide = this.peptide_pair_index2type2peptides.get(peptide_pair_index).get("target").get(0);
        }else{
            // no paired peptide
            System.err.println("No target peptide exists for "+peptide);
            System.exit(1);
        }
        return paired_peptide;
    }
}
