package main.java.FDR;

import com.compomics.util.experiment.identification.protein_sequences.digestion.ExtendedPeptide;
import com.compomics.util.parameters.identification.search.DigestionParameters;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.concurrent.ConcurrentHashMap;

import static main.java.FDR.FDREval.digest_protein;

public final class DigestProteinWorker4FDR implements Runnable{

    private ConcurrentHashMap<String, HashSet<String>> pro2pepSeq;
    private String proteinID = "";
    private String proteinSequence = "";
    private final DigestionParameters digestionParameters;


    public DigestProteinWorker4FDR(String proID, String proSeq, DigestionParameters digestionParameters, ConcurrentHashMap<String, HashSet<String>> pro2pepMap){
        this.proteinID = proID;
        this.proteinSequence = proSeq;
        this.digestionParameters = digestionParameters;
        this.pro2pepSeq = pro2pepMap;
    }


    @Override
    public void run() {
        ArrayList<ExtendedPeptide> peptides = digest_protein(this.proteinSequence,this.digestionParameters);
        for(ExtendedPeptide pep: peptides){
            if(!pro2pepSeq.containsKey(this.proteinID)){
                pro2pepSeq.put(this.proteinID,new HashSet<>());
            }
            pro2pepSeq.get(this.proteinID).add(pep.peptide.getSequence());
        }
    }
}
