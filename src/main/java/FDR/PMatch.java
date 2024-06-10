package main.java.FDR;

public class PMatch {

    public String peptide;
    public String mod_peptide;
    public int charge;
    public String protein;
    public String modification;

    /**
     * A score used to rank the hits.
     */
    public double rank_score = Double.NEGATIVE_INFINITY;
    public boolean score_higher_is_better = true;

    // A unique ID for each match
    // Precursor level: peptide + modification +charge
    // Peptide level: peptide;
    // Protein level: protein
    public String id;
    public boolean is_entrapment = false;
}
