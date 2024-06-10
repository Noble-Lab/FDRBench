package main.java.FDR;


public enum FDRType {
    protein("protein level FDR"),
    peptide("peptide level FDR"),
    precursor("precursor level FDR"),
    psm("PSM level FDR");

    /**
     * The description.
     */
    public final String description;

    /**
     * Constructor.
     *
     * @param description the description of PsmType
     */
    FDRType(String description) {
        this.description = description;
    }
}