package main.java.FDR;

import com.compomics.util.experiment.biology.enzymes.Enzyme;
import com.compomics.util.experiment.biology.enzymes.EnzymeFactory;
import com.compomics.util.parameters.identification.search.DigestionParameters;
import com.compomics.util.pride.CvTerm;
import main.java.util.Cloger;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

import static java.util.stream.Collectors.toList;

public class DBGear {

    public boolean I2L = false;

    public DBGear() {

    }

    public HashSet<String> digest_protein(Enzyme enzyme, String proteinSequence){
        proteinSequence = proteinSequence.toUpperCase();
        proteinSequence = proteinSequence.replaceAll("\\*$", "");
        // needed ?
        if(this.I2L) {
            proteinSequence = proteinSequence.replaceAll("I", "L");
        }
        HashSet<String> peptides = enzyme.digest(proteinSequence, CParameter.maxMissedCleavages, CParameter.minPeptideLength, CParameter.maxPeptideLength);
        if(CParameter.clip_nTerm_M && proteinSequence.startsWith("M")){
            List<String> n_term_peptides = peptides.stream().filter(proteinSequence::startsWith).filter(pep -> pep.length() >= (CParameter.minPeptideLength+1)).map(pep -> pep.substring(1)).collect(toList());
            if(n_term_peptides.size()>=1){
                peptides.addAll(n_term_peptides);
            }
        }
        return peptides;
    }

    public static DigestionParameters getDigestionPreferences(){
        return(getDigestionPreferences("Trypsin",2));
    }


    /**
     * Get the enzyme object based on an enzyme index
     * @param ind enzyme index, 0-based index.
     * @return An enzyme object
     */
    public static Enzyme getEnzymeByIndex(int ind){

        ArrayList<Enzyme> enzymes = new ArrayList<>();

        // 0 non-specific digestion
        Enzyme enzyme = new Enzyme("NoEnzyme");
        String all_aas = "ABCDEFGHIKLMNPQRSTUVWXY";
        for(int i=0;i<all_aas.length();i++){
            enzyme.addAminoAcidBefore(all_aas.charAt(i));
        }
        enzyme.setCvTerm(new CvTerm("PSI-MS", "MS:1001956", "NoEnzyme", null));
        enzymes.add(enzyme);


        enzyme = new Enzyme("Trypsin");
        enzyme.addAminoAcidBefore('R');
        enzyme.addAminoAcidBefore('K');
        enzyme.addRestrictionAfter('P');
        enzyme.setCvTerm(new CvTerm("PSI-MS", "MS:1001251", "Trypsin", null));
        enzymes.add(enzyme);

        enzyme = new Enzyme("Trypsin (no P rule)");
        enzyme.addAminoAcidBefore('R');
        enzyme.addAminoAcidBefore('K');
        enzyme.setCvTerm(new CvTerm("PSI-MS", "MS:1001313", "Trypsin/P", null));
        enzymes.add(enzyme);

        enzyme = new Enzyme("Arg-C");
        enzyme.addAminoAcidBefore('R');
        enzyme.addRestrictionAfter('P');
        enzyme.setCvTerm(new CvTerm("PSI-MS", "MS:1001303", "Arg-C", null));
        enzymes.add(enzyme);

        enzyme = new Enzyme("Arg-C (no P rule)");
        enzyme.addAminoAcidBefore('R');
        enzymes.add(enzyme);

        enzyme = new Enzyme("Arg-N");
        enzyme.addAminoAcidAfter('R');
        enzymes.add(enzyme);

        enzyme = new Enzyme("Glu-C");
        enzyme.addAminoAcidBefore('E');
        enzyme.setCvTerm(new CvTerm("PSI-MS", "MS:1001917", "glutamyl endopeptidase", null));
        enzymes.add(enzyme);

        enzyme = new Enzyme("Lys-C");
        enzyme.addAminoAcidBefore('K');
        enzyme.addRestrictionAfter('P');
        enzyme.setCvTerm(new CvTerm("PSI-MS", "MS:1001309", "Lys-C", null));
        enzymes.add(enzyme);

        if(ind < 0 || ind > enzymes.size()){
            System.err.println("Please provide a valid enzyme number:"+ind);
            System.exit(0);
        }
        Cloger.getInstance().logger.info("Use enzyme:"+enzymes.get(ind).getName());
        return(enzymes.get(ind));
    }


    public static DigestionParameters getDigestionPreferences(String enzymeName, int enzymeMissedCleavages) {
        DigestionParameters digestionParameters = new DigestionParameters();
        if(EnzymeFactory.getInstance().enzymeLoaded(enzymeName)){
            digestionParameters.setCleavageParameter(DigestionParameters.CleavageParameter.enzyme);
            Enzyme enzyme = EnzymeFactory.getInstance().getEnzyme(enzymeName);
            digestionParameters.addEnzyme(enzyme);
            digestionParameters.setnMissedCleavages(enzymeName, enzymeMissedCleavages);
        }else if(enzymeName.equalsIgnoreCase("NoEnzyme")){
            digestionParameters.setCleavageParameter(DigestionParameters.CleavageParameter.unSpecific);
            Enzyme enzyme = getEnzymeByIndex(0);
            digestionParameters.addEnzyme(enzyme);
            digestionParameters.setnMissedCleavages(enzymeName, 1000);
        }else{
            System.err.println("No valid enzyme:"+enzymeName);
            System.exit(0);
        }

        return digestionParameters;
    }

}
