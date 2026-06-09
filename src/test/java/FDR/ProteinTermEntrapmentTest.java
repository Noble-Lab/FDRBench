package test.java.FDR;

import main.java.FDR.CParameter;
import main.java.FDR.DBGear;
import main.java.FDR.FDREval;

import org.junit.After;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertTrue;

/**
 * Tests for the protein N-/C-terminus fixing feature ({@code -fix_protein_nc}) in
 * {@link FDREval#generate_protein_entrapment_database}.
 *
 * <p>Each test drives the real generation path on a single protein and inspects the
 * entrapment protein in the output FASTA. The two sequences are the README example
 * heavy-chain protein with and without the leading initiator Met, so both an
 * {@code M}-start and a non-{@code M} start protein N-terminus are exercised.</p>
 */
public class ProteinTermEntrapmentTest {

    /** README example protein (heavy chain), starts with the initiator Met. */
    private static final String SEQ_M =
            "MDCTWRILLLVAAATGTHAEVQLVQSGAEVKKPGATVKISCKVSGYTFTDYYMHWVQQAPGKGLEWMGLVDPEDGETIYAEKFQGRVTITADTSTDTAYMELSSLRSEDTAVYYCAT";

    /** Same protein with the leading Met removed: N-terminus is now {@code D}. */
    private static final String SEQ_NO_M =
            "DCTWRILLLVAAATGTHAEVQLVQSGAEVKKPGATVKISCKVSGYTFTDYYMHWVQQAPGKGLEWMGLVDPEDGETIYAEKFQGRVTITADTSTDTAYMELSSLRSEDTAVYYCAT";

    @BeforeClass
    public static void initEnzymes() {
        DBGear.init_enzymes();
    }

    /** Restore the mutated static config to its documented defaults so tests can't leak state. */
    @After
    public void resetStaticConfig() {
        boolean[] d = FDREval.resolveFixProteinNc(null);
        FDREval.fix_protein_n_term = d[0];
        FDREval.fix_protein_c_term = d[1];
        FDREval.fix_nc_aa = FDREval.Fix_NC.NC;
        FDREval.add_decoy = true;
        FDREval.I2L = false;
        FDREval.for_diann = false;
        FDREval.for_uniprot = false;
        FDREval.global_random_seed = 2000;
        CParameter.enzyme = 1;
        CParameter.minPeptideLength = 7;
        CParameter.maxPeptideLength = 45;
    }

    @Test
    public void proteinNTermPreservedForMetStart() throws IOException {
        assertProteinNTermPreserved(SEQ_M);
    }

    @Test
    public void proteinNTermPreservedForNonMetStart() throws IOException {
        assertProteinNTermPreserved(SEQ_NO_M);
    }

    /**
     * With {@code -fix_nc c -fix_protein_nc n} the protein N-terminal residue is pinned,
     * the per-peptide C-term stays fixed (so the protein C-term is preserved too), and the
     * entrapment is still an equal-length permutation that differs from the {@code n} off run.
     */
    private void assertProteinNTermPreserved(String target) throws IOException {
        String off = generateEntrapment(target, FDREval.Fix_NC.C, false, false); // -fix_nc c
        String on  = generateEntrapment(target, FDREval.Fix_NC.C, true,  false); // -fix_nc c -fix_protein_nc n

        assertEquals("length preserved (off)", target.length(), off.length());
        assertEquals("length preserved (on)",  target.length(), on.length());
        assertTrue("entrapment is a permutation of target (off)", isPermutation(target, off));
        assertTrue("entrapment is a permutation of target (on)",  isPermutation(target, on));
        assertNotEquals("-fix_protein_nc n changes the entrapment", off, on);
        assertEquals("protein N-term residue preserved with -fix_protein_nc n",
                str(target.charAt(0)), str(on.charAt(0)));
        assertEquals("protein C-term residue preserved by -fix_nc c",
                str(lastChar(target)), str(lastChar(on)));
    }

    /**
     * With {@code -fix_nc n} the peptide C-terminus is free, so only {@code -fix_protein_nc c}
     * pins the protein C-terminal residue.
     */
    @Test
    public void proteinCTermPreservedWhenFixNcN() throws IOException {
        String target = SEQ_M; // ends in ...YYCAT -> C-term 'T'
        String off = generateEntrapment(target, FDREval.Fix_NC.N, false, false); // -fix_nc n (C free)
        String on  = generateEntrapment(target, FDREval.Fix_NC.N, false, true);  // -fix_nc n -fix_protein_nc c

        assertEquals("length preserved", target.length(), on.length());
        assertTrue("entrapment is a permutation of target", isPermutation(target, on));
        assertNotEquals("-fix_protein_nc c changes the entrapment", off, on);
        assertEquals("protein C-term residue preserved with -fix_protein_nc c",
                str(lastChar(target)), str(lastChar(on)));
    }

    /** {@code -fix_protein_nc nc} pins both protein termini. */
    @Test
    public void bothProteinTerminiPreservedWithNc() throws IOException {
        String target = SEQ_M;
        String on = generateEntrapment(target, FDREval.Fix_NC.N, true, true);

        assertEquals("protein N-term preserved", str(target.charAt(0)), str(on.charAt(0)));
        assertEquals("protein C-term preserved", str(lastChar(target)), str(lastChar(on)));
        assertTrue("entrapment is a permutation of target", isPermutation(target, on));
    }

    /** N-term fixing explicitly disabled: still an equal-length permutation (no regression). */
    @Test
    public void nTermFixingDisabledIsEqualLengthPermutation() throws IOException {
        String target = SEQ_NO_M;
        String off = generateEntrapment(target, FDREval.Fix_NC.C, false, false);

        assertEquals("length preserved", target.length(), off.length());
        assertTrue("entrapment is a permutation of target", isPermutation(target, off));
    }

    /**
     * The {@code -fix_protein_nc} parser maps each keyword (and the {@code null} default) to the
     * right {@code {nTermFixed, cTermFixed}} pair. Because the field defaults derive from
     * {@code resolveFixProteinNc(null)}, this also locks down the N-on-by-default decision.
     */
    @Test
    public void fixProteinNcParsingMapsKeywordsToFlags() {
        assertArrayEquals("default (null) = N only", new boolean[]{true,  false}, FDREval.resolveFixProteinNc(null));
        assertArrayEquals("n = N only",              new boolean[]{true,  false}, FDREval.resolveFixProteinNc("n"));
        assertArrayEquals("c = C only",              new boolean[]{false, true},  FDREval.resolveFixProteinNc("c"));
        assertArrayEquals("nc = both",               new boolean[]{true,  true},  FDREval.resolveFixProteinNc("nc"));
        assertArrayEquals("cn = both",               new boolean[]{true,  true},  FDREval.resolveFixProteinNc("cn"));
        assertArrayEquals("off = neither",           new boolean[]{false, false}, FDREval.resolveFixProteinNc("off"));
    }

    @Test(expected = IllegalArgumentException.class)
    public void fixProteinNcParserRejectsInvalidValue() {
        FDREval.resolveFixProteinNc("xyz");
    }

    /**
     * Builds a 1-fold, decoy-free protein entrapment database for a single protein and
     * returns the entrapment protein sequence.
     */
    private static String generateEntrapment(String proteinSeq, FDREval.Fix_NC fixNc,
                                             boolean fixProteinN, boolean fixProteinC) throws IOException {
        // All public static knobs set explicitly; everything else stays at its default.
        FDREval.fix_nc_aa = fixNc;
        FDREval.fix_protein_n_term = fixProteinN;
        FDREval.fix_protein_c_term = fixProteinC;
        FDREval.add_decoy = false;        // keep output to target + entrapment only
        FDREval.I2L = false;
        FDREval.for_diann = false;
        FDREval.for_uniprot = false;
        FDREval.global_random_seed = 2000;
        CParameter.enzyme = 2;            // Trypsin (no P rule), the -db default
        CParameter.minPeptideLength = 7;
        CParameter.maxPeptideLength = 35;

        File db = File.createTempFile("fdrbench_db", ".fasta");
        db.deleteOnExit();
        try (BufferedWriter w = new BufferedWriter(new FileWriter(db))) {
            w.write(">TEST\n" + proteinSeq + "\n");
        }
        File out = File.createTempFile("fdrbench_out", ".fasta");
        out.deleteOnExit();

        FDREval.generate_protein_entrapment_database(db.getAbsolutePath(), out.getAbsolutePath(), 1);

        String entrapment = null;
        try (BufferedReader r = new BufferedReader(new FileReader(out))) {
            String header = null, line;
            while ((line = r.readLine()) != null) {
                if (line.startsWith(">")) {
                    header = line.substring(1);
                } else if (header != null) {
                    if (header.contains(FDREval.entrapment_label)) {
                        entrapment = line;
                    }
                    header = null;
                }
            }
        }
        if (entrapment == null) {
            throw new IllegalStateException("No entrapment protein found in output for: " + proteinSeq);
        }
        return entrapment;
    }

    private static char lastChar(String s) {
        return s.charAt(s.length() - 1);
    }

    private static String str(char c) {
        return String.valueOf(c);
    }

    /** True when {@code a} and {@code b} contain exactly the same multiset of characters. */
    private static boolean isPermutation(String a, String b) {
        if (a.length() != b.length()) {
            return false;
        }
        char[] ca = a.toCharArray();
        char[] cb = b.toCharArray();
        Arrays.sort(ca);
        Arrays.sort(cb);
        return Arrays.equals(ca, cb);
    }
}
