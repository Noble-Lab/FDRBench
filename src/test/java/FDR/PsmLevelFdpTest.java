package test.java.FDR;

import main.java.FDR.FDREval;
import main.java.FDR.FDRType;

import org.junit.Test;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * Tests for PSM-level FDP estimation via the combined entrapment method
 * ({@link FDREval#calc_fdp_fast_combined_entrapment_method}), the fix for issue #15.
 *
 * <p>PSM-level input has one row per PSM, so the same peptide can appear in several rows.
 * The combined entrapment method counts every row, so the cumulative target/entrapment
 * counts must reflect each PSM independently (no collapsing to unique peptides).</p>
 */
public class PsmLevelFdpTest {

    /**
     * Five PSMs over three peptides (two peptides repeat) drive the combined entrapment
     * method at the PSM level. With r=1 the combined FDP is 2 * n_p / (n_p + n_t), and the
     * cumulative counts at every threshold must add up to the number of PSMs seen so far,
     * proving rows are counted per-PSM rather than per-unique-peptide.
     */
    @Test
    public void psmLevelCountsEveryPsm() throws IOException {
        // peptide_type: TARGETA and TARGETC are targets, ENTRAPB is an entrapment hit.
        File pep = tempFile("fdrbench_pep", ".txt");
        try (BufferedWriter w = new BufferedWriter(new FileWriter(pep))) {
            w.write("sequence\tdecoy\tproteins\tpeptide_type\tpeptide_pair_index\n");
            w.write("TARGETA\tNo\tprotA\ttarget\t0\n");
            w.write("ENTRAPB\tNo\tprotB_p_target\tp_target\t0\n");
            w.write("TARGETC\tNo\tprotC\ttarget\t1\n");
        }

        // Five PSMs, deliberately out of q_value order; TARGETA and ENTRAPB each appear twice.
        File psm = tempFile("fdrbench_psm", ".tsv");
        try (BufferedWriter w = new BufferedWriter(new FileWriter(psm))) {
            w.write("peptide\tq_value\tscore\n");
            w.write("ENTRAPB\t0.005\t6\n"); // sorted position 4
            w.write("TARGETA\t0.001\t10\n"); // sorted position 0
            w.write("TARGETC\t0.004\t7\n"); // sorted position 3
            w.write("TARGETA\t0.002\t9\n"); // sorted position 1
            w.write("ENTRAPB\t0.003\t8\n"); // sorted position 2
        }

        File out = tempFile("fdrbench_out", ".csv");

        FDREval fdrEval = new FDREval();
        fdrEval.fdp_level = FDRType.psm;
        fdrEval.score_column_name = "-"; // rank by q_value only
        fdrEval.calc_fdp_fast_combined_entrapment_method(
                psm.getAbsolutePath(), pep.getAbsolutePath(), out.getAbsolutePath(), 1.0);

        List<HashMap<String, String>> rows = readCsv(out);
        assertEquals("one output row per PSM", 5, rows.size());

        // Rows are written in q_value ascending order:
        // 0: TARGETA  -> n_t=1 n_p=0
        // 1: TARGETA  -> n_t=2 n_p=0   (duplicate peptide still counted)
        // 2: ENTRAPB  -> n_t=2 n_p=1
        // 3: TARGETC  -> n_t=3 n_p=1
        // 4: ENTRAPB  -> n_t=3 n_p=2
        assertEquals("TARGETA", rows.get(0).get("peptide"));
        assertEquals("ENTRAPB", rows.get(2).get("peptide"));
        assertEquals("ENTRAPB", rows.get(4).get("peptide"));

        // Per-PSM counting: at the last threshold every PSM is counted (3 + 2 = 5),
        // whereas peptide-level deduplication would yield only 3 unique peptides.
        int nt = Integer.parseInt(rows.get(4).get("n_t"));
        int np = Integer.parseInt(rows.get(4).get("n_p"));
        assertEquals("targets counted per PSM", 3, nt);
        assertEquals("entrapments counted per PSM", 2, np);
        assertEquals("every PSM is counted", 5, nt + np);

        // combined_fdp = 2 * n_p / (n_p + n_t) at r = 1.
        assertEquals(2.0 * 1 / 3, parse(rows.get(2).get("combined_fdp")), 1e-9);
        assertEquals(2.0 * 2 / 5, parse(rows.get(4).get("combined_fdp")), 1e-9);

        // Duplicate target PSM before any entrapment is seen -> FDP still 0.
        assertEquals(0.0, parse(rows.get(1).get("combined_fdp")), 1e-9);
    }

    private static File tempFile(String prefix, String suffix) throws IOException {
        File f = File.createTempFile(prefix, suffix);
        f.deleteOnExit();
        return f;
    }

    private static double parse(String s) {
        return Double.parseDouble(s);
    }

    /** Reads a comma-separated file with a header row into a list of column-name -> value maps. */
    private static List<HashMap<String, String>> readCsv(File file) throws IOException {
        List<HashMap<String, String>> rows = new ArrayList<>();
        try (BufferedReader r = new BufferedReader(new FileReader(file))) {
            String header = r.readLine();
            assertTrue("output has a header", header != null);
            String[] cols = header.split(",");
            String line;
            while ((line = r.readLine()) != null) {
                String[] vals = line.split(",");
                HashMap<String, String> row = new HashMap<>();
                for (int i = 0; i < cols.length && i < vals.length; i++) {
                    row.put(cols[i], vals[i]);
                }
                rows.add(row);
            }
        }
        return rows;
    }
}
