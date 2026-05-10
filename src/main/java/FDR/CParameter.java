package main.java.FDR;

import java.io.IOException;
import java.io.InputStream;
import java.util.Properties;

public class CParameter {
    public static boolean clip_nTerm_M = true;
    public static String decoy_prefix = "rev_";
    public static int minPeptideLength = 7;
    public static int maxPeptideLength = 45;
    public static int maxMissedCleavages = 2;
    public static int enzyme = 1;
    public static int cpu = 0;

    public static String getVersion() {
        Properties properties = new Properties();
        try (InputStream in = CParameter.class.getClassLoader()
                .getResourceAsStream("project.properties")) {
            if (in != null) {
                properties.load(in);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return properties.getProperty("version");
    }
}
