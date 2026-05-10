package main.java.gui;

import java.util.Arrays;

import main.java.FDR.FDREval;

/**
 * Unified entry point for FDRBench:
 * <ul>
 *   <li>No arguments (or {@code --gui} / {@code -gui}) → launch the Swing GUI.
 *       This is what the {@code jpackage}-built {@code FDRBench.exe} hits when
 *       the user double-clicks it.</li>
 *   <li>Any other arguments → forward to the CLI ({@link FDREval#main}).</li>
 * </ul>
 *
 * The jar manifest's {@code Main-Class} (set in {@code pom.xml}) points here,
 * so {@code java -jar fdrbench-X.Y.Z.jar} and {@code FDRBench.exe} both end
 * up invoking this class — keeping the desktop launcher and the CLI on the
 * same single jar.
 */
public class FDRBenchLauncher {

    private static final String BANNER =
            "FDRBench - FDR Control Evaluation Tool for Proteomics\n" +
            "https://github.com/Noble-Lab/FDRBench\n";

    public static void main(String[] args) {
        if (args.length == 0 || containsGuiFlag(args)) {
            launchGUI(filterGuiFlags(args));
        } else {
            launchCLI(args);
        }
    }

    private static boolean containsGuiFlag(String[] args) {
        for (String a : args) {
            if ("--gui".equalsIgnoreCase(a) || "-gui".equalsIgnoreCase(a)) {
                return true;
            }
        }
        return false;
    }

    private static String[] filterGuiFlags(String[] args) {
        return Arrays.stream(args)
                .filter(a -> !"--gui".equalsIgnoreCase(a) && !"-gui".equalsIgnoreCase(a))
                .toArray(String[]::new);
    }

    private static void launchGUI(String[] args) {
        // GUI startup writes to stdout via FDRBenchGUI.main; that's harmless
        // when run from a console and ignored when run from the .exe.
        System.out.println(BANNER);
        System.out.println("Launching FDRBench GUI...\n");
        FDRBenchGUI.main(args);
    }

    private static void launchCLI(String[] args) {
        System.out.println(BANNER);
        try {
            FDREval.main(args);
        } catch (Exception e) {
            System.err.println("Error running FDRBench: " + e.getMessage());
            e.printStackTrace();
            System.exit(1);
        }
    }
}
