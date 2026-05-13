package main.java.gui;

import com.formdev.flatlaf.FlatDarkLaf;
import com.formdev.flatlaf.FlatLaf;
import com.formdev.flatlaf.FlatLightLaf;
import main.java.FDR.CParameter;
import main.java.FDR.DBGear;
import javax.swing.*;
import javax.swing.filechooser.FileNameExtensionFilter;
import java.awt.*;
import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.Locale;
import java.util.prefs.Preferences;

/**
 * FDRBench GUI - FDR Control Evaluation Tool
 * Two workflows: Entrapment Database Generation and FDP Estimation
 */
public class FDRBenchGUI extends JFrame {

    // Workflow selection
    private static final int WORKFLOW_DB_GENERATION = 0;
    private static final int WORKFLOW_FDP_ESTIMATION = 1;
    private int currentWorkflow = WORKFLOW_DB_GENERATION;

    // Brand colors
    private static final Color ACCENT_COLOR = new Color(46, 204, 113);

    // Layout constants
    private static final int ROW_SPACING = 4;
    private static final int COL_SPACING = 8;
    private static final Insets DEFAULT_INSETS = new Insets(ROW_SPACING, COL_SPACING, ROW_SPACING, COL_SPACING);
    private static final int DEFAULT_WIDTH = 700;
    private static final int DEFAULT_HEIGHT = 750;
    private static final int MIN_WIDTH = 700;
    private static final int MIN_HEIGHT = 750;
    // Shared height for all "input boxes" (text fields, spinners, combos) so
    // every row in the parameters panel lines up cleanly.
    // Package-private so FdpPlotPanel can use the same
    // value for its FDP-file input row.
    static final int COMPONENT_HEIGHT = 28;

    // Workflow dropdown
    private JComboBox<String> workflowCombo;

    // Input fields - Workflow 1: DB Generation
    private JTextField dbFileField;
    private JTextField foreignSpeciesField;
    // Backing list for the Foreign Species multi-file picker, when more than one file is selected, the
    // text field renders a "(N files selected)" hyperlink and the actual
    // paths live here. Used by buildCommandArgs to assemble the -ms argument.
    private final java.util.List<String> foreignSpeciesFiles = new java.util.ArrayList<>();
    private JTextField dbOutputField;
    private JComboBox<String> dbLevelCombo;
    private JComboBox<String> enzymeCombo;
    private JSpinner missCleavageSpinner;
    private JSpinner minLengthSpinner;
    private JSpinner maxLengthSpinner;
    private JSpinner foldSpinner;
    private JCheckBox clipNmCheckbox;
    private JComboBox<String> fixNcCombo;
    private JCheckBox i2lCheckbox;
    private JCheckBox decoyCheckbox;
    private JCheckBox diannCheckbox;
    private JCheckBox uniprotCheckbox;
    private JCheckBox exportDbCheckbox;
    private JCheckBox checkDuplicatesCheckbox;
    private JCheckBox noSharedCheckbox;

    // Top-level selectors (hoisted above the workflow cards). The Level combo
    // is workflow-dependent so we keep two combos side-by-side in a CardLayout
    // panel. Sequence Generation is shared between both workflows.
    private JComboBox<String> seqMethodCombo;
    private CardLayout levelCardLayout;
    private JPanel levelCardPanel;
    private List<JComponent> foreignSpeciesRowComponents;
    // Tracks the R Ratio row so it can be hidden in the FDP workflow when
    // Sequence Generation is set to Random Shuffling (R only matters when
    // entrapments come from foreign species).
    private List<JComponent> rRatioRowComponents;
    // Tracks the FDP Fold row so it can be hidden when method = Foreign
    // Species — k-fold FDP only applies to random-shuffling entrapment;
    // foreign-species runs use the R Ratio instead.
    private List<JComponent> fdpFoldRowComponents;
    // Tracks the FDP Peptide Pair File row so it can be hidden when the
    // calculation level is "protein" (pair tracking is peptide-level only).
    private List<JComponent> pepPairRowComponents;
    // Tracks the Decoy Label and Decoy Label Pos rows so they can be hidden
    // when the Add Decoys checkbox is off (decoy values aren't sent to the
    // CLI in that case).
    private List<JComponent> decoyLabelRowComponents;
    // Tracks the Export Protein DB row so it can be hidden at protein level
    // (where the main output is already a protein FASTA, making the export
    // redundant).
    private List<JComponent> exportDbRowComponents;
    // Tracks the Check Duplicates row so it can be hidden unless we're in
    // protein-level random-shuffling mode (the only path FDREval actually
    // honors -check on; see generate_protein_entrapment_database).
    private List<JComponent> checkDupRowComponents;
    // Tracks the No Shared Peptides row so it can be hidden unless we're in
    // protein-level foreign-species mode (the only path FDREval actually
    // honors -ns on; see generate_protein_entrapment_database_from_multiple_species_data).
    private List<JComponent> noSharedRowComponents;
    // Tracks the Add Decoys row so the whole option can be hidden when
    // method = Foreign Species — generated foreign-species entrapments
    // don't get decoyed, so the flag has no effect there.
    private List<JComponent> addDecoysRowComponents;
    // Rows hidden whenever Level = protein (Missed Cleavages, Peptide Length,
    // Clip N-terminal M) — these have no meaning when the output is whole
    // proteins.
    private List<JComponent> peptideOnlyRowComponents;
    // Rows needed for digestion (Enzyme, Fix N/C Terminal). Visible at peptide
    // level always, and at protein level only when the user picked Random
    // Shuffling — protein-level shuffling still digests the target before
    // shuffling and reassembling.
    private List<JComponent> digestionRowComponents;
    // Fix N/C Terminal — separate because it's only used by the entrapment
    // shuffler, NOT by the dedup-only digest path that runs at protein +
    // foreign species + No Shared Peptides.
    private List<JComponent> fixNcRowComponents;
    // Reference kept so we can hide the whole Digestion Settings section when
    // every row inside it would be hidden (protein + Foreign Species).
    private JPanel dbDigestionPanel;
    private static final int SEQ_METHOD_SHUFFLING = 0;
    private static final int SEQ_METHOD_FOREIGN = 1;

    // Input fields - Workflow 2: FDP Estimation
    private JTextField inputFileField;
    private JTextField pepPairFileField;
    private JTextField fdpOutputField;
    private JComboBox<String> fdpLevelCombo;
    private JComboBox<String> scoreColumnCombo;
    // Sentinel item meaning "do not pass -score" so FDRBench ranks by q_value
    // alone. Must not collide with any real column name a user might choose.
    private static final String NO_SCORE_ITEM = "none";
    private JComboBox<String> scoreDirectionCombo;
    private JSpinner fdpFoldSpinner;
    private JComboBox<String> pickMethodCombo;
    private JTextField rRatioField;
    private JTextField fdpEntrapmentLabelField;
    private JComboBox<String> fdpEntrapmentPosCombo;

    // Advanced settings (shared)
    private JSpinner seedSpinner;
    private JTextField entrapmentLabelField;
    private JTextField decoyLabelField;
    private JComboBox<String> entrapmentPosCombo;
    private JComboBox<String> decoyPosCombo;
    private JCheckBox debugCheckbox;

    // Console and UI components
    private JTextArea consoleArea;
    private JProgressBar progressBar;
    private JButton runButton;
    private JButton stopButton;
    private JTabbedPane tabbedPane;
    private JLabel statusLabel;

    // Plot tab — only attached to the tabbed pane while the FDP Estimation
    // workflow is selected.
    private FdpPlotPanel fdpPlotPanel;
    private static final String PLOT_TAB_TITLE = "Plot";

    // Header components
    private JPanel headerPanel;
    private JLabel headerTitleLabel;
    private JLabel headerSubtitleLabel;
    private JLabel headerVersionLabel;
    private JToggleButton darkModeToggle;
    private JToggleButton particleToggle;

    // Panel references needed across methods
    private JPanel parametersContentPanel;
    private JPanel workflowCardsPanel;
    private CardLayout workflowCardLayout;

    // Execution
    private Process currentProcess;
    private volatile boolean isRunning = false;
    // Mirrors the on-screen Console output to fdrbench_log.txt inside the
    // chosen output folder for the duration of a run.
    private BufferedWriter logWriter;

    // Preferences
    private static final Preferences prefs = Preferences.userNodeForPackage(FDRBenchGUI.class);
    private static final String PREF_LAST_DIR = "lastDirectory";
    private static final String PREF_DARK_MODE = "darkMode";

    public FDRBenchGUI() {
        setTitle("FDRBench - FDR Control Evaluation Tool");
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        setMinimumSize(new Dimension(MIN_WIDTH, MIN_HEIGHT));
        setPreferredSize(new Dimension(DEFAULT_WIDTH, DEFAULT_HEIGHT));

        // Hide title bar icon (FlatLaf-painted title bar) — the same image is
        // still used for the taskbar / dock and inside the header below.
        getRootPane().putClientProperty("JRootPane.titleBarShowIcon", false);

        // Load Application Icon for the OS taskbar / FlatLaf title bar.
        // (java.awt.Taskbar — for the macOS dock — is Java 9+ only and this
        // project targets 1.8, so we stick with setIconImage which already
        // handles Windows / Linux taskbars correctly.)
        try {
            java.net.URL iconUrl = getClass().getResource("/fdrbench-icon.png");
            if (iconUrl != null) {
                setIconImage(new ImageIcon(iconUrl).getImage());
            }
        } catch (Exception e) {
            System.err.println("Failed to load application icon: " + e.getMessage());
        }

        // Load persisted theme preference
        boolean dark = prefs.getBoolean(PREF_DARK_MODE, false);

        try {
            if (dark) {
                FlatDarkLaf.setup();
            } else {
                FlatLightLaf.setup();
            }
            customizeUIDefaults();
        } catch (Exception e) {
            System.err.println("Theme setup failed: " + e.getMessage());
        }

        initComponents();
        updateHeaderForegrounds();

        pack();

        // sizing strategy: on monitors taller than 1080px give
        // the window a little extra room, capped at the screen size; on smaller
        // monitors fall back to the larger of the packed size and the minimum.
        Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
        int newWidth;
        int newHeight;
        if (screenSize.height > 1080) {
            int targetWidth = DEFAULT_WIDTH + 100;
            int targetHeight = DEFAULT_HEIGHT + 100;
            newWidth = Math.max(MIN_WIDTH, Math.min(targetWidth, screenSize.width));
            newHeight = Math.max(MIN_HEIGHT, Math.min(targetHeight, screenSize.height));
        } else {
            Dimension packedSize = getSize();
            Dimension minSize = getMinimumSize();
            newWidth = Math.max(packedSize.width, minSize.width);
            newHeight = Math.max(packedSize.height, minSize.height);
        }
        setSize(newWidth, newHeight);
        setLocationRelativeTo(null);
    }

    /**
     * Returns a Font derived from the LAF's default with the given style and
     * size. Avoids the Windows-only "Segoe UI" hardcoding so the GUI renders
     * cleanly on Linux/macOS too.
     */
    private static Font derivedFont(int style, float size) {
        Font base = UIManager.getFont("defaultFont");
        if (base == null) {
            base = UIManager.getFont("Label.font");
        }
        if (base == null) {
            base = new Font(Font.SANS_SERIF, Font.PLAIN, 13);
        }
        return base.deriveFont(style, size);
    }

    private static void customizeUIDefaults() {
        // Derive from the LAF default so we don't depend on a Windows-only font.
        Font defaultFont = UIManager.getFont("Label.font");
        if (defaultFont != null) {
            UIManager.put("defaultFont", defaultFont.deriveFont(13f));
        } else {
            UIManager.put("defaultFont", new Font(Font.SANS_SERIF, Font.PLAIN, 13));
        }
        UIManager.put("Button.arc", 10);
        UIManager.put("Component.arc", 10);
        UIManager.put("ProgressBar.arc", 10);
        UIManager.put("TextComponent.arc", 8);
        UIManager.put("Component.focusWidth", 1);
        UIManager.put("Component.innerFocusWidth", 0);
        UIManager.put("Component.hideMnemonics", true);
        UIManager.put("TabbedPane.showTabSeparators", true);
        UIManager.put("TabbedPane.tabInsets", new Insets(8, 14, 8, 14));
        UIManager.put("ScrollBar.width", 12);
        UIManager.put("ScrollBar.thumbArc", 999);
        UIManager.put("ScrollBar.thumbInsets", new Insets(2, 2, 2, 2));
        UIManager.put("CheckBox.arc", 4);
        UIManager.put("Spinner.arc", 8);
        ToolTipManager.sharedInstance().setDismissDelay(30000);
    }

    private void initComponents() {
        setLayout(new BorderLayout());

        // Header
        add(createHeader(), BorderLayout.NORTH);

        // Tabbed pane
        tabbedPane = new JTabbedPane();
        tabbedPane.setTabLayoutPolicy(JTabbedPane.SCROLL_TAB_LAYOUT);
        // Font inherits from defaultFont (set in customizeUIDefaults).
        tabbedPane.addTab("Workflow", createParametersPanel());
        // The Plot tab is owned by FdpPlotPanel; it gets attached / detached
        // dynamically when the user toggles between workflows.
        fdpPlotPanel = new FdpPlotPanel();
        tabbedPane.addTab("Console", createConsolePanel());

        JPanel mainPanel = new JPanel(new BorderLayout(10, 10));
        mainPanel.setBorder(BorderFactory.createEmptyBorder(10, 15, 10, 15));
        mainPanel.add(tabbedPane, BorderLayout.CENTER);

        add(mainPanel, BorderLayout.CENTER);

        // Footer
        add(createFooter(), BorderLayout.SOUTH);

        // Enter key triggers Run
        getRootPane().setDefaultButton(runButton);

        // Initial visibility
        updateWorkflowPanelVisibility();
    }

    // ==================== HEADER ====================

    // Constants pulled out to the outer class because Java 8 forbids static
    // declarations in non-static inner classes (DynamicHeaderPanel needs to be
    // non-static so its updateUI() override can reach updateHeaderForegrounds).
    private static final int PARTICLE_COUNT = 80;
    private static final double CONNECTION_THRESHOLD = 130.0;

    private class DynamicHeaderPanel extends JPanel {
        private final List<Particle> particles = new ArrayList<>();
        private final javax.swing.Timer timer;
        private final Color[] WHITE_ALPHA;
        private boolean animationEnabled = true;

        @Override
        public void updateUI() {
            super.updateUI();
            // After a LAF switch, re-derive header text colors from the new
            // theme and force a repaint so the gradient picks up the new
            // FlatLaf.isLafDark() value.
            updateHeaderForegrounds();
            repaint();
        }

        DynamicHeaderPanel() {
            super(new BorderLayout());

            // Pre-build alpha-colour table for fast particle rendering
            WHITE_ALPHA = new Color[256];
            for (int i = 0; i < 256; i++) {
                WHITE_ALPHA[i] = new Color(255, 255, 255, i);
            }

            for (int i = 0; i < PARTICLE_COUNT; i++) {
                particles.add(new Particle());
            }

            timer = new javax.swing.Timer(33, e -> {
                if (!animationEnabled)
                    return;
                int w = getWidth() > 0 ? getWidth() : 800;
                int h = getHeight() > 0 ? getHeight() : 150;
                for (Particle p : particles)
                    p.update(w, h);
                repaint();
            });
            timer.start();

            // Stop the animation loop when the panel isn't actually drawing
            // (window minimized, tab hidden) so we don't burn CPU in the
            // background.
            addHierarchyListener(e -> {
                if ((e.getChangeFlags() & java.awt.event.HierarchyEvent.SHOWING_CHANGED) != 0) {
                    if (isShowing() && animationEnabled) {
                        timer.start();
                    } else {
                        timer.stop();
                    }
                }
            });
        }

        void setAnimationEnabled(boolean enabled) {
            this.animationEnabled = enabled;
            if (enabled)
                timer.start();
            else
                timer.stop();
            repaint();
        }

        @Override
        protected void paintComponent(Graphics g) {
            super.paintComponent(g);
            int w = getWidth(), h = getHeight();
            if (w <= 0 || h <= 0)
                return;

            Graphics2D g2 = (Graphics2D) g.create();
            g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

            boolean dark = FlatLaf.isLafDark();
            Color base = new Color(0x2F82B7);
            Color top    = dark ? adjust(base, -40) : adjust(base, 60);
            Color mid    = dark ? adjust(base, -20) : adjust(base, 35);
            Color bottom = dark ? adjust(base, -10) : adjust(base, 15);

            g2.setPaint(new LinearGradientPaint(0f, 0f, 0f, h,
                    new float[] { 0f, 0.5f, 1f }, new Color[] { top, mid, bottom }));
            g2.fillRect(0, 0, w, h);

            if (animationEnabled) {
                drawParticles(g2);
            }

            g2.dispose();
        }

        private void drawParticles(Graphics2D g2) {
            double thresholdSq = CONNECTION_THRESHOLD * CONNECTION_THRESHOLD;
            g2.setStroke(new BasicStroke(1.0f));

            for (int i = 0; i < particles.size(); i++) {
                Particle p1 = particles.get(i);
                for (int j = i + 1; j < particles.size(); j++) {
                    Particle p2 = particles.get(j);
                    double distSq = p1.distanceSq(p2);
                    if (distSq < thresholdSq) {
                        double dist = Math.sqrt(distSq);
                        int alpha = (int) ((1.0 - (dist / CONNECTION_THRESHOLD)) * 80);
                        g2.setColor(WHITE_ALPHA[Math.min(255, Math.max(0, alpha))]);
                        g2.drawLine((int) p1.x, (int) p1.y, (int) p2.x, (int) p2.y);
                    }
                }
            }

            for (Particle p : particles) {
                int alphaInt = (int) (p.alpha * 255);
                g2.setColor(WHITE_ALPHA[Math.min(255, Math.max(0, alphaInt))]);
                int size = (int) p.size;
                g2.fillOval((int) (p.x - size / 2.0), (int) (p.y - size / 2.0), size, size);
            }
        }

        class Particle {
            double x, y, vx, vy, size, alpha;

            Particle() {
                reset(800, 150);
            }

            void reset(int w, int h) {
                x = Math.random() * (w > 0 ? w : 800);
                y = Math.random() * (h > 0 ? h : 150);
                double speed = 0.35 + Math.random() * 0.55;
                double angle = Math.random() * 2 * Math.PI;
                vx = Math.cos(angle) * speed;
                vy = Math.sin(angle) * speed;
                size = 2.0 + Math.random() * 2.5;
                alpha = 0.2 + Math.random() * 0.4;
            }

            double distanceSq(Particle o) {
                double dx = x - o.x, dy = y - o.y;
                return dx * dx + dy * dy;
            }

            void update(int w, int h) {
                x += vx;
                y += vy;
                
                // Keep particles within slightly larger bounds to avoid edges
                double margin = 20;
                if (x < -margin) x = w + margin;
                else if (x > w + margin) x = -margin;
                if (y < -margin) y = h + margin;
                else if (y > h + margin) y = -margin;
            }
        }
    }

    private JPanel createHeader() {
        DynamicHeaderPanel dhPanel = new DynamicHeaderPanel();
        headerPanel = dhPanel;
        headerPanel.setOpaque(false);
        headerPanel.setBorder(BorderFactory.createEmptyBorder(20, 25, 20, 25));

        JPanel titlePanel = new JPanel(new FlowLayout(FlowLayout.LEFT, 15, 0));
        titlePanel.setOpaque(false);

        // Header icon — loaded from /fdrbench-icon.png and rendered at 48×48
        // with a rounded-rectangle clip. The
        // 128×128 scaled image keeps it sharp on HiDPI displays without
        // blowing up paint cost.
        JLabel headerIconLabel = new JLabel("");
        try {
            java.net.URL iconUrl = getClass().getResource("/fdrbench-icon.png");
            if (iconUrl != null) {
                ImageIcon originalIcon = new ImageIcon(iconUrl);
                final Image highResImage = originalIcon.getImage()
                        .getScaledInstance(128, 128, Image.SCALE_SMOOTH);
                headerIconLabel.setIcon(new javax.swing.Icon() {
                    @Override
                    public void paintIcon(Component c, Graphics g, int x, int y) {
                        Graphics2D g2 = (Graphics2D) g.create();
                        try {
                            g2.setRenderingHint(RenderingHints.KEY_INTERPOLATION,
                                    RenderingHints.VALUE_INTERPOLATION_BICUBIC);
                            g2.setRenderingHint(RenderingHints.KEY_RENDERING,
                                    RenderingHints.VALUE_RENDER_QUALITY);
                            g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
                                    RenderingHints.VALUE_ANTIALIAS_ON);
                            g2.setClip(new java.awt.geom.RoundRectangle2D.Float(
                                    x, y, 48, 48, 16, 16));
                            g2.drawImage(highResImage, x, y, 48, 48, null);
                        } finally {
                            g2.dispose();
                        }
                    }
                    @Override public int getIconWidth()  { return 48; }
                    @Override public int getIconHeight() { return 48; }
                });
            }
        } catch (Exception e) {
            // Header icon is non-essential — leave the label blank if loading fails.
        }

        JPanel textPanel = new JPanel();
        textPanel.setOpaque(false);
        textPanel.setLayout(new BoxLayout(textPanel, BoxLayout.Y_AXIS));

        headerTitleLabel = new JLabel("FDRBench");
        headerTitleLabel.setFont(derivedFont(Font.BOLD, 28f));

        headerSubtitleLabel = new JLabel("FDR Control Evaluation Tool for Proteomics");
        headerSubtitleLabel.setFont(derivedFont(Font.PLAIN, 13f));

        textPanel.add(headerTitleLabel);
        textPanel.add(Box.createVerticalStrut(3));
        textPanel.add(headerSubtitleLabel);
        titlePanel.add(headerIconLabel);
        titlePanel.add(textPanel);

        JPanel rightPanel = new JPanel(new FlowLayout(FlowLayout.RIGHT, 15, 0));
        rightPanel.setOpaque(false);

        particleToggle = createHeaderToggleButton("Effects On", true, e -> {
            boolean sel = particleToggle.isSelected();
            particleToggle.setText(sel ? "Effects On" : "Effects Off");
            dhPanel.setAnimationEnabled(sel);
        });

        boolean dark = FlatLaf.isLafDark();
        darkModeToggle = createHeaderToggleButton(dark ? "Light Mode" : "Dark Mode", dark,
                e -> toggleDarkMode(darkModeToggle.isSelected()));

        rightPanel.add(particleToggle);
        rightPanel.add(darkModeToggle);

        String version = CParameter.getVersion();
        headerVersionLabel = new JLabel(version != null ? "v" + version : "");
        headerVersionLabel.setFont(derivedFont(Font.PLAIN, 12f));
        rightPanel.add(headerVersionLabel);

        headerPanel.add(titlePanel, BorderLayout.WEST);
        headerPanel.add(rightPanel, BorderLayout.EAST);

        // Keep the toggle buttons clickable even when the subtitle text would
        // otherwise overlap them at narrow widths.
        headerPanel.setComponentZOrder(rightPanel, 0);

        return headerPanel;
    }

    private JToggleButton createHeaderToggleButton(String text, boolean selected,
            java.awt.event.ActionListener action) {
        JToggleButton btn = new JToggleButton(text) {
            @Override
            protected void paintComponent(Graphics g) {
                Graphics2D g2 = (Graphics2D) g.create();
                g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
                boolean dark = FlatLaf.isLafDark();
                Color bg = dark ? new Color(0, 0, 0, 60) : new Color(255, 255, 255, 60);
                if (getModel().isRollover())
                    bg = withAlpha(bg, 100);
                g2.setColor(bg);
                g2.fillRoundRect(0, 0, getWidth(), getHeight(), 15, 15);
                g2.setColor(withAlpha(getForeground(), 80));
                g2.drawRoundRect(0, 0, getWidth() - 1, getHeight() - 1, 15, 15);
                g2.dispose();
                super.paintComponent(g);
            }
        };
        btn.setSelected(selected);
        btn.setContentAreaFilled(false);
        btn.setBorderPainted(false);
        btn.setOpaque(false);
        btn.setFocusPainted(false);
        btn.setMargin(new Insets(6, 16, 6, 16));
        btn.setFont(btn.getFont().deriveFont(12f));
        btn.setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
        btn.setPreferredSize(new Dimension(100, 30));
        btn.addActionListener(action);
        return btn;
    }

    private void updateHeaderForegrounds() {
        boolean dark = FlatLaf.isLafDark();
        Color base = new Color(0x2F82B7);
        Color bgSample = dark ? adjust(base, -20) : adjust(base, 35);
        Color fgPrimary = pickOnColor(bgSample);

        if (headerTitleLabel != null)
            headerTitleLabel.setForeground(fgPrimary);
        if (headerSubtitleLabel != null)
            headerSubtitleLabel.setForeground(withAlpha(fgPrimary, 210));
        if (headerVersionLabel != null)
            headerVersionLabel.setForeground(withAlpha(fgPrimary, 180));
        if (darkModeToggle != null) {
            darkModeToggle.setSelected(dark);
            darkModeToggle.setText(dark ? "Light Mode" : "Dark Mode");
            darkModeToggle.setForeground(fgPrimary);
        }
        if (particleToggle != null)
            particleToggle.setForeground(fgPrimary);
        if (headerPanel != null)
            headerPanel.repaint();
    }

    // ==================== PARAMETERS PANEL ====================

    private JPanel createParametersPanel() {
        JPanel parametersPanel = new JPanel(new BorderLayout());
        parametersPanel.setBorder(BorderFactory.createEmptyBorder(15, 15, 15, 15));

        // Workflow + DB-only selectors at top of the parameters panel
        JPanel workflowPanel = new JPanel(new GridBagLayout());
        GridBagConstraints wgbc = new GridBagConstraints();
        wgbc.insets = new Insets(0, 5, 8, 5);
        wgbc.anchor = GridBagConstraints.WEST;
        wgbc.fill = GridBagConstraints.HORIZONTAL;

        // Row 0: Workflow
        wgbc.gridx = 0;
        wgbc.gridy = 0;
        wgbc.weightx = 0;
        workflowPanel.add(createLabel("Workflow:", "Select the processing workflow"), wgbc);

        String[] workflows = {
                "1. Entrapment Database Generation",
                "2. FDP Estimation / FDR Control Evaluation"
        };
        workflowCombo = new JComboBox<>(workflows);
        styleComboBox(workflowCombo);
        workflowCombo.addActionListener(e -> {
            currentWorkflow = workflowCombo.getSelectedIndex();
            updateWorkflowPanelVisibility();
        });
        wgbc.gridx = 1;
        wgbc.weightx = 1;
        workflowPanel.add(workflowCombo, wgbc);

        // Row 1: Level — single label, but the combo swaps via CardLayout
        // depending on the active workflow. DB Generation gets protein/peptide;
        // FDP Estimation gets precursor/peptide/protein/psm.
        wgbc.gridx = 0;
        wgbc.gridy = 1;
        wgbc.weightx = 0;
        workflowPanel.add(createLabel("Level:",
                "Generation/calculation level — varies by workflow"), wgbc);
        dbLevelCombo = new JComboBox<>(new String[] { "protein", "peptide" });
        styleComboBox(dbLevelCombo);
        dbLevelCombo.addActionListener(e -> updateDbLevelVisibility());
        fdpLevelCombo = new JComboBox<>(new String[] { "precursor", "peptide", "protein", "psm" });
        styleComboBox(fdpLevelCombo);
        fdpLevelCombo.addActionListener(e -> updatePepPairVisibility());
        levelCardLayout = new CardLayout();
        levelCardPanel = new JPanel(levelCardLayout);
        levelCardPanel.setOpaque(false);
        levelCardPanel.add(dbLevelCombo, "db");
        levelCardPanel.add(fdpLevelCombo, "fdp");
        wgbc.gridx = 1;
        wgbc.weightx = 1;
        workflowPanel.add(levelCardPanel, wgbc);

        // Row 2: Sequence Generation Method — shared by both workflows.
        // For DB Generation it controls Foreign-Species file visibility; for
        // FDP Estimation it indicates how the supplied entrapment was built.
        wgbc.gridx = 0;
        wgbc.gridy = 2;
        wgbc.insets = new Insets(0, 5, 15, 5);
        wgbc.weightx = 0;
        workflowPanel.add(createLabel("Sequence Generation:",
                "How entrapment sequences are produced. \"Random Shuffling\" "
                        + "shuffles target peptides; \"Foreign Species\" uses one "
                        + "or more foreign-species FASTA files as the entrapment "
                        + "source."), wgbc);
        seqMethodCombo = new JComboBox<>(new String[] { "Random Shuffling", "Foreign Species" });
        styleComboBox(seqMethodCombo);
        seqMethodCombo.addActionListener(e -> {
            updateForeignSpeciesVisibility();
            // The Enzyme / Fix N/C rows depend on the method too at protein
            // level, so re-run the digestion-row visibility logic.
            updateDbLevelVisibility();
        });
        wgbc.gridx = 1;
        wgbc.weightx = 1;
        workflowPanel.add(seqMethodCombo, wgbc);

        parametersPanel.add(workflowPanel, BorderLayout.NORTH);

        // Create workflow panels as local variables
        JPanel dbInputPanel      = createDbInputPanel();
        dbDigestionPanel = createDbDigestionPanel();
        JPanel dbOutputFormatPanel = createDbOutputFormatPanel();
        JPanel fdpInputPanel     = createFdpInputPanel();
        JPanel fdpSettingsPanel  = createFdpSettingsPanel();

        // Ensure left-alignment for BoxLayout Y_AXIS children
        dbInputPanel.setAlignmentX(Component.LEFT_ALIGNMENT);
        dbDigestionPanel.setAlignmentX(Component.LEFT_ALIGNMENT);
        dbOutputFormatPanel.setAlignmentX(Component.LEFT_ALIGNMENT);
        fdpInputPanel.setAlignmentX(Component.LEFT_ALIGNMENT);
        fdpSettingsPanel.setAlignmentX(Component.LEFT_ALIGNMENT);

        JPanel dbWorkflowStack = new JPanel();
        dbWorkflowStack.setLayout(new BoxLayout(dbWorkflowStack, BoxLayout.Y_AXIS));
        dbWorkflowStack.setAlignmentX(Component.LEFT_ALIGNMENT);
        dbWorkflowStack.add(dbInputPanel);
        dbWorkflowStack.add(Box.createVerticalStrut(10));
        dbWorkflowStack.add(dbDigestionPanel);
        dbWorkflowStack.add(Box.createVerticalStrut(10));
        dbWorkflowStack.add(dbOutputFormatPanel);

        JPanel fdpWorkflowStack = new JPanel();
        fdpWorkflowStack.setLayout(new BoxLayout(fdpWorkflowStack, BoxLayout.Y_AXIS));
        fdpWorkflowStack.setAlignmentX(Component.LEFT_ALIGNMENT);
        fdpWorkflowStack.add(fdpInputPanel);
        fdpWorkflowStack.add(Box.createVerticalStrut(10));
        fdpWorkflowStack.add(fdpSettingsPanel);

        workflowCardLayout = new CardLayout();
        workflowCardsPanel = new JPanel(workflowCardLayout) {
            @Override
            public Dimension getPreferredSize() {
                for (Component component : getComponents()) {
                    if (component.isVisible()) {
                        return component.getPreferredSize();
                    }
                }
                return super.getPreferredSize();
            }

            @Override
            public Dimension getMinimumSize() {
                return getPreferredSize();
            }
        };
        workflowCardsPanel.setAlignmentX(Component.LEFT_ALIGNMENT);
        workflowCardsPanel.add(dbWorkflowStack, "db");
        workflowCardsPanel.add(fdpWorkflowStack, "fdp");

        parametersContentPanel = new JPanel();
        parametersContentPanel.setLayout(new BoxLayout(parametersContentPanel, BoxLayout.Y_AXIS));
        parametersContentPanel.add(workflowCardsPanel);

        JScrollPane scrollPane = new JScrollPane(parametersContentPanel);
        // Use an explicit empty border (not null and not a UIResource) so the
        // LAF's installBorder() leaves it alone on theme switch — otherwise
        // updateComponentTreeUI() would reinstall the default border, causing
        // a visible flash and an extra outline around the parameters panel.
        scrollPane.setBorder(BorderFactory.createEmptyBorder());
        scrollPane.setViewportBorder(BorderFactory.createEmptyBorder());
        scrollPane.getVerticalScrollBar().setUnitIncrement(16);
        parametersPanel.add(scrollPane, BorderLayout.CENTER);

        return parametersPanel;
    }

    /**
     * Border that re-resolves {@code Component.borderColor} from
     * {@link UIManager} on every paint, so the line color tracks LAF changes
     * (light <-> dark) instead of freezing the value captured at build time.
     *
     * <p>{@code BorderFactory.createLineBorder(UIManager.getColor(...))}
     * snapshots the {@link Color} reference; after a LAF switch the
     * panels still hold the old reference, which under FlatLaf collapses to
     * black for the parameter section borders.</p>
     */
    private static javax.swing.border.Border lafEdgeBorder(
            int top, int left, int bottom, int right) {
        return new javax.swing.border.AbstractBorder() {
            @Override
            public void paintBorder(Component c, Graphics g,
                                    int x, int y, int width, int height) {
                Color color = UIManager.getColor("Component.borderColor");
                if (color == null) color = new Color(128, 128, 128);
                Color saved = g.getColor();
                g.setColor(color);
                if (top > 0)    g.fillRect(x, y, width, top);
                if (left > 0)   g.fillRect(x, y, left, height);
                if (bottom > 0) g.fillRect(x, y + height - bottom, width, bottom);
                if (right > 0)  g.fillRect(x + width - right, y, right, height);
                g.setColor(saved);
            }
            @Override
            public Insets getBorderInsets(Component c) {
                return new Insets(top, left, bottom, right);
            }
            @Override
            public Insets getBorderInsets(Component c, Insets insets) {
                insets.set(top, left, bottom, right);
                return insets;
            }
        };
    }

    private JPanel createSectionPanel(String title) {

        JPanel panel = new JPanel(new GridBagLayout());
        javax.swing.border.TitledBorder titledBorder = BorderFactory.createTitledBorder(
                lafEdgeBorder(1, 1, 1, 1), title);
        titledBorder.setTitleFont(derivedFont(Font.BOLD, 12f));
        panel.setBorder(BorderFactory.createCompoundBorder(
                titledBorder,
                BorderFactory.createEmptyBorder(4, 8, 6, 8)));
        return panel;
    }

    private JPanel createDbInputPanel() {
        JPanel panel = createSectionPanel("Input/Output Files");
        GridBagConstraints gbc = new GridBagConstraints();
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = DEFAULT_INSETS;
        int row = 0;

        gbc.gridx = 0;
        gbc.gridy = row;
        gbc.weightx = 0;
        panel.add(createLabel("Protein Database:", "Protein database file in FASTA format"), gbc);
        dbFileField = createTextField("Path to FASTA database");
        gbc.gridx = 1;
        gbc.weightx = 1;
        panel.add(dbFileField, gbc);
        gbc.gridx = 2;
        gbc.weightx = 0;
        panel.add(createDbButtonsPanel(dbFileField), gbc);
        row++;

        JLabel foreignSpeciesLabel = createLabel("Foreign Species:",
                "FASTA file(s) of foreign species used as the entrapment "
                        + "source. Use Browse to pick one or more files, or "
                        + "Folder to use every FASTA inside a directory.");
        gbc.gridx = 0;
        gbc.gridy = row;
        gbc.weightx = 0;
        panel.add(foreignSpeciesLabel, gbc);
        foreignSpeciesField = createTextField("Path(s) to foreign species FASTA files");
        gbc.gridx = 1;
        gbc.weightx = 1;
        panel.add(foreignSpeciesField, gbc);
        gbc.gridx = 2;
        gbc.weightx = 0;
        JPanel foreignSpeciesButtons = createMultiFileButtonsPanel(
                foreignSpeciesField, foreignSpeciesFiles,
                new String[] { "fasta", "fa" });
        panel.add(foreignSpeciesButtons, gbc);
        foreignSpeciesRowComponents = java.util.Arrays.asList(
                foreignSpeciesLabel, foreignSpeciesField, foreignSpeciesButtons);
        row++;

        gbc.gridx = 0;
        gbc.gridy = row;
        gbc.weightx = 0;
        panel.add(createLabel("Output Folder:",
                "Folder for the generated entrapment database. The output "
                        + "filename is auto-derived from the protein database "
                        + "name and the selected level."), gbc);
        dbOutputField = createTextField("Path to output folder");
        gbc.gridx = 1;
        gbc.weightx = 1;
        panel.add(dbOutputField, gbc);
        gbc.gridx = 2;
        gbc.weightx = 0;
        panel.add(createFolderButton(dbOutputField), gbc);

        return panel;
    }

    private static final String DEFAULT_ENZYME_NAME = "Trypsin (no P rule)";

    private static int defaultEnzymeIndex() {
        int idx = DBGear.getEnzymeNames().indexOf(DEFAULT_ENZYME_NAME);
        return idx >= 0 ? idx : 0;
    }

    private JPanel createDbDigestionPanel() {
        JPanel panel = createSectionPanel("Digestion Settings");
        GridBagConstraints gbc = new GridBagConstraints();
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = DEFAULT_INSETS;
        int row = 0;

        JLabel enzymeLabel = createLabel("Enzyme:", "Enzyme for protein digestion");
        gbc.gridx = 0;
        gbc.gridy = row;
        gbc.weightx = 0;
        panel.add(enzymeLabel, gbc);
        List<String> enzymeNames = DBGear.getEnzymeNames();
        String[] enzymeItems = new String[enzymeNames.size()];
        for (int i = 0; i < enzymeNames.size(); i++) {
            enzymeItems[i] = i + ": " + enzymeNames.get(i);
        }
        enzymeCombo = new JComboBox<>(enzymeItems);
        enzymeCombo.setSelectedIndex(defaultEnzymeIndex());
        styleComboBox(enzymeCombo);
        gbc.gridx = 1;
        gbc.weightx = 1;
        panel.add(enzymeCombo, gbc);
        row++;

        JLabel missCleavageLabel = createLabel("Missed Cleavages:", "Maximum missed cleavages");
        gbc.gridx = 0;
        gbc.gridy = row;
        gbc.weightx = 0;
        panel.add(missCleavageLabel, gbc);
        missCleavageSpinner = createSpinner(1, 0, 5, 1);
        gbc.gridx = 1;
        gbc.weightx = 1;
        panel.add(missCleavageSpinner, gbc);
        row++;

        JLabel peptideLengthLabel = createLabel("Peptide Length:", "Min and max peptide length");
        gbc.gridx = 0;
        gbc.gridy = row;
        gbc.weightx = 0;
        panel.add(peptideLengthLabel, gbc);
        // GridBag with weightx=1 on each spinner so they share col 1's width
        // equally — no trailing gap to the right of the max spinner, while
        // the min spinner still aligns flush left with the rows above/below.
        JPanel lengthPanel = new JPanel(new GridBagLayout());
        lengthPanel.setOpaque(false);
        minLengthSpinner = createSpinner(7, 1, 50, 1);
        maxLengthSpinner = createSpinner(35, 1, 100, 1);
        GridBagConstraints lgc = new GridBagConstraints();
        lgc.fill = GridBagConstraints.HORIZONTAL;
        lgc.gridy = 0;
        lgc.gridx = 0; lgc.weightx = 1; lgc.insets = new Insets(0, 0, 0, 0);
        lengthPanel.add(minLengthSpinner, lgc);
        lgc.gridx = 1; lgc.weightx = 0; lgc.insets = new Insets(0, 5, 0, 5);
        lengthPanel.add(new JLabel("-"), lgc);
        lgc.gridx = 2; lgc.weightx = 1; lgc.insets = new Insets(0, 0, 0, 0);
        lengthPanel.add(maxLengthSpinner, lgc);
        gbc.gridx = 1;
        gbc.weightx = 1;
        panel.add(lengthPanel, gbc);
        row++;

        JLabel clipNmLabel = createLabel("Clip N-terminal M:",
                "Include peptides with and without leading methionine");
        gbc.gridx = 0;
        gbc.gridy = row;
        gbc.weightx = 0;
        panel.add(clipNmLabel, gbc);
        clipNmCheckbox = new JCheckBox();
        gbc.gridx = 1;
        gbc.weightx = 1;
        panel.add(clipNmCheckbox, gbc);
        row++;

        JLabel fixNcLabel = createLabel("Fix N/C Terminal:", "Fix N and/or C terminal amino acids");
        gbc.gridx = 0;
        gbc.gridy = row;
        gbc.weightx = 0;
        panel.add(fixNcLabel, gbc);
        fixNcCombo = new JComboBox<>(new String[] { "C only", "N only", "Both N and C" });
        fixNcCombo.setSelectedIndex(0);
        styleComboBox(fixNcCombo);
        gbc.gridx = 1;
        gbc.weightx = 1;
        panel.add(fixNcCombo, gbc);

        // Rows that only ever apply at peptide level — hide whenever the user
        // is generating proteins.
        peptideOnlyRowComponents = java.util.Arrays.asList(
                missCleavageLabel, missCleavageSpinner,
                peptideLengthLabel, lengthPanel,
                clipNmLabel, clipNmCheckbox);

        // Enzyme is needed by every path that digests a protein — peptide
        // level, protein-level random shuffling, AND the protein + Foreign
        // Species + No Shared Peptides dedup check.
        digestionRowComponents = java.util.Arrays.asList(
                enzymeLabel, enzymeCombo);

        // Fix N/C Terminal is only consulted by the entrapment shuffler (it
        // pins residues during shuffling). Plain digest-for-dedup (-ns)
        // doesn't use it, so it must NOT show in that case.
        fixNcRowComponents = java.util.Arrays.asList(fixNcLabel, fixNcCombo);

        return panel;
    }

    private JPanel createDbOutputFormatPanel() {
        JPanel panel = createSectionPanel("Generation Options");
        GridBagConstraints gbc = new GridBagConstraints();
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = DEFAULT_INSETS;
        int row = 0;

        gbc.gridx = 0;
        gbc.gridy = row;
        gbc.weightx = 0;
        panel.add(createLabel("Fold:", "Number of folds for entrapment"), gbc);
        foldSpinner = createSpinner(1, 1, 10, 1);
        gbc.gridx = 1;
        gbc.weightx = 1;
        panel.add(foldSpinner, gbc);
        row++;

        gbc.gridx = 0;
        gbc.gridy = row;
        gbc.weightx = 0;
        panel.add(createLabel("Convert I->L:", "Convert isoleucine to leucine"), gbc);
        i2lCheckbox = new JCheckBox();
        i2lCheckbox.setSelected(true);
        gbc.gridx = 1;
        gbc.weightx = 1;
        panel.add(i2lCheckbox, gbc);
        row++;

        gbc.gridx = 0;
        gbc.gridy = row;
        gbc.weightx = 0;
        panel.add(createLabel("DIA-NN Format:", "Output compatible with DIA-NN"), gbc);
        diannCheckbox = new JCheckBox();
        diannCheckbox.setSelected(true);
        gbc.gridx = 1;
        gbc.weightx = 1;
        panel.add(diannCheckbox, gbc);
        row++;

        gbc.gridx = 0;
        gbc.gridy = row;
        gbc.weightx = 0;
        panel.add(createLabel("UniProt Format:", "Use UniProt protein ID format"), gbc);
        uniprotCheckbox = new JCheckBox();
        uniprotCheckbox.setSelected(true);
        gbc.gridx = 1;
        gbc.weightx = 1;
        panel.add(uniprotCheckbox, gbc);
        row++;

        JLabel exportDbLbl = createLabel("Export Protein DB:",
                "Also write a peptide-grouped protein FASTA next to the "
                        + "main output (peptide-level only — at protein level "
                        + "the main output is already a protein FASTA).");
        gbc.gridx = 0;
        gbc.gridy = row;
        gbc.weightx = 0;
        panel.add(exportDbLbl, gbc);
        exportDbCheckbox = new JCheckBox();
        gbc.gridx = 1;
        gbc.weightx = 1;
        panel.add(exportDbCheckbox, gbc);
        exportDbRowComponents = java.util.Arrays.asList(exportDbLbl, exportDbCheckbox);
        row++;

        JLabel checkDupLbl = createLabel("Check Duplicates:",
                "For protein-level random shuffling: when the "
                        + "global-uniqueness shuffle can't produce a valid "
                        + "entrapment for a peptide, select a random peptide "
                        + "without the dedup constraint.");
        gbc.gridx = 0;
        gbc.gridy = row;
        gbc.weightx = 0;
        panel.add(checkDupLbl, gbc);
        checkDuplicatesCheckbox = new JCheckBox();
        // Default ON — matches the README example for protein-level random
        // shuffling (`-check`), which is the only path that honors it.
        checkDuplicatesCheckbox.setSelected(true);
        gbc.gridx = 1;
        gbc.weightx = 1;
        panel.add(checkDuplicatesCheckbox, gbc);
        checkDupRowComponents = java.util.Arrays.asList(checkDupLbl, checkDuplicatesCheckbox);
        row++;

        JLabel noSharedLbl = createLabel("No Shared Peptides:",
                "Drop any foreign-species protein whose digestion produces a "
                        + "peptide that already exists in the target database, "
                        + "so the entrapment set is fully disjoint from the "
                        + "targets (protein-level + Foreign Species only).");
        gbc.gridx = 0;
        gbc.gridy = row;
        gbc.weightx = 0;
        panel.add(noSharedLbl, gbc);
        noSharedCheckbox = new JCheckBox();
        // Toggling -ns at protein + foreign species needs to reveal/hide the
        // Digestion section, since the dedup check digests both target and
        // entrapment proteins and uses Enzyme / Peptide Length / etc.
        noSharedCheckbox.addActionListener(e -> updateDbLevelVisibility());
        gbc.gridx = 1;
        gbc.weightx = 1;
        panel.add(noSharedCheckbox, gbc);
        noSharedRowComponents = java.util.Arrays.asList(noSharedLbl, noSharedCheckbox);
        row++;

        // -------- Formerly the Advanced Settings section, merged in here. --

        gbc.gridx = 0;
        gbc.gridy = row;
        gbc.weightx = 0;
        panel.add(createLabel("Random Seed:", "Random seed for reproducibility"), gbc);
        seedSpinner = createSpinner(2000, 0, 999999, 1);
        gbc.gridx = 1;
        gbc.weightx = 1;
        panel.add(seedSpinner, gbc);
        row++;

        // Entrapment Label + Pos share a row: text field grows, inline
        // "Pos:" label and the position dropdown sit flush right.
        gbc.gridx = 0;
        gbc.gridy = row;
        gbc.weightx = 0;
        panel.add(createLabel("Entrapment Label:", "Label for entrapment peptides"), gbc);
        entrapmentLabelField = createTextField("_p_target");
        entrapmentLabelField.setText("_p_target");
        entrapmentPosCombo = new JComboBox<>(new String[] { "End", "Start" });
        entrapmentPosCombo.setSelectedIndex(0);
        styleComboBox(entrapmentPosCombo);
        // Field and combo share col 1 with equal widths: 3-cell GridBag,
        // weightx=1 on both, preferred widths forced to 0 so the weightx
        // distribution dominates (otherwise the field's 200-px default
        // preferred width would make it noticeably wider than the combo).
        entrapmentLabelField.setPreferredSize(new Dimension(0,
                entrapmentLabelField.getPreferredSize().height));
        entrapmentPosCombo.setPreferredSize(new Dimension(0,
                entrapmentPosCombo.getPreferredSize().height));
        JPanel entrapmentRow = new JPanel(new GridBagLayout());
        entrapmentRow.setOpaque(false);
        GridBagConstraints egc = new GridBagConstraints();
        egc.fill = GridBagConstraints.HORIZONTAL;
        egc.gridy = 0;
        egc.gridx = 0; egc.weightx = 1; egc.insets = new Insets(0, 0, 0, 8);
        entrapmentRow.add(entrapmentLabelField, egc);
        egc.gridx = 1; egc.weightx = 0; egc.insets = new Insets(0, 0, 0, 4);
        entrapmentRow.add(createLabel("Label Position:",
                "Position of the entrapment label"), egc);
        egc.gridx = 2; egc.weightx = 1; egc.insets = new Insets(0, 0, 0, 0);
        entrapmentRow.add(entrapmentPosCombo, egc);
        gbc.gridx = 1;
        gbc.weightx = 1;
        panel.add(entrapmentRow, gbc);
        row++;

        // Add Decoys + Decoy Label + Pos share one row. When the checkbox is
        // unchecked, the inline decoy label/pos controls hide themselves but
        // the row (checkbox) remains visible. The whole row is itself only
        // meaningful for Random Shuffling generation — Foreign Species runs
        // don't generate decoys for the entrapment proteins they import.
        JLabel addDecoysLbl = createLabel("Add Decoys:",
                "Add decoy sequences to the generated output (Random "
                        + "Shuffling generation only).");
        decoyCheckbox = new JCheckBox();
        decoyLabelField = createTextField("rev_");
        decoyLabelField.setText("rev_");
        decoyPosCombo = new JComboBox<>(new String[] { "Start", "End" });
        decoyPosCombo.setSelectedIndex(0);
        styleComboBox(decoyPosCombo);
        JLabel decoyLabelLbl = createLabel("Decoy Label:", "Label for decoy sequences");
        JLabel decoyPosLbl   = createLabel("Label Position:", "Position of the decoy label");

        // Pin checkbox height so the row height stays constant whether the
        // decoy details are visible or hidden — otherwise the row would
        // collapse to the checkbox's smaller natural height.
        Dimension cbSize = decoyCheckbox.getPreferredSize();
        decoyCheckbox.setPreferredSize(new Dimension(
                cbSize.width, decoyLabelField.getPreferredSize().height));
        // Force field/combo preferred widths to 0 so weightx=1 distribution
        // drives the actual widths — both halves end up equal regardless of
        // natural preferred widths.
        decoyLabelField.setPreferredSize(new Dimension(0,
                decoyLabelField.getPreferredSize().height));
        decoyPosCombo.setPreferredSize(new Dimension(0,
                decoyPosCombo.getPreferredSize().height));

        // Left half = [checkbox][Decoy Label:][field]. Checkbox lives here so
        // it sits flush at the very left of col 1 (aligned with the other
        // panel checkboxes), and stays visible even when the rest of the
        // decoy details are hidden.
        //
        // Use nested BorderLayouts (not GridBag): when the inner controls
        // are hidden, GridBagLayout has no weighted visible cells and
        // centers the lone visible component within its container — which
        // pushed the checkbox to the middle of col 1. BorderLayout never
        // auto-centers, so WEST stays anchored to the left.
        JPanel decoyInner = new JPanel(new BorderLayout(4, 0));
        decoyInner.setOpaque(false);
        decoyInner.add(decoyLabelLbl, BorderLayout.WEST);
        decoyInner.add(decoyLabelField, BorderLayout.CENTER);

        JPanel decoyLeftHalf = new JPanel(new BorderLayout(8, 0));
        decoyLeftHalf.setOpaque(false);
        decoyLeftHalf.add(decoyCheckbox, BorderLayout.WEST);
        decoyLeftHalf.add(decoyInner, BorderLayout.CENTER);
        // Override the half's preferred width to 0 so the outer GridBag's
        // weightx=1 distribution makes both halves the same width.
        decoyLeftHalf.setPreferredSize(new Dimension(0,
                decoyLabelField.getPreferredSize().height));

        // Outer 3-cell GridBag: [leftHalf | "Label Position:" | combo].
        // Equal weightx and 0 preferreds → identical halves; "Label Position:"
        // sits exactly at the visible seam, mirroring the entrapment row.
        JPanel decoyRowControls = new JPanel(new GridBagLayout());
        decoyRowControls.setOpaque(false);
        GridBagConstraints ddgc = new GridBagConstraints();
        ddgc.fill = GridBagConstraints.HORIZONTAL;
        ddgc.gridy = 0;
        ddgc.gridx = 0; ddgc.weightx = 1; ddgc.insets = new Insets(0, 0, 0, 8);
        decoyRowControls.add(decoyLeftHalf, ddgc);
        ddgc.gridx = 1; ddgc.weightx = 0; ddgc.insets = new Insets(0, 0, 0, 4);
        decoyRowControls.add(decoyPosLbl, ddgc);
        ddgc.gridx = 2; ddgc.weightx = 1; ddgc.insets = new Insets(0, 0, 0, 0);
        decoyRowControls.add(decoyPosCombo, ddgc);

        gbc.gridx = 0;
        gbc.gridy = row;
        gbc.weightx = 0;
        panel.add(addDecoysLbl, gbc);
        gbc.gridx = 1;
        gbc.weightx = 1;
        panel.add(decoyRowControls, gbc);
        row++;

        addDecoysRowComponents = java.util.Arrays.asList(addDecoysLbl, decoyRowControls);
        // Hide just the inner controls (label/field/pos label/combo) when
        // Add Decoys is unchecked. Checkbox is NOT in this list so it stays
        // visible at the left of leftHalf.
        decoyLabelRowComponents = java.util.Arrays.asList(
                decoyLabelLbl, decoyLabelField, decoyPosLbl, decoyPosCombo);
        decoyCheckbox.addActionListener(e -> updateDecoyLabelVisibility());
        updateDecoyLabelVisibility();

        gbc.gridx = 0;
        gbc.gridy = row;
        gbc.weightx = 0;
        panel.add(createLabel("Debug Logging:", "Print detailed information for debugging"), gbc);
        debugCheckbox = new JCheckBox();
        gbc.gridx = 1;
        gbc.weightx = 1;
        panel.add(debugCheckbox, gbc);

        return panel;
    }

    private JPanel createFdpInputPanel() {
        JPanel panel = createSectionPanel("Input/Output Files");
        GridBagConstraints gbc = new GridBagConstraints();
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = DEFAULT_INSETS;
        int row = 0;

        // Input file
        gbc.gridx = 0;
        gbc.gridy = row;
        gbc.weightx = 0;
        panel.add(createLabel("Input File:", "PSM/peptide/precursor/protein file"), gbc);
        inputFileField = createTextField("Path to input file");
        // Refresh the Score Column dropdown whenever the input file path
        // changes — pulls numeric columns from the file's first data row so
        // the user can pick from a real list instead of typing a name.
        inputFileField.getDocument().addDocumentListener(new javax.swing.event.DocumentListener() {
            private void refresh() {
                String path = inputFileField.getText().trim();
                File f = path.isEmpty() ? null : new File(path);
                refreshScoreColumns(f != null && f.isFile() ? f : null);
            }
            @Override public void insertUpdate(javax.swing.event.DocumentEvent e) { refresh(); }
            @Override public void removeUpdate(javax.swing.event.DocumentEvent e) { refresh(); }
            @Override public void changedUpdate(javax.swing.event.DocumentEvent e) { refresh(); }
        });
        gbc.gridx = 1;
        gbc.weightx = 1;
        panel.add(inputFileField, gbc);
        gbc.gridx = 2;
        gbc.weightx = 0;
        panel.add(createBrowseAndViewPanel(inputFileField, "tsv,csv,txt"), gbc);
        row++;

        // Peptide pair file (hidden when level = protein — pair tracking is
        // peptide/precursor/PSM specific).
        JLabel pepPairLabel = createLabel("Peptide Pair File:",
                "Optional. Peptide pair file used for "
                        + "peptide / precursor / PSM-level FDP. Leave empty "
                        + "to skip pair tracking.");
        gbc.gridx = 0;
        gbc.gridy = row;
        gbc.weightx = 0;
        panel.add(pepPairLabel, gbc);
        pepPairFileField = createTextField(
                "Optional — path to peptide pair file");
        gbc.gridx = 1;
        gbc.weightx = 1;
        panel.add(pepPairFileField, gbc);
        gbc.gridx = 2;
        gbc.weightx = 0;
        JPanel pepPairButtons = createBrowseAndViewPanel(pepPairFileField, "tsv,csv,txt");
        panel.add(pepPairButtons, gbc);
        pepPairRowComponents = java.util.Arrays.asList(
                pepPairLabel, pepPairFileField, pepPairButtons);
        row++;

        // Output folder
        gbc.gridx = 0;
        gbc.gridy = row;
        gbc.weightx = 0;
        panel.add(createLabel("Output Folder:",
                "Folder for the FDP estimation result. The output filename "
                        + "is auto-derived from the input file name and the "
                        + "selected level."), gbc);
        fdpOutputField = createTextField("Path to output folder");
        gbc.gridx = 1;
        gbc.weightx = 1;
        panel.add(fdpOutputField, gbc);
        gbc.gridx = 2;
        gbc.weightx = 0;
        panel.add(createFolderButton(fdpOutputField), gbc);

        return panel;
    }

    private JPanel createFdpSettingsPanel() {
        JPanel panel = createSectionPanel("FDP Calculation Settings");
        GridBagConstraints gbc = new GridBagConstraints();
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = DEFAULT_INSETS;
        int row = 0;

        // Score Column + Score Direction share a row: the column dropdown
        // grows; "Direction:" label and direction dropdown sit flush right.
        gbc.gridx = 0;
        gbc.gridy = row;
        gbc.weightx = 0;
        panel.add(createLabel("Score Column:",
                "Column used for ranking. Auto-populated from the input "
                        + "file's numeric columns; \"none\" ranks by q_value "
                        + "alone."), gbc);
        scoreColumnCombo = new JComboBox<>(new String[] { NO_SCORE_ITEM });
        scoreColumnCombo.setEditable(false);
        scoreColumnCombo.setSelectedItem(NO_SCORE_ITEM);
        styleComboBox(scoreColumnCombo);
        scoreDirectionCombo = new JComboBox<>(new String[] { "Lower is better", "Higher is better" });
        styleComboBox(scoreDirectionCombo);
        // 3-cell GridBag: combo | "Direction:" | combo. Both combos use
        // weightx=1 AND have their preferred widths forced to 0 so the
        // weightx distribution drives the final widths — they end up equal
        // regardless of differences in the items' natural widths.
        scoreColumnCombo.setPreferredSize(new Dimension(0,
                scoreColumnCombo.getPreferredSize().height));
        scoreDirectionCombo.setPreferredSize(new Dimension(0,
                scoreDirectionCombo.getPreferredSize().height));
        JLabel scoreDirInlineLbl = createLabel("Direction:",
                "Whether lower or higher scores are better");
        JPanel scoreRow = new JPanel(new GridBagLayout());
        scoreRow.setOpaque(false);
        GridBagConstraints sgc = new GridBagConstraints();
        sgc.fill = GridBagConstraints.HORIZONTAL;
        sgc.gridy = 0;
        sgc.gridx = 0; sgc.weightx = 1; sgc.insets = new Insets(0, 0, 0, 8);
        scoreRow.add(scoreColumnCombo, sgc);
        sgc.gridx = 1; sgc.weightx = 0; sgc.insets = new Insets(0, 0, 0, 4);
        scoreRow.add(scoreDirInlineLbl, sgc);
        sgc.gridx = 2; sgc.weightx = 1; sgc.insets = new Insets(0, 0, 0, 0);
        scoreRow.add(scoreDirectionCombo, sgc);
        gbc.gridx = 1;
        gbc.weightx = 1;
        panel.add(scoreRow, gbc);
        row++;

        JLabel fdpFoldLabel = createLabel("Fold:",
                "Number of entrapment folds (Random Shuffling only — for "
                        + "Foreign Species use the R Ratio instead).");
        gbc.gridx = 0;
        gbc.gridy = row;
        gbc.weightx = 0;
        panel.add(fdpFoldLabel, gbc);
        fdpFoldSpinner = createSpinner(1, 1, 10, 1);
        gbc.gridx = 1;
        gbc.weightx = 1;
        panel.add(fdpFoldSpinner, gbc);
        fdpFoldRowComponents = java.util.Arrays.asList(fdpFoldLabel, fdpFoldSpinner);
        row++;

        JLabel rRatioLabel = createLabel("R Ratio:",
                "#entrapment/#target ratio for combined entrapment "
                        + "(used when entrapments come from foreign species)");
        gbc.gridx = 0;
        gbc.gridy = row;
        gbc.weightx = 0;
        panel.add(rRatioLabel, gbc);
        rRatioField = createTextField("Optional — numeric, e.g. 1.0");
        gbc.gridx = 1;
        gbc.weightx = 1;
        panel.add(rRatioField, gbc);
        rRatioRowComponents = java.util.Arrays.asList(rRatioLabel, rRatioField);
        row++;

        gbc.gridx = 0;
        gbc.gridy = row;
        gbc.weightx = 0;
        panel.add(createLabel("Pick Method:", "How to pick one representative protein from a protein group.\n" +
                "Proteins in the same protein group are separated by ';'.\n"+
                "This is only used when (1) determine protein type (target or entrapment) for evaluating protein level FDR control;\n" +
                "(2) determine peptide type (target or entrapment) for evaluating \n" +
                "PSM/precursor/peptide level FDR control without \n" +
                "providing peptide type mapping file."), gbc);
        pickMethodCombo = new JComboBox<>(new String[] { "first", "last", "random", "as_is" });
        styleComboBox(pickMethodCombo);
        gbc.gridx = 1;
        gbc.weightx = 1;
        panel.add(pickMethodCombo, gbc);
        row++;

        // Entrapment Label + Pos share a row, mirroring the DB-gen panel.
        gbc.gridx = 0;
        gbc.gridy = row;
        gbc.weightx = 0;
        panel.add(createLabel("Entrapment Label:",
                "Suffix/prefix that marks entrapment entries in the input "
                        + "file's protein column. Must match the label used "
                        + "when generating the entrapment database."), gbc);
        fdpEntrapmentLabelField = createTextField("_p_target");
        fdpEntrapmentLabelField.setText("_p_target");
        fdpEntrapmentPosCombo = new JComboBox<>(new String[] { "End", "Start" });
        fdpEntrapmentPosCombo.setSelectedIndex(0);
        styleComboBox(fdpEntrapmentPosCombo);
        // Field and combo share col 1 with equal widths — same trick as the
        // DB-gen entrapment row.
        fdpEntrapmentLabelField.setPreferredSize(new Dimension(0,
                fdpEntrapmentLabelField.getPreferredSize().height));
        fdpEntrapmentPosCombo.setPreferredSize(new Dimension(0,
                fdpEntrapmentPosCombo.getPreferredSize().height));
        JLabel labelPosInlineLbl = createLabel("Label Position:",
                "Position of the entrapment label on the protein ID.");
        // Align the two combined rows: pad the narrower inline label so its
        // width matches the wider one ("Label Position:" vs "Direction:").
        // Without this, the right-hand combos start at different X positions
        // on the Score row and the Entrapment row.
        int inlineLabelWidth = Math.max(
                scoreDirInlineLbl.getPreferredSize().width,
                labelPosInlineLbl.getPreferredSize().width);
        scoreDirInlineLbl.setPreferredSize(new Dimension(inlineLabelWidth,
                scoreDirInlineLbl.getPreferredSize().height));
        labelPosInlineLbl.setPreferredSize(new Dimension(inlineLabelWidth,
                labelPosInlineLbl.getPreferredSize().height));
        JPanel fdpEntrapmentRow = new JPanel(new GridBagLayout());
        fdpEntrapmentRow.setOpaque(false);
        GridBagConstraints fgc = new GridBagConstraints();
        fgc.fill = GridBagConstraints.HORIZONTAL;
        fgc.gridy = 0;
        fgc.gridx = 0; fgc.weightx = 1; fgc.insets = new Insets(0, 0, 0, 8);
        fdpEntrapmentRow.add(fdpEntrapmentLabelField, fgc);
        fgc.gridx = 1; fgc.weightx = 0; fgc.insets = new Insets(0, 0, 0, 4);
        fdpEntrapmentRow.add(labelPosInlineLbl, fgc);
        fgc.gridx = 2; fgc.weightx = 1; fgc.insets = new Insets(0, 0, 0, 0);
        fdpEntrapmentRow.add(fdpEntrapmentPosCombo, fgc);
        gbc.gridx = 1;
        gbc.weightx = 1;
        panel.add(fdpEntrapmentRow, gbc);

        return panel;
    }

    private void updateWorkflowPanelVisibility() {
        boolean isDbGen = (currentWorkflow == WORKFLOW_DB_GENERATION);
        workflowCardLayout.show(workflowCardsPanel, isDbGen ? "db" : "fdp");
        // Swap the Level combo to the one that fits the active workflow.
        if (levelCardLayout != null && levelCardPanel != null) {
            levelCardLayout.show(levelCardPanel, isDbGen ? "db" : "fdp");
        }
        updateForeignSpeciesVisibility();
        updateDbLevelVisibility();
        updatePlotTabVisibility();
        // Revalidate up the chain so the scroll pane re-measures the current
        // card's preferred size and collapses the blank space beneath it.
        workflowCardsPanel.revalidate();
        workflowCardsPanel.repaint();
        parametersContentPanel.revalidate();
        parametersContentPanel.repaint();
        Container viewport = workflowCardsPanel.getParent();
        if (viewport != null) {
            viewport.revalidate();
            viewport.repaint();
        }
    }

    /**
     * Add the Plot tab when the FDP Estimation workflow is selected, remove
     * it when DB Generation is selected. The tab sits between Workflow and
     * Console, so its index is always the second tab (index 1) when present.
     */
    private void updatePlotTabVisibility() {
        if (tabbedPane == null || fdpPlotPanel == null) {
            return;
        }
        int existing = tabbedPane.indexOfComponent(fdpPlotPanel);
        boolean shouldShow = (currentWorkflow == WORKFLOW_FDP_ESTIMATION);
        if (shouldShow && existing < 0) {
            // Insert before the Console tab so order is Workflow / Plot / Console.
            int consoleIdx = tabbedPane.indexOfTab("Console");
            int insertAt = consoleIdx >= 0 ? consoleIdx : tabbedPane.getTabCount();
            tabbedPane.insertTab(PLOT_TAB_TITLE, null, fdpPlotPanel,
                    "FDP plot for the most recent FDP-Estimation run", insertAt);
        } else if (!shouldShow && existing >= 0) {
            tabbedPane.remove(existing);
        }
    }

    private void updateForeignSpeciesVisibility() {
        if (foreignSpeciesRowComponents == null || seqMethodCombo == null) {
            return;
        }
        boolean show = (currentWorkflow == WORKFLOW_DB_GENERATION)
                && seqMethodCombo.getSelectedIndex() == SEQ_METHOD_FOREIGN;
        for (JComponent c : foreignSpeciesRowComponents) {
            c.setVisible(show);
        }
        Container parent = foreignSpeciesRowComponents.get(0).getParent();
        if (parent != null) {
            parent.revalidate();
            parent.repaint();
        }
        updateRRatioVisibility();
        updateFdpFoldVisibility();
        updatePepPairVisibility();
        updateCheckDupVisibility();
        updateNoSharedVisibility();
        updateAddDecoysVisibility();
    }

    /**
     * Show the Check Duplicates row only when the only generation path that
     * actually honors -check is selected — protein level + random shuffling.
     * In every other workflow combination the flag is a no-op (verified by
     * grepping FDREval.java: it's only referenced inside
     * generate_protein_entrapment_database).
     */
    private void updateCheckDupVisibility() {
        if (checkDupRowComponents == null || dbLevelCombo == null) {
            return;
        }
        boolean isProtein   = "protein".equals(dbLevelCombo.getSelectedItem());
        boolean isShuffling = seqMethodCombo == null
                || seqMethodCombo.getSelectedIndex() == SEQ_METHOD_SHUFFLING;
        boolean show = isProtein && isShuffling
                && currentWorkflow == WORKFLOW_DB_GENERATION;
        for (JComponent c : checkDupRowComponents) {
            c.setVisible(show);
        }
        Container parent = checkDupRowComponents.get(0).getParent();
        if (parent != null) {
            parent.revalidate();
            parent.repaint();
        }
    }

    /**
     * Show the No Shared Peptides row only when the only generation path
     * that actually honors -ns is selected — protein level + foreign species
     * (only generate_protein_entrapment_database_from_multiple_species_data
     * references it).
     */
    private void updateNoSharedVisibility() {
        if (noSharedRowComponents == null || dbLevelCombo == null
                || seqMethodCombo == null) {
            return;
        }
        boolean isProtein = "protein".equals(dbLevelCombo.getSelectedItem());
        boolean isForeign = seqMethodCombo.getSelectedIndex() == SEQ_METHOD_FOREIGN;
        boolean show = isProtein && isForeign
                && currentWorkflow == WORKFLOW_DB_GENERATION;
        for (JComponent c : noSharedRowComponents) {
            c.setVisible(show);
        }
        Container parent = noSharedRowComponents.get(0).getParent();
        if (parent != null) {
            parent.revalidate();
            parent.repaint();
        }
    }

    /**
     * Hide the Decoy Label / Decoy Label Pos rows when Add Decoys is off,
     * when the row itself is hidden (Foreign Species method), or when we're
     * in the FDP workflow.
     */
    private void updateDecoyLabelVisibility() {
        if (decoyLabelRowComponents == null || decoyCheckbox == null) {
            return;
        }
        boolean show = decoyCheckbox.isSelected() && isAddDecoysApplicable();
        for (JComponent c : decoyLabelRowComponents) {
            c.setVisible(show);
        }
        Container parent = decoyLabelRowComponents.get(0).getParent();
        if (parent != null) {
            parent.revalidate();
            parent.repaint();
        }
    }

    /**
     * Hide the whole Add Decoys row outside DB Generation + Random Shuffling.
     * Foreign-species runs don't decoy the imported entrapment proteins, so
     * the {@code -decoy} flag is a no-op in that context.
     */
    private void updateAddDecoysVisibility() {
        if (addDecoysRowComponents == null) {
            return;
        }
        boolean show = isAddDecoysApplicable();
        for (JComponent c : addDecoysRowComponents) {
            c.setVisible(show);
        }
        Container parent = addDecoysRowComponents.get(0).getParent();
        if (parent != null) {
            parent.revalidate();
            parent.repaint();
        }
        // Cascade — the decoy-label rows depend on this row too.
        updateDecoyLabelVisibility();
    }

    private boolean isAddDecoysApplicable() {
        boolean isShuffling = seqMethodCombo == null
                || seqMethodCombo.getSelectedIndex() == SEQ_METHOD_SHUFFLING;
        return currentWorkflow == WORKFLOW_DB_GENERATION && isShuffling;
    }

    /**
     * Hide the Peptide Pair File row in the FDP workflow when the calculation
     * level is "protein" — the pair file tracks peptide-to-target pairings
     * which has no meaning for protein-level FDR control.
     */
    private void updatePepPairVisibility() {
        if (pepPairRowComponents == null || fdpLevelCombo == null) {
            return;
        }
        boolean show = currentWorkflow == WORKFLOW_FDP_ESTIMATION
                && !"protein".equals(fdpLevelCombo.getSelectedItem());
        for (JComponent c : pepPairRowComponents) {
            c.setVisible(show);
        }
        Container parent = pepPairRowComponents.get(0).getParent();
        if (parent != null) {
            parent.revalidate();
            parent.repaint();
        }
    }

    /**
     * In the FDP workflow, Fold and R Ratio are mutually exclusive:
     * <ul>
     *   <li>Random Shuffling → show <b>Fold</b> (k-fold FDP), hide R Ratio.</li>
     *   <li>Foreign Species → show <b>R Ratio</b> (combined-entrapment ratio),
     *       hide Fold.</li>
     * </ul>
     * In other workflows both rows stay hidden.
     */
    private void updateFdpFoldVisibility() {
        if (fdpFoldRowComponents == null) {
            return;
        }
        boolean show = currentWorkflow == WORKFLOW_FDP_ESTIMATION
                && (seqMethodCombo == null
                        || seqMethodCombo.getSelectedIndex() == SEQ_METHOD_SHUFFLING);
        for (JComponent c : fdpFoldRowComponents) {
            c.setVisible(show);
        }
        Container parent = fdpFoldRowComponents.get(0).getParent();
        if (parent != null) {
            parent.revalidate();
            parent.repaint();
        }
    }

    /**
     * Hide the R Ratio row in the FDP workflow when Sequence Generation is
     * Random Shuffling — R is the #entrapment/#target ratio used for combined
     * entrapment with foreign species, so it has no meaning for shuffling.
     */
    private void updateRRatioVisibility() {
        if (rRatioRowComponents == null || seqMethodCombo == null) {
            return;
        }
        boolean show = (currentWorkflow == WORKFLOW_FDP_ESTIMATION)
                && seqMethodCombo.getSelectedIndex() == SEQ_METHOD_FOREIGN;
        for (JComponent c : rRatioRowComponents) {
            c.setVisible(show);
        }
        Container parent = rRatioRowComponents.get(0).getParent();
        if (parent != null) {
            parent.revalidate();
            parent.repaint();
        }
    }

    /**
     * Toggle the digestion-related rows in the DB Generation workflow.
     *
     * <ul>
     *   <li>Pure peptide-only rows (Missed Cleavages, Peptide Length, Clip
     *       N-terminal M) hide whenever Level = protein.</li>
     *   <li>Digestion rows (Enzyme, Fix N/C Terminal) hide only at protein
     *       level with Sequence Generation = Foreign Species — protein-level
     *       random shuffling still digests the target before shuffling.</li>
     * </ul>
     */
    private void updateDbLevelVisibility() {
        if (dbLevelCombo == null) {
            return;
        }
        boolean isPeptide = "peptide".equals(dbLevelCombo.getSelectedItem());
        boolean isShuffling = seqMethodCombo == null
                || seqMethodCombo.getSelectedIndex() == SEQ_METHOD_SHUFFLING;
        // At protein + Foreign Species, digestion settings are normally
        // hidden — but if -ns is enabled the algorithm digests both target
        // and entrapment proteins to compare peptides, so the user needs to
        // see (and tune) Enzyme / Missed Cleavages / Peptide Length / Clip
        // N-term M / Fix N/C Terminal in that case.
        boolean nsActive = !isPeptide && !isShuffling
                && currentWorkflow == WORKFLOW_DB_GENERATION
                && noSharedCheckbox != null && noSharedCheckbox.isSelected();
        boolean showPeptideOnly = isPeptide || nsActive;
        boolean showDigestion   = isPeptide || isShuffling || nsActive;

        if (peptideOnlyRowComponents != null) {
            for (JComponent c : peptideOnlyRowComponents) {
                c.setVisible(showPeptideOnly);
            }
        }
        if (digestionRowComponents != null) {
            for (JComponent c : digestionRowComponents) {
                c.setVisible(showDigestion);
            }
        }
        // Fix N/C Terminal is consulted ONLY inside the entrapment shuffler
        // (it pins residues during shuffle/swap). At peptide + Foreign
        // Species the only place fix_c/fix_n appear is inside an
        // `if (add_decoy)` block, and Add Decoys is itself hidden in that
        // combination — so the option has no effect there. Bottom line: show
        // the row only when method = Random Shuffling.
        if (fixNcRowComponents != null) {
            for (JComponent c : fixNcRowComponents) {
                c.setVisible(isShuffling);
            }
        }
        // Export Protein DB only matters at peptide level — at protein level
        // the main -o output already IS a protein FASTA.
        if (exportDbRowComponents != null) {
            for (JComponent c : exportDbRowComponents) {
                c.setVisible(isPeptide);
            }
        }
        // Check Duplicates only takes effect in protein-level random
        // shuffling (FDREval.generate_protein_entrapment_database).
        updateCheckDupVisibility();
        // No Shared Peptides only takes effect in protein-level foreign
        // species (FDREval.generate_protein_entrapment_database_from_multiple_species_data).
        updateNoSharedVisibility();
        // Hide the whole Digestion Settings section when every row inside is
        // hidden (the protein + Foreign Species combination), so the empty
        // titled border doesn't dangle.
        if (dbDigestionPanel != null) {
            dbDigestionPanel.setVisible(isPeptide || showDigestion);
            Container parent = dbDigestionPanel.getParent();
            if (parent != null) {
                parent.revalidate();
                parent.repaint();
            }
        }
    }

    // ==================== CONSOLE PANEL ====================

    private JPanel createConsolePanel() {
        JPanel panel = new JPanel(new BorderLayout());
        panel.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));

        JPanel headerWrapper = new JPanel(new BorderLayout());
        headerWrapper.setOpaque(false);
        headerWrapper.setBorder(BorderFactory.createEmptyBorder(0, 0, 8, 0));

        JLabel consoleLabel = new JLabel("[>] Console Output");
        consoleLabel.setFont(derivedFont(Font.BOLD, 13f));
        headerWrapper.add(consoleLabel, BorderLayout.WEST);

        JButton copyButton = new JButton("Copy Output");
        styleButton(copyButton);
        copyButton.addActionListener(e -> {
            String content = consoleArea.getText();
            if (content != null && !content.isEmpty()) {
                java.awt.datatransfer.StringSelection sel = new java.awt.datatransfer.StringSelection(content);
                Toolkit.getDefaultToolkit().getSystemClipboard().setContents(sel, null);
                copyButton.setText("Copied!");
                new javax.swing.Timer(1500, ev -> copyButton.setText("Copy Output")).start();
            }
        });
        headerWrapper.add(copyButton, BorderLayout.EAST);
        panel.add(headerWrapper, BorderLayout.NORTH);

        consoleArea = new JTextArea();
        consoleArea.setEditable(false);
        // line wrapping is off so high-volume process output
        // doesn't trigger constant relayout while running.
        consoleArea.setLineWrap(false);
        consoleArea.setWrapStyleWord(false);

        JScrollPane consoleScrollPane = new JScrollPane(consoleArea);
        consoleScrollPane.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
        consoleScrollPane.setBorder(lafEdgeBorder(1, 1, 1, 1));
        panel.add(consoleScrollPane, BorderLayout.CENTER);

        return panel;
    }

    // ==================== FOOTER ====================

    private JPanel createFooter() {
        JPanel footer = new JPanel(new BorderLayout());
        footer.setBorder(lafEdgeBorder(1, 0, 0, 0));

        progressBar = new JProgressBar();
        progressBar.setIndeterminate(false);
        progressBar.setStringPainted(true);
        progressBar.setString("Ready");
        progressBar.setFont(derivedFont(Font.PLAIN, 11f));
        progressBar.setBorder(BorderFactory.createEmptyBorder(6, 12, 6, 12));
        footer.add(progressBar, BorderLayout.NORTH);

        // 6 buttons (Run / Stop / Preview / Clear / Reset / Help) total wider
        // than the default 700px window's content area, so they wrap to a
        // second row. Plain FlowLayout reports one-row preferred height even
        // when laying out two rows, which makes BorderLayout.CENTER clip the
        // wrapped row (the Help button disappears). WrapLayout reports the
        // actual wrapped height so the row stays visible.
        JPanel buttonsPanel = new JPanel(new WrapLayout(FlowLayout.CENTER, 6, 8));

        runButton = createPrimaryButton("Run FDRBench", ACCENT_COLOR);
        runButton.setMnemonic(java.awt.event.KeyEvent.VK_R);
        runButton.addActionListener(e -> runFDRBench());
        buttonsPanel.add(runButton);

        stopButton = createPrimaryButton("Stop", new Color(231, 76, 60));
        stopButton.setEnabled(false);
        stopButton.addActionListener(e -> stopFDRBench());
        buttonsPanel.add(stopButton);

        JButton previewButton = createSecondaryButton("Preview Command");
        previewButton.addActionListener(e -> previewCommand());
        buttonsPanel.add(previewButton);

        JButton clearButton = createSecondaryButton("Clear Console");
        clearButton.addActionListener(e -> consoleArea.setText(""));
        buttonsPanel.add(clearButton);

        JButton resetButton = createSecondaryButton("Reset Defaults");
        resetButton.addActionListener(e -> resetToDefaults());
        buttonsPanel.add(resetButton);

        JButton helpButton = createSecondaryButton("Help");
        helpButton.addActionListener(e -> showHelp());
        buttonsPanel.add(helpButton);

        footer.add(buttonsPanel, BorderLayout.CENTER);

        JPanel statusBar = new JPanel(new BorderLayout());
        statusBar.setBorder(BorderFactory.createEmptyBorder(5, 15, 5, 15));
        statusLabel = new JLabel("Ready — configure parameters in the Workflow tab");
        statusLabel.setFont(derivedFont(Font.PLAIN, 11f));
        statusBar.add(statusLabel, BorderLayout.WEST);

        JLabel javaLabel = new JLabel("Java: " + System.getProperty("java.version"));
        javaLabel.setFont(derivedFont(Font.PLAIN, 11f));
        statusBar.add(javaLabel, BorderLayout.EAST);

        footer.add(statusBar, BorderLayout.SOUTH);

        return footer;
    }

    // ==================== COMMAND BUILDING ====================

    private List<String> buildCommandArgs() {
        List<String> cmd = new ArrayList<>();

        // Resolve the launcher we should re-spawn for the analysis subprocess:
        //
        //   - Dev / `java -jar fdrbench.jar` runs:  ProcessHandle returns the
        //     real java.exe -> use "<java> -jar <jar> <args>".
        //   - jpackage `FDRBench.exe` runs:         ProcessHandle returns the
        //     bundled launcher -> just re-spawn it with the args; the runtime
        //     does NOT ship java.exe (only javaw.exe), so going through the
        //     launcher avoids the "java.exe: CreateProcess error=2" failure.
        //     FDRBenchLauncher.main() routes args -> CLI for us.
        String javaExec = getJavaExecutable();
        boolean exeLaunch = false;
        if (javaExec.endsWith("java.exe") || javaExec.endsWith("java")) {
            // use as-is
        } else if (javaExec.toLowerCase().endsWith("fdrbench.exe")) {
            exeLaunch = true;
        } else {
            // unknown launcher -> fall back to PATH java
            javaExec = "java";
        }

        cmd.add(javaExec);
        if (!exeLaunch) {
            cmd.add("-jar");
            cmd.add(resolveJarPath());
        }

        if (currentWorkflow == WORKFLOW_DB_GENERATION) {
            addOption(cmd, "-db", dbFileField.getText());
            // Only forward -ms when the user explicitly chose Foreign Species,
            // so a stale path left in the field doesn't sneak into a Random
            // Shuffling run.
            if (seqMethodCombo != null
                    && seqMethodCombo.getSelectedIndex() == SEQ_METHOD_FOREIGN) {
                addOption(cmd, "-ms", foreignSpeciesArgValue());
            }
            File dbOut = deriveDbOutputFile();
            if (dbOut != null) {
                addOption(cmd, "-o", dbOut.getAbsolutePath());
            }
            addOption(cmd, "-level", String.valueOf(dbLevelCombo.getSelectedItem()));
            addOption(cmd, "-enzyme", String.valueOf(enzymeCombo.getSelectedIndex()));
            addOption(cmd, "-miss_c", String.valueOf(missCleavageSpinner.getValue()));
            addOption(cmd, "-minLength", String.valueOf(minLengthSpinner.getValue()));
            addOption(cmd, "-maxLength", String.valueOf(maxLengthSpinner.getValue()));
            addOption(cmd, "-fold", String.valueOf(foldSpinner.getValue()));
            addFlag(cmd, "-clip_n_m", clipNmCheckbox.isSelected());
            addOption(cmd, "-fix_nc", getFixNcValue());
            // -decoy only takes effect for Random Shuffling generation —
            // skip it for Foreign Species so a stale checked state doesn't
            // sneak into the run.
            boolean addDecoysActive = isAddDecoysApplicable()
                    && decoyCheckbox.isSelected();
            addFlag(cmd, "-decoy", addDecoysActive);
            addFlag(cmd, "-I2L", i2lCheckbox.isSelected());
            addFlag(cmd, "-diann", diannCheckbox.isSelected());
            addFlag(cmd, "-uniprot", uniprotCheckbox.isSelected());
            // -export_db only makes sense at peptide level; at protein level
            // the main -o output is already a protein FASTA, so skip the
            // flag (and ignore any stale checkbox state from a previous
            // peptide-level session).
            boolean isPeptideLevel = "peptide".equals(dbLevelCombo.getSelectedItem());
            addFlag(cmd, "-export_db", isPeptideLevel && exportDbCheckbox.isSelected());
            // -check is a no-op outside protein-level random shuffling
            // (only generate_protein_entrapment_database honors it). Skip it
            // elsewhere so a stale checked state doesn't show up in the
            // Preview Command output.
            boolean isShufflingMethod = seqMethodCombo == null
                    || seqMethodCombo.getSelectedIndex() == SEQ_METHOD_SHUFFLING;
            addFlag(cmd, "-check", !isPeptideLevel && isShufflingMethod
                    && checkDuplicatesCheckbox.isSelected());
            // -ns is a no-op outside protein-level foreign species (only
            // generate_protein_entrapment_database_from_multiple_species_data
            // honors it). Skip it elsewhere so a stale checked state doesn't
            // sneak into a different run.
            addFlag(cmd, "-ns", !isPeptideLevel && !isShufflingMethod
                    && noSharedCheckbox.isSelected());
        } else {
            addOption(cmd, "-i", inputFileField.getText());
            // -pep is peptide-level pair-tracking; skip it for protein-level
            // FDP so a stale path can't sneak into the run.
            if (!"protein".equals(fdpLevelCombo.getSelectedItem())) {
                addOption(cmd, "-pep", pepPairFileField.getText());
            }
            File fdpOut = deriveFdpOutputFile();
            if (fdpOut != null) {
                addOption(cmd, "-o", fdpOut.getAbsolutePath());
            }
            addOption(cmd, "-level", String.valueOf(fdpLevelCombo.getSelectedItem()));
            String scoreCol = getSelectedScoreColumn();
            if (scoreCol != null && !scoreCol.isEmpty()) {
                addOption(cmd, "-score", scoreCol + ":" + (scoreDirectionCombo.getSelectedIndex() == 0 ? "0" : "1"));
            }
            // -fold drives the k-fold path; -r drives the combined-entrapment
            // path. Foreign Species always uses -r. Random Shuffling normally
            // uses -fold, but the k-fold path requires a peptide pair file at
            // peptide-tier levels (peptide/precursor/psm); when no pair file
            // is provided, reroute the Fold value to -r so the CLI takes the
            // combined-entrapment path (which derives entrapment status from
            // the protein column instead).
            boolean isShuffling = seqMethodCombo == null
                    || seqMethodCombo.getSelectedIndex() == SEQ_METHOD_SHUFFLING;
            String fdpLevel = String.valueOf(fdpLevelCombo.getSelectedItem());
            boolean isPeptideTierLevel = "peptide".equals(fdpLevel)
                    || "precursor".equals(fdpLevel)
                    || "psm".equals(fdpLevel);
            boolean rerouteToR = isShuffling
                    && isPeptideTierLevel
                    && pepPairFileField.getText().trim().isEmpty();
            if (isShuffling && !rerouteToR) {
                addOption(cmd, "-fold", String.valueOf(fdpFoldSpinner.getValue()));
            }
            addOption(cmd, "-pick", String.valueOf(pickMethodCombo.getSelectedItem()));
            if (seqMethodCombo != null
                    && seqMethodCombo.getSelectedIndex() == SEQ_METHOD_FOREIGN) {
                addOption(cmd, "-r", rRatioField.getText());
            } else if (rerouteToR) {
                addOption(cmd, "-r", String.valueOf(fdpFoldSpinner.getValue()));
            }
        }

        addOption(cmd, "-seed", String.valueOf(seedSpinner.getValue()));
        // The entrapment label/pos for the FDP workflow lives in the FDP
        // settings panel so it can be set to match an externally-generated
        // database; the DB-generation panel keeps its own copy.
        boolean isFdpWorkflow = (currentWorkflow == WORKFLOW_FDP_ESTIMATION);
        String entrapmentLabel = isFdpWorkflow
                ? fdpEntrapmentLabelField.getText()
                : entrapmentLabelField.getText();
        int entrapmentPosIndex = isFdpWorkflow
                ? fdpEntrapmentPosCombo.getSelectedIndex()
                : entrapmentPosCombo.getSelectedIndex();
        addOptionIfChanged(cmd, "-entrapment_label", entrapmentLabel, "_p_target");
        addOptionIfChanged(cmd, "-decoy_label", decoyLabelField.getText(), "rev_");
        if (entrapmentPosIndex == 1) {
            addOption(cmd, "-entrapment_pos", "0");
        }
        if (decoyPosCombo.getSelectedIndex() == 1) {
            addOption(cmd, "-decoy_pos", "1");
        }
        addFlag(cmd, "-debug", debugCheckbox.isSelected());

        return cmd;
    }

    private String buildCommandPreview() {
        StringBuilder preview = new StringBuilder();
        for (String part : buildCommandArgs()) {
            if (preview.length() > 0) {
                preview.append(' ');
            }
            preview.append(quoteForPreview(part));
        }
        return preview.toString();
    }

    private void previewCommand() {
        String validationError = validateInputs();
        if (validationError != null) {
            JOptionPane.showMessageDialog(this, validationError, "Invalid Settings", JOptionPane.ERROR_MESSAGE);
            return;
        }

        String command = buildCommandPreview();
        JTextArea cmdArea = new JTextArea(command);
        cmdArea.setFont(new Font(Font.MONOSPACED, Font.PLAIN, 12));
        cmdArea.setEditable(false);
        cmdArea.setLineWrap(true);
        cmdArea.setWrapStyleWord(true);

        JScrollPane scrollPane = new JScrollPane(cmdArea);
        scrollPane.setPreferredSize(new Dimension(600, 200));

        int result = JOptionPane.showOptionDialog(this, scrollPane, "Command Preview",
                JOptionPane.DEFAULT_OPTION, JOptionPane.INFORMATION_MESSAGE, null,
                new String[] { "Copy to Clipboard", "Close" }, "Close");

        if (result == 0) {
            java.awt.datatransfer.StringSelection sel = new java.awt.datatransfer.StringSelection(command);
            Toolkit.getDefaultToolkit().getSystemClipboard().setContents(sel, null);
        }
    }

    private void showHelp() {
        String helpHtml = "<html>"
                + "<body style=\"font-family: 'Segoe UI', sans-serif; font-size: 12pt; padding: 10px;\">"
                + "<h2 style=\"margin-top: 0;\">FDRBench &mdash; FDR Control Evaluation Tool</h2>"
                + "<p>FDRBench is a tool for false discovery rate (FDR) control evaluation in "
                + "proteomics. It supports two workflows:</p>"
                + "<ol style=\"margin-left: 20px; padding-left: 0;\">"
                + "<li><b>Entrapment Database Generation</b> &mdash; build entrapment databases "
                + "either by randomly shuffling target sequences or by using sequences from "
                + "foreign species.</li>"
                + "<li><b>FDP Estimation / FDR Control Evaluation</b> &mdash; estimate the "
                + "false discovery proportion (FDP) using the lower-bound, combined, and "
                + "paired methods.</li>"
                + "</ol>"
                + "<h3>Quick Start:</h3>"
                + "<p><b>Entrapment database:</b></p>"
                + "<ul style=\"margin-left: 20px; padding-left: 0;\">"
                + "<li>Provide a target protein database (FASTA).</li>"
                + "<li>Pick the generation level (peptide or protein).</li>"
                + "<li>Choose Sequence Generation: Random Shuffling or Foreign Species "
                + "(provide foreign-species FASTA file(s) for the latter; the Download "
                + "button can fetch a proteome from UniProt).</li>"
                + "<li>Set the output file and click Run FDRBench.</li>"
                + "</ul>"
                + "<p><b>FDP estimation:</b></p>"
                + "<ul style=\"margin-left: 20px; padding-left: 0;\">"
                + "<li>Provide a PSM/peptide/precursor/protein file.</li>"
                + "<li>Pick the calculation level and the score column.</li>"
                + "<li>Configure fold or R ratio as needed and click Run FDRBench.</li>"
                + "</ul>"
                + "<p>For detailed documentation and examples, visit:<br/>"
                + "<a href=\"https://github.com/Noble-Lab/FDRBench\">"
                + "https://github.com/Noble-Lab/FDRBench</a></p>"
                + "<h3>How to cite:</h3>"
                + "<p>Wen, B., Freestone, J., Riffle, M. <i>et al.</i> "
                + "<a href=\"https://doi.org/10.1038/s41592-025-02719-x\">Assessment of false "
                + "discovery rate control in tandem mass spectrometry analysis using entrapment</a>.<br/>"
                + "<i>Nat Methods</i> <b>22</b>, 1454&ndash;1463 (2025).</p>"
                + "</body></html>";

        JEditorPane editorPane = new JEditorPane("text/html", helpHtml);
        editorPane.setEditable(false);
        editorPane.setOpaque(false);
        editorPane.putClientProperty(JEditorPane.HONOR_DISPLAY_PROPERTIES, Boolean.TRUE);

        editorPane.addHyperlinkListener(e -> {
            if (e.getEventType() == javax.swing.event.HyperlinkEvent.EventType.ACTIVATED) {
                try {
                    Desktop.getDesktop().browse(e.getURL().toURI());
                } catch (Exception ex) {
                    // Silently ignore — the dialog itself isn't blocked by a failed browse.
                }
            }
        });

        JScrollPane scrollPane = new JScrollPane(editorPane);
        scrollPane.setPreferredSize(new Dimension(560, 460));
        scrollPane.setBorder(BorderFactory.createEmptyBorder());
        scrollPane.setViewportBorder(BorderFactory.createEmptyBorder());

        JOptionPane.showMessageDialog(this, scrollPane, "FDRBench Help",
                JOptionPane.INFORMATION_MESSAGE);
    }

    // ==================== EXECUTION ====================

    private void runFDRBench() {
        if (isRunning)
            return;

        String validationError = validateInputs();
        if (validationError != null) {
            JOptionPane.showMessageDialog(this, validationError, "Invalid Settings", JOptionPane.ERROR_MESSAGE);
            return;
        }

        List<String> command = buildCommandArgs();
        // Open the run-log file before any logToConsole — that way the
        // header lines below are mirrored into fdrbench_log.txt too.
        openLogWriter();
        String startStamp = new java.text.SimpleDateFormat("yyyy-MM-dd HH:mm:ss")
                .format(new java.util.Date());
        logToConsole("Starting FDRBench...\n");
        logToConsole("Date: " + startStamp + "\n");
        logToConsole("Workflow: " + (currentWorkflow == WORKFLOW_DB_GENERATION
                ? "Entrapment Database Generation"
                : "FDP Estimation") + "\n");
        logToConsole("Command: " + buildCommandPreview() + "\n\n");

        isRunning = true;
        runButton.setEnabled(false);
        stopButton.setEnabled(true);
        progressBar.setIndeterminate(true);
        progressBar.setString("Running…");
        statusLabel.setText("Running FDRBench…");

        // Switch to the Console tab — find it by name since the Plot tab may
        // shift the Console index.
        int consoleIdx = tabbedPane.indexOfTab("Console");
        if (consoleIdx >= 0) {
            tabbedPane.setSelectedIndex(consoleIdx);
        }

        Thread workerThread = new Thread(() -> {
            try {
                ProcessBuilder pb = new ProcessBuilder(command);
                pb.redirectErrorStream(true);
                pb.directory(new File(System.getProperty("user.dir")));
                currentProcess = pb.start();

                try (BufferedReader reader = new BufferedReader(
                        new InputStreamReader(currentProcess.getInputStream()))) {
                    String line;
                    while ((line = reader.readLine()) != null) {
                        final String l = line;
                        SwingUtilities.invokeLater(() -> logToConsole(l + "\n"));
                    }
                }

                int exitCode = currentProcess.waitFor();
                SwingUtilities.invokeLater(() -> {
                    logToConsole("\nProcess finished with exit code: " + exitCode + "\n");
                    finishExecution(exitCode == 0);
                });

            } catch (Exception e) {
                SwingUtilities.invokeLater(() -> {
                    logToConsole("\nError: " + e.getMessage() + "\n");
                    finishExecution(false);
                });
            }
        });
        workerThread.setName("fdrbench-runner");
        workerThread.setDaemon(true);
        workerThread.start();
    }

    private void stopFDRBench() {
        if (!isRunning) {
            return;
        }
        Process proc = currentProcess;
        if (proc != null && proc.isAlive()) {
            proc.destroyForcibly();
            logToConsole("\n[Stopped by user]\n");
        }
        // The worker thread will call finishExecution() once the stream
        // closes and waitFor() returns, so we don't need to call it here.
    }

    private void finishExecution(boolean success) {
        if (!isRunning) {
            // Already finalized — avoid double-firing the auto-reset timer.
            return;
        }
        isRunning = false;
        currentProcess = null;
        runButton.setEnabled(true);
        stopButton.setEnabled(false);
        progressBar.setIndeterminate(false);
        String msg = success ? "Completed successfully" : "Stopped";
        progressBar.setString(msg);
        statusLabel.setText(msg);

        // Stamp the log with the result then flush + close the file.
        String endStamp = new java.text.SimpleDateFormat("yyyy-MM-dd HH:mm:ss")
                .format(new java.util.Date());
        logToConsole("\n[" + endStamp + "] " + msg + "\n");
        closeLogWriter();

        if (success) {
            // Both workflows now configure an output folder; the actual
            // output file is auto-derived per workflow (deriveDbOutputFile /
            // deriveFdpOutputFile).
            File outFile = (currentWorkflow == WORKFLOW_DB_GENERATION)
                    ? deriveDbOutputFile()
                    : deriveFdpOutputFile();
            if (outFile != null && outFile.exists()) {
                // Load the FDP plot and switch focus to it before the
                // success dialog opens, so the user lands on the chart.
                if (currentWorkflow == WORKFLOW_FDP_ESTIMATION
                        && fdpPlotPanel != null) {
                    fdpPlotPanel.loadFromCsv(outFile);
                    int plotIdx = tabbedPane.indexOfTab(PLOT_TAB_TITLE);
                    if (plotIdx >= 0) {
                        tabbedPane.setSelectedIndex(plotIdx);
                    }
                }
                showSuccessDialog(outFile);
            }
        }

        // Auto-reset UI to idle after 4 seconds
        javax.swing.Timer resetTimer = new javax.swing.Timer(4000, e -> {
            if (!isRunning) {
                progressBar.setString("");
                statusLabel.setText("Ready — configure parameters in the Workflow tab");
            }
        });
        resetTimer.setRepeats(false);
        resetTimer.start();
    }

    private void showSuccessDialog(File outFile) {
        java.util.List<File> generated = findGeneratedFiles(outFile);
        if (generated.isEmpty() && outFile != null) {
            generated.add(outFile);
        }

        JDialog dialog = new JDialog(this, "Success", true);
        JPanel content = new JPanel(new BorderLayout(0, 10));
        content.setBorder(BorderFactory.createEmptyBorder(15, 15, 10, 15));
        content.add(new JLabel("Process completed successfully!"), BorderLayout.NORTH);

        // If more than one related file was produced, expose them in a
        // dropdown so the user can View / Open Folder against any of them.
        // Otherwise fall back to a non-editable path text field.
        JComboBox<File> combo = null;
        JTextField pathField = null;
        if (generated.size() > 1) {
            combo = new JComboBox<>(generated.toArray(new File[0]));
            combo.setRenderer(new DefaultListCellRenderer() {
                @Override
                public Component getListCellRendererComponent(JList<?> list,
                        Object value, int index, boolean isSelected,
                        boolean cellHasFocus) {
                    super.getListCellRendererComponent(list, value, index,
                            isSelected, cellHasFocus);
                    if (value instanceof File) {
                        setText(((File) value).getName());
                    }
                    return this;
                }
            });
            content.add(combo, BorderLayout.CENTER);
        } else {
            pathField = new JTextField(outFile != null
                    ? outFile.getAbsolutePath() : "");
            pathField.setEditable(false);
            content.add(pathField, BorderLayout.CENTER);
        }

        JPanel buttons = new JPanel(new FlowLayout(FlowLayout.RIGHT, 8, 0));
        JButton viewButton = new JButton("View");
        styleButton(viewButton);
        JButton openFolderButton = new JButton("Open Folder");
        styleButton(openFolderButton);
        JButton closeButton = new JButton("Close");
        styleButton(closeButton);
        buttons.add(viewButton);
        buttons.add(openFolderButton);
        buttons.add(closeButton);
        content.add(buttons, BorderLayout.SOUTH);

        final JComboBox<File> finalCombo = combo;
        final File defaultFile = outFile;
        java.util.function.Supplier<File> selected = () ->
                finalCombo != null
                        ? (File) finalCombo.getSelectedItem()
                        : defaultFile;

        viewButton.addActionListener(e -> {
            File f = selected.get();
            if (f != null && f.exists()) {
                viewFile(f);
            }
        });
        openFolderButton.addActionListener(e -> {
            File f = selected.get();
            File folder = f != null ? f.getParentFile() : null;
            if (folder != null && folder.exists()) {
                try {
                    Desktop.getDesktop().open(folder);
                } catch (Exception ex) {
                    JOptionPane.showMessageDialog(dialog,
                            "Could not open folder: " + ex.getMessage());
                }
            }
        });
        closeButton.addActionListener(e -> dialog.dispose());

        dialog.setContentPane(content);
        dialog.getRootPane().setDefaultButton(closeButton);
        dialog.pack();
        dialog.setLocationRelativeTo(this);
        dialog.setVisible(true);
    }

    /**
     * Enumerate files in the output folder that were produced by the run.
     * Strategy: take the primary output's name, strip its extension and any
     * level suffix (<code>_pep</code>, <code>_pro</code>, <code>_protein</code>),
     * and grab everything else in the same folder that starts with that
     * prefix — typically the FASTA companion, the exported protein DB, etc.
     * Falls back to just the primary output if nothing else matches.
     */
    private static java.util.List<File> findGeneratedFiles(File primaryOutput) {
        java.util.List<File> result = new ArrayList<>();
        if (primaryOutput == null) {
            return result;
        }
        if (primaryOutput.exists()) {
            result.add(primaryOutput);
        }
        File folder = primaryOutput.getParentFile();
        if (folder == null || !folder.isDirectory()) {
            return result;
        }
        String name = primaryOutput.getName();
        int dot = name.lastIndexOf('.');
        String base = dot > 0 ? name.substring(0, dot) : name;
        String prefix = base;
        if (base.endsWith("_pep") || base.endsWith("_pro")) {
            prefix = base.substring(0, base.length() - 4);
        } else if (base.endsWith("_protein")) {
            prefix = base.substring(0, base.length() - "_protein".length());
        }
        String prefixLower = prefix.toLowerCase(Locale.ROOT);
        File[] files = folder.listFiles();
        if (files == null) {
            return result;
        }
        for (File f : files) {
            if (f.isFile() && !f.equals(primaryOutput)
                    && f.getName().toLowerCase(Locale.ROOT).startsWith(prefixLower)) {
                result.add(f);
            }
        }
        // Primary first, then alphabetical.
        result.sort((a, b) -> {
            if (a.equals(primaryOutput)) return -1;
            if (b.equals(primaryOutput)) return 1;
            return a.getName().compareToIgnoreCase(b.getName());
        });
        return result;
    }

    /**
     * View dispatch — tabular preview for {@code .csv}/{@code .tsv}/{@code .txt}
     * files; plain-text preview for FASTA so the entries don't get mangled
     * into a one-column table.
     */
    private void viewFile(File file) {
        if (file == null || !file.exists()) {
            return;
        }
        String n = file.getName().toLowerCase(Locale.ROOT);
        if (n.endsWith(".fasta") || n.endsWith(".fa") || n.endsWith(".faa")) {
            showTextPreview(file);
        } else {
            showFilePreview(file.getAbsolutePath());
        }
    }

    /**
     * FASTA-record-aware preview: keeps the first {@value #PREVIEW_FASTA_ENTRIES}
     * records and the last {@value #PREVIEW_FASTA_ENTRIES} (sliding window) so
     * the user can spot-check both ends of large protein databases without
     * loading the whole file. A "…" separator marks the gap when the file
     * has more than 2 × PREVIEW_FASTA_ENTRIES records.
     */
    private static final int PREVIEW_FASTA_ENTRIES = 100;

    private void showTextPreview(File file) {
        java.util.List<java.util.List<String>> topRecs = new ArrayList<>();
        java.util.Deque<java.util.List<String>> tailRecs = new java.util.ArrayDeque<>();
        long totalRecs = 0;
        try (BufferedReader r = new BufferedReader(
                new InputStreamReader(new java.io.FileInputStream(file),
                        java.nio.charset.StandardCharsets.UTF_8))) {
            java.util.List<String> current = null;
            String line;
            while ((line = r.readLine()) != null) {
                if (line.startsWith(">")) {
                    if (current != null) {
                        commitFastaRecord(current, topRecs, tailRecs);
                        totalRecs++;
                    }
                    current = new ArrayList<>();
                    current.add(line);
                } else if (current != null) {
                    current.add(line);
                }
                // Lines before the first ">" header (rare; comments etc.) are
                // dropped — they aren't part of any record.
            }
            if (current != null) {
                commitFastaRecord(current, topRecs, tailRecs);
                totalRecs++;
            }
        } catch (IOException ex) {
            JOptionPane.showMessageDialog(this,
                    "Could not read file: " + ex.getMessage(),
                    "Read error", JOptionPane.ERROR_MESSAGE);
            return;
        }

        StringBuilder sb = new StringBuilder();
        for (java.util.List<String> rec : topRecs) {
            for (String l : rec) sb.append(l).append('\n');
        }
        boolean truncated = !tailRecs.isEmpty();
        if (truncated) {
            sb.append("…\n");
            for (java.util.List<String> rec : tailRecs) {
                for (String l : rec) sb.append(l).append('\n');
            }
        }

        JTextArea area = new JTextArea(sb.toString());
        area.setFont(new Font(Font.MONOSPACED, Font.PLAIN, 12));
        area.setEditable(false);
        area.setLineWrap(false);
        area.setCaretPosition(0);
        JScrollPane sp = new JScrollPane(area);
        sp.setPreferredSize(new Dimension(800, 500));
        sp.setBorder(BorderFactory.createEmptyBorder());

        String shown = truncated
                ? "showing first " + topRecs.size()
                        + " + last " + tailRecs.size()
                : "showing all " + topRecs.size();
        String title = "Preview — " + file.getName()
                + " (" + shown + " of " + totalRecs + " entries)";
        JOptionPane.showMessageDialog(this, sp, title, JOptionPane.PLAIN_MESSAGE);
    }

    private static void commitFastaRecord(java.util.List<String> rec,
            java.util.List<java.util.List<String>> topRecs,
            java.util.Deque<java.util.List<String>> tailRecs) {
        if (topRecs.size() < PREVIEW_FASTA_ENTRIES) {
            topRecs.add(rec);
        } else {
            if (tailRecs.size() == PREVIEW_FASTA_ENTRIES) {
                tailRecs.pollFirst();
            }
            tailRecs.addLast(rec);
        }
    }

    /**
     * Writes a message to the on-screen Console tab and, when a run is in
     * progress, mirrors it into {@code fdrbench_log.txt} inside the chosen
     * output folder.
     */
    private synchronized void logToConsole(String message) {
        consoleArea.append(message);
        consoleArea.setCaretPosition(consoleArea.getDocument().getLength());
        if (logWriter != null) {
            try {
                logWriter.write(message);
                logWriter.flush();
            } catch (IOException e) {
                // Don't recurse — print to stderr if the log itself is broken.
                System.err.println("Failed to write to log: " + e.getMessage());
            }
        }
    }

    /** Open fdrbench_log.txt in the active workflow's output folder. */
    private void openLogWriter() {
        closeLogWriter();
        String folder = currentWorkflow == WORKFLOW_DB_GENERATION
                ? (dbOutputField != null ? dbOutputField.getText().trim() : "")
                : (fdpOutputField != null ? fdpOutputField.getText().trim() : "");
        if (folder.isEmpty()) {
            return;
        }
        File dir = new File(folder);
        if (!dir.exists()) {
            dir.mkdirs();
        }
        try {
            logWriter = new BufferedWriter(new FileWriter(
                    new File(dir, "fdrbench_log.txt")));
        } catch (IOException e) {
            logToConsole("Failed to create log file: " + e.getMessage() + "\n");
        }
    }

    private void closeLogWriter() {
        if (logWriter != null) {
            try {
                logWriter.close();
            } catch (IOException e) {
                System.err.println("Failed to close log: " + e.getMessage());
            }
            logWriter = null;
        }
    }

    private void resetToDefaults() {
        if (isRunning) return;

        int choice = JOptionPane.showConfirmDialog(this,
                "Reset all parameters to their default values?",
                "Reset Defaults", JOptionPane.YES_NO_OPTION);
        if (choice != JOptionPane.YES_OPTION) return;

        // Workflow 1: DB Generation
        dbFileField.setText("");
        foreignSpeciesField.setText("");
        foreignSpeciesFiles.clear();
        updateFileFieldState(foreignSpeciesField, foreignSpeciesFiles);
        dbOutputField.setText("");
        dbLevelCombo.setSelectedIndex(0);
        enzymeCombo.setSelectedIndex(defaultEnzymeIndex());
        missCleavageSpinner.setValue(1);
        minLengthSpinner.setValue(7);
        maxLengthSpinner.setValue(35);
        foldSpinner.setValue(1);
        clipNmCheckbox.setSelected(false);
        fixNcCombo.setSelectedIndex(0);
        i2lCheckbox.setSelected(true);
        decoyCheckbox.setSelected(false);
        diannCheckbox.setSelected(true);
        uniprotCheckbox.setSelected(true);
        exportDbCheckbox.setSelected(false);
        checkDuplicatesCheckbox.setSelected(true);
        noSharedCheckbox.setSelected(false);
        if (seqMethodCombo != null) {
            seqMethodCombo.setSelectedIndex(SEQ_METHOD_SHUFFLING);
        }
        updateForeignSpeciesVisibility();

        // Workflow 2: FDP Estimation
        inputFileField.setText("");
        pepPairFileField.setText("");
        fdpOutputField.setText("");
        fdpLevelCombo.setSelectedIndex(0);
        scoreColumnCombo.setModel(new javax.swing.DefaultComboBoxModel<>(
                new String[] { NO_SCORE_ITEM }));
        scoreColumnCombo.setSelectedItem(NO_SCORE_ITEM);
        scoreDirectionCombo.setSelectedIndex(0);
        fdpFoldSpinner.setValue(1);
        pickMethodCombo.setSelectedIndex(0);
        rRatioField.setText("");
        fdpEntrapmentLabelField.setText("_p_target");
        fdpEntrapmentPosCombo.setSelectedIndex(0);

        // Advanced
        seedSpinner.setValue(2000);
        entrapmentLabelField.setText("_p_target");
        decoyLabelField.setText("rev_");
        entrapmentPosCombo.setSelectedIndex(0);
        decoyPosCombo.setSelectedIndex(0);
        debugCheckbox.setSelected(false);
        updateDecoyLabelVisibility();

        logToConsole("All parameters reset to defaults.\n");
    }

    private String validateInputs() {
        if (currentWorkflow == WORKFLOW_DB_GENERATION) {
            if ((Integer) minLengthSpinner.getValue() > (Integer) maxLengthSpinner.getValue()) {
                return "Minimum peptide length cannot be greater than maximum peptide length.";
            }

            String dbPath = dbFileField.getText().trim();
            if (dbPath.isEmpty()) {
                return "Protein database is required for entrapment database generation.";
            }
            if (!new File(dbPath).exists()) {
                return "The specified Protein Database file does not exist:\n" + dbPath;
            }

            String outPath = dbOutputField.getText().trim();
            if (outPath.isEmpty()) {
                return "Output folder is required for entrapment database generation.";
            }
            File outFolder = new File(outPath);
            if (!outFolder.exists()) {
                return "The output folder does not exist:\n" + outPath;
            }
            if (!outFolder.isDirectory()) {
                return "The output path is not a folder:\n" + outPath;
            }

            if (seqMethodCombo != null
                    && seqMethodCombo.getSelectedIndex() == SEQ_METHOD_FOREIGN) {
                String foreignPaths = foreignSpeciesArgValue().trim();
                if (foreignPaths.isEmpty()) {
                    return "Foreign species FASTA file(s) are required when "
                            + "Sequence Generation is set to \"Foreign Species\".";
                }
                for (String path : foreignPaths.split(",")) {
                    String p = path.trim();
                    if (!p.isEmpty() && !new File(p).exists()) {
                        return "Foreign species file does not exist:\n" + p;
                    }
                }
            }

        } else if (currentWorkflow == WORKFLOW_FDP_ESTIMATION) {
            String inPath = inputFileField.getText().trim();
            if (inPath.isEmpty()) {
                return "Input file is required for FDP estimation.";
            }
            if (!new File(inPath).exists()) {
                return "The specified Input File does not exist:\n" + inPath;
            }

            String outPath = fdpOutputField.getText().trim();
            if (outPath.isEmpty()) {
                return "Output folder is required for FDP estimation.";
            }
            File outFolder = new File(outPath);
            if (!outFolder.exists()) {
                return "The output folder does not exist:\n" + outPath;
            }
            if (!outFolder.isDirectory()) {
                return "The output path is not a folder:\n" + outPath;
            }

            // Peptide Pair File is optional — only validate the path when the
            // user actually typed something AND the row is currently
            // applicable (it's hidden at protein level).
            if (!"protein".equals(fdpLevelCombo.getSelectedItem())) {
                String pepPath = pepPairFileField.getText().trim();
                if (!pepPath.isEmpty() && !new File(pepPath).exists()) {
                    return "The specified Peptide Pair File does not exist:\n"
                            + pepPath;
                }
            }

            // R Ratio is only meaningful for Foreign Species combined
            // entrapment, so skip the numeric check when Random Shuffling is
            // selected (the field is hidden in that case).
            if (seqMethodCombo != null
                    && seqMethodCombo.getSelectedIndex() == SEQ_METHOD_FOREIGN) {
                String rRatio = rRatioField.getText().trim();
                if (!rRatio.isEmpty()) {
                    try {
                        Double.parseDouble(rRatio);
                    } catch (NumberFormatException e) {
                        return "R ratio must be numeric when provided.";
                    }
                }
            }
        }

        return null;
    }

    private void addOption(List<String> cmd, String option, String value) {
        if (value != null) {
            String trimmed = value.trim();
            if (!trimmed.isEmpty()) {
                cmd.add(option);
                cmd.add(trimmed);
            }
        }
    }

    /**
     * Returns the user's chosen Score column, or {@code null} when the user
     * picked {@link #NO_SCORE_ITEM} so FDRBench should rank on q_value alone.
     */
    private String getSelectedScoreColumn() {
        if (scoreColumnCombo == null) return null;
        Object sel = scoreColumnCombo.getSelectedItem();
        if (sel == null) return null;
        String s = sel.toString().trim();
        if (s.isEmpty() || NO_SCORE_ITEM.equals(s)) return null;
        return s;
    }

    /**
     * Refresh the Score Column dropdown from the given input file. Reads the
     * header line plus the first data row, keeps columns whose first-row value
     * parses as a number, drops {@code q_value} (which FDRBench already uses
     * for primary ranking), and surfaces {@code score} at the top when present.
     * If {@code file} is null/unreadable, falls back to the static defaults.
     */
    private void refreshScoreColumns(File file) {
        if (scoreColumnCombo == null) return;
        String previous = scoreColumnCombo.getSelectedItem() == null
                ? null
                : scoreColumnCombo.getSelectedItem().toString();

        java.util.List<String> items = new java.util.ArrayList<>();
        items.add(NO_SCORE_ITEM);

        if (file != null && file.isFile()) {
            try (BufferedReader r = new BufferedReader(new InputStreamReader(
                    new java.io.FileInputStream(file), java.nio.charset.StandardCharsets.UTF_8))) {
                String headerLine = r.readLine();
                String firstData = headerLine == null ? null : r.readLine();
                if (headerLine != null && firstData != null) {
                    char delim = headerLine.indexOf('\t') >= 0 ? '\t'
                            : (headerLine.indexOf(',') >= 0 ? ',' : '\t');
                    String[] headers = headerLine.split(java.util.regex.Pattern.quote(String.valueOf(delim)), -1);
                    String[] values  = firstData.split(java.util.regex.Pattern.quote(String.valueOf(delim)), -1);
                    java.util.List<String> numeric = new java.util.ArrayList<>();
                    for (int i = 0; i < headers.length && i < values.length; i++) {
                        String name = headers[i].trim();
                        if (name.isEmpty()
                                || "q_value".equalsIgnoreCase(name)
                                || "charge".equalsIgnoreCase(name)) {
                            continue;
                        }
                        String v = values[i].trim();
                        if (v.isEmpty()) continue;
                        try {
                            Double.parseDouble(v);
                            numeric.add(name);
                        } catch (NumberFormatException ignored) { }
                    }
                    // Promote "score" to the top of the actionable list when present.
                    int scoreIdx = -1;
                    for (int i = 0; i < numeric.size(); i++) {
                        if ("score".equalsIgnoreCase(numeric.get(i))) { scoreIdx = i; break; }
                    }
                    if (scoreIdx >= 0) {
                        items.add(numeric.remove(scoreIdx));
                    }
                    items.addAll(numeric);
                }
            } catch (IOException ignored) {
                // Fall through to defaults below.
            }
        }

        scoreColumnCombo.setModel(new javax.swing.DefaultComboBoxModel<>(
                items.toArray(new String[0])));
        // Preserve the user's selection if still valid; otherwise prefer
        // "score", else fall back to "none".
        if (previous != null && items.contains(previous)) {
            scoreColumnCombo.setSelectedItem(previous);
        } else if (items.contains("score")) {
            scoreColumnCombo.setSelectedItem("score");
        } else {
            scoreColumnCombo.setSelectedItem(NO_SCORE_ITEM);
        }
    }

    private void addOptionIfChanged(List<String> cmd, String option, String value, String defaultValue) {
        if (value != null) {
            String trimmed = value.trim();
            if (!trimmed.isEmpty() && !trimmed.equals(defaultValue)) {
                cmd.add(option);
                cmd.add(trimmed);
            }
        }
    }

    private void addFlag(List<String> cmd, String option, boolean enabled) {
        if (enabled) {
            cmd.add(option);
        }
    }

    private String getFixNcValue() {
        if (fixNcCombo.getSelectedIndex() == 0) {
            return "c";
        }
        if (fixNcCombo.getSelectedIndex() == 1) {
            return "n";
        }
        return "nc";
    }

    /**
     * Find the launcher used to start this JVM.
     * {@code getJavaExecutable()}.
     *
     * <p>{@link ProcessHandle#current()} (Java 9+) returns the actual binary
     * the OS spawned, which is:
     * <ul>
     *   <li>{@code .../bin/java(.exe)} for {@code java -jar fdrbench.jar} runs
     *       (developer / CLI), and</li>
     *   <li>{@code .../FDRBench.exe} for jpackage app-image runs.</li>
     * </ul>
     * Returning the launcher itself lets {@link #buildCommandArgs()} re-spawn
     * the bundled {@code .exe} for analyses — necessary because jpackage's
     * stripped runtime ships {@code javaw.exe} only (no {@code java.exe}).
     */
    @SuppressWarnings("unchecked")
    private String getJavaExecutable() {
        // Reflective ProcessHandle.current().info().command() — the class is
        // JDK 9+, but pom.xml targets 1.8, so we can't reference it directly.
        // The jpackage runtime is JDK 17+, so this branch always succeeds in
        // the packaged build.
        try {
            Class<?> phCls = Class.forName("java.lang.ProcessHandle");
            Class<?> phInfoCls = Class.forName("java.lang.ProcessHandle$Info");
            Object current = phCls.getMethod("current").invoke(null);
            Object info = phCls.getMethod("info").invoke(current);
            java.util.Optional<String> cmd =
                    (java.util.Optional<String>) phInfoCls.getMethod("command").invoke(info);
            if (cmd != null && cmd.isPresent()) {
                return cmd.get();
            }
        } catch (Throwable ignored) {
            // Pre-JDK-9 or restricted; fall through.
        }
        String javaHome = System.getProperty("java.home");
        String sep = File.separator;
        return javaHome + sep + "bin" + sep
                + (System.getProperty("os.name").toLowerCase().contains("win")
                        ? "java.exe" : "java");
    }

    private String resolveJarPath() {
        // First: ask the ClassLoader where this class actually came from.
        // When launched via "java -jar /abs/path/foo.jar", the protection
        // domain reports the JAR's absolute URL — which is what we must pass
        // to the subprocess so it works regardless of the user's CWD. The
        // CWD heuristics below only help for IDE runs (classes/, not a JAR)
        // and shouldn't shadow the truth from the ClassLoader.
        try {
            java.security.CodeSource src = FDRBenchGUI.class
                    .getProtectionDomain().getCodeSource();
            if (src != null && src.getLocation() != null) {
                File f = new File(src.getLocation().toURI());
                if (f.isFile() && f.getName().toLowerCase(Locale.ROOT).endsWith(".jar")) {
                    return f.getAbsolutePath();
                }
            }
        } catch (Exception ignored) {
            // Fall through to CWD-based discovery (e.g. IDE runs from target/classes).
        }

        // Try to find the JAR file in target/ or the same folder as the app
        File currentDir = new File(System.getProperty("user.dir"));
        List<File> jarCandidates = new ArrayList<>();

        // Check current directory
        collectJarCandidates(currentDir, jarCandidates, 1);

        // Check target/
        File targetDir = new File(currentDir, "target");
        if (targetDir.exists()) {
            collectJarCandidates(targetDir, jarCandidates, 2);
        }

        if (!jarCandidates.isEmpty()) {
            // Sort by last modified to get the newest one
            jarCandidates.sort(Comparator.comparingLong(File::lastModified).reversed());
            return jarCandidates.get(0).getAbsolutePath();
        }

        // Fallback to a default name
        return "fdrbench.jar";
    }

    private void collectJarCandidates(File dir, List<File> jarCandidates, int depth) {
        if (dir == null || !dir.exists() || depth < 0) {
            return;
        }
        File[] files = dir.listFiles();
        if (files == null) {
            return;
        }
        Arrays.sort(files, Comparator.comparing(File::getName));
        for (File file : files) {
            if (file.isDirectory()) {
                collectJarCandidates(file, jarCandidates, depth - 1);
            } else if (file.getName().startsWith("fdrbench-") && file.getName().endsWith(".jar")) {
                jarCandidates.add(file);
            }
        }
    }

    private String quoteForPreview(String value) {
        if (value.indexOf(' ') >= 0 || value.indexOf(',') >= 0) {
            return "\"" + value.replace("\"", "\\\"") + "\"";
        }
        return value;
    }

    // ==================== UI HELPERS ====================

    private JLabel createLabel(String text, String tooltip) {
        // No explicit font: inherits the 13pt defaultFont set in
        // customizeUIDefaults().
        JLabel label = new JLabel(text);
        if (tooltip != null)
            label.setToolTipText(tooltip);
        return label;
    }

    private JTextField createTextField(String placeholder) {
        JTextField field = new JTextField();
        // Pin all three size hints to COMPONENT_HEIGHT so layout managers
        // can't grow the field past it (e.g. BorderLayout.CENTER inside the
        // Plot tab) and FlatLaf's intrinsic preferred height can't shrink it
        // below it.
        field.setPreferredSize(new Dimension(200, COMPONENT_HEIGHT));
        field.setMinimumSize(new Dimension(0, COMPONENT_HEIGHT));
        field.setMaximumSize(new Dimension(Integer.MAX_VALUE, COMPONENT_HEIGHT));
        field.setToolTipText(placeholder);
        field.putClientProperty("JTextField.placeholderText", placeholder);
        // Keep the tooltip in sync with the field's content so long paths
        // remain inspectable even when truncated by the field's width.
        field.getDocument().addDocumentListener(new javax.swing.event.DocumentListener() {
            private void update() {
                String text = field.getText();
                field.setToolTipText(text == null || text.isEmpty() ? placeholder : text);
            }
            @Override public void insertUpdate(javax.swing.event.DocumentEvent e) { update(); }
            @Override public void removeUpdate(javax.swing.event.DocumentEvent e) { update(); }
            @Override public void changedUpdate(javax.swing.event.DocumentEvent e) { update(); }
        });
        return field;
    }

    private JSpinner createSpinner(int value, int min, int max, int step) {
        JSpinner spinner = new JSpinner(new SpinnerNumberModel(value, min, max, step));
        // Pin the editor to a fixed column count and our shared component
        // height so all input rows line up. Width is left to the LAF.
        ((JSpinner.DefaultEditor) spinner.getEditor()).getTextField().setColumns(5);
        Dimension prefSize = spinner.getPreferredSize();
        spinner.setPreferredSize(new Dimension(prefSize.width, COMPONENT_HEIGHT));
        spinner.setMinimumSize(new Dimension(60, COMPONENT_HEIGHT));
        return spinner;
    }

    private void styleComboBox(JComboBox<?> combo) {
        // Match the shared component height so combos line up with text
        // fields and spinners; preserve the LAF's natural preferred width so
        // the popup still fits the longest item.
        Dimension prefSize = combo.getPreferredSize();
        combo.setPreferredSize(new Dimension(prefSize.width, COMPONENT_HEIGHT));
    }

    private void styleButton(JButton button) {
        button.setFont(button.getFont().deriveFont(Font.PLAIN, 12f));
        button.setFocusPainted(false);
        button.setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
        button.setMargin(new Insets(6, 12, 6, 12));
        button.putClientProperty("JButton.buttonType", "roundRect");
        button.putClientProperty("JButton.hoverBackground",
                UIManager.getColor("Button.hoverBackground"));
    }

    /**
     * Browse + Download buttons for a protein-database row. The Download button opens the UniProt
     * download dialog and writes the resulting FASTA path back into the field.
     */
    private JPanel createDbButtonsPanel(JTextField targetField) {
        JPanel panel = new JPanel(new GridLayout(1, 0, 5, 0));
        panel.setOpaque(false);

        panel.add(createBrowseButton(targetField, "fasta"));

        JButton downloadButton = new JButton("Download");
        styleButton(downloadButton);
        downloadButton.setToolTipText("Download protein database from UniProt");
        downloadButton.addActionListener(e -> {
            // Use the configured Output Folder as the default download
            // directory; the field already holds a directory path now.
            String defaultDir = dbOutputField != null
                    ? dbOutputField.getText().trim()
                    : "";
            if (defaultDir.isEmpty()) {
                defaultDir = prefs.get(PREF_LAST_DIR, "");
            }
            UniProtDownloadDialog dialog = new UniProtDownloadDialog(this, targetField, defaultDir);
            dialog.showDialog();
        });
        panel.add(downloadButton);

        return panel;
    }

    /**
     * Folder + Open buttons for the entrapment-DB Output Folder row: the user picks a directory; "Open" reveals it in the
     * system file explorer (handy after a successful run).
     */
    private JPanel createFolderButton(JTextField targetField) {
        JPanel panel = new JPanel(new GridLayout(1, 0, 5, 0));
        panel.setOpaque(false);

        JButton folderButton = new JButton("Folder");
        styleButton(folderButton);
        folderButton.setToolTipText("Choose an output folder");
        folderButton.addActionListener(e -> {
            JFileChooser chooser = new JFileChooser();
            chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
            String lastDir = prefs.get(PREF_LAST_DIR, System.getProperty("user.home"));
            chooser.setCurrentDirectory(new File(lastDir));
            if (chooser.showOpenDialog(this) == JFileChooser.APPROVE_OPTION) {
                File f = chooser.getSelectedFile();
                targetField.setText(f.getAbsolutePath());
                prefs.put(PREF_LAST_DIR, f.getAbsolutePath());
            }
        });
        panel.add(folderButton);

        JButton openButton = new JButton("Open");
        styleButton(openButton);
        openButton.setToolTipText("Open the output folder in file explorer");
        openButton.addActionListener(e -> {
            String path = targetField.getText().trim();
            if (path.isEmpty()) {
                JOptionPane.showMessageDialog(this,
                        "Please choose an output folder first.",
                        "No folder", JOptionPane.WARNING_MESSAGE);
                return;
            }
            File folder = new File(path);
            if (!folder.exists() || !folder.isDirectory()) {
                JOptionPane.showMessageDialog(this,
                        "Folder does not exist:\n" + path,
                        "Folder not found", JOptionPane.ERROR_MESSAGE);
                return;
            }
            try {
                Desktop.getDesktop().open(folder);
            } catch (Exception ex) {
                JOptionPane.showMessageDialog(this,
                        "Could not open folder: " + ex.getMessage(),
                        "Open error", JOptionPane.ERROR_MESSAGE);
            }
        });
        panel.add(openButton);

        return panel;
    }

    /**
     * Auto-derive the entrapment DB output file from the chosen folder, the
     * input protein database name, and the selected level. Mirrors the naming
     * convention used in the FDRBench README:
     * <ul>
     *   <li>protein level → {@code <db>_entrapment_pro.fasta}</li>
     *   <li>peptide level → {@code <db>_entrapment_pep.txt}</li>
     * </ul>
     * If no input DB is set yet the prefix falls back to {@code entrapment}.
     */
    private File deriveDbOutputFile() {
        String folder = dbOutputField != null ? dbOutputField.getText().trim() : "";
        if (folder.isEmpty()) {
            return null;
        }
        String level = dbLevelCombo != null
                ? String.valueOf(dbLevelCombo.getSelectedItem())
                : "peptide";
        String dbBase = "entrapment";
        if (dbFileField != null) {
            String dbPath = dbFileField.getText().trim();
            if (!dbPath.isEmpty()) {
                String name = new File(dbPath).getName();
                int dot = name.lastIndexOf('.');
                dbBase = (dot > 0 ? name.substring(0, dot) : name) + "_entrapment";
            }
        }
        String fileName = "protein".equals(level)
                ? dbBase + "_pro.fasta"
                : dbBase + "_pep.txt";
        return new File(folder, fileName);
    }

    /**
     * Resolve the comma-separated path string to forward via {@code -ms}.
     * Prefers the multi-file backing list (populated by Browse / Folder /
     * the edit dialog); falls back to the field's raw text — which the user
     * may have typed manually in single-file mode.
     */
    private String foreignSpeciesArgValue() {
        if (!foreignSpeciesFiles.isEmpty()) {
            return String.join(",", foreignSpeciesFiles);
        }
        return foreignSpeciesField != null ? foreignSpeciesField.getText() : "";
    }

    /**
     * Auto-derive the FDP-Estimation output file from the chosen folder, the
     * input PSM/peptide/protein file name, and the selected calculation level.
     * Returns {@code <folder>/<input>_fdp_<level>.csv}, or {@code null} if
     * no output folder is set yet.
     */
    private File deriveFdpOutputFile() {
        String folder = fdpOutputField != null ? fdpOutputField.getText().trim() : "";
        if (folder.isEmpty()) {
            return null;
        }
        String level = fdpLevelCombo != null
                ? String.valueOf(fdpLevelCombo.getSelectedItem())
                : "precursor";
        String inputBase = "fdp";
        if (inputFileField != null) {
            String inputPath = inputFileField.getText().trim();
            if (!inputPath.isEmpty()) {
                String name = new File(inputPath).getName();
                int dot = name.lastIndexOf('.');
                inputBase = (dot > 0 ? name.substring(0, dot) : name);
            }
        }
        return new File(folder, inputBase + "_fdp_" + level + ".csv");
    }

    /**
     * Browse + View buttons for an FDP-Estimation file row. The View button
     * loads the file's first {@value PREVIEW_ROWS} rows into a JTable so the
     * user can confirm the file's columns/structure without leaving the GUI.
     */
    private JPanel createBrowseAndViewPanel(JTextField targetField, String extensions) {
        JPanel panel = new JPanel(new GridLayout(1, 0, 5, 0));
        panel.setOpaque(false);
        panel.add(createBrowseButton(targetField, extensions));

        JButton viewButton = new JButton("View");
        styleButton(viewButton);
        viewButton.setToolTipText("Show the first " + PREVIEW_ROWS
                + " rows of the file in a table");
        viewButton.addActionListener(e -> showFilePreview(targetField.getText().trim()));
        panel.add(viewButton);
        return panel;
    }

    private static final int PREVIEW_ROWS = 100;

    private void showFilePreview(String path) {
        if (path.isEmpty()) {
            JOptionPane.showMessageDialog(this,
                    "Please select a file first.",
                    "No file", JOptionPane.WARNING_MESSAGE);
            return;
        }
        File file = new File(path);
        if (!file.exists()) {
            JOptionPane.showMessageDialog(this,
                    "File does not exist:\n" + path,
                    "File not found", JOptionPane.ERROR_MESSAGE);
            return;
        }

        String[] header;
        java.util.List<String[]> topRows = new ArrayList<>();
        // Sliding window of the last PREVIEW_ROWS lines — bounded memory even
        // for very large files (only 2N rows kept around).
        java.util.Deque<String[]> tailWindow = new java.util.ArrayDeque<>();
        long totalDataRows = 0;
        try (BufferedReader r = new BufferedReader(
                new InputStreamReader(new java.io.FileInputStream(file),
                        java.nio.charset.StandardCharsets.UTF_8))) {
            String line = r.readLine();
            if (line == null) {
                JOptionPane.showMessageDialog(this,
                        "File is empty.",
                        "Empty file", JOptionPane.WARNING_MESSAGE);
                return;
            }
            char delim = line.indexOf('\t') >= 0 ? '\t'
                    : (line.indexOf(',') >= 0 ? ',' : '\t');
            header = splitDelim(line, delim);
            while ((line = r.readLine()) != null) {
                String[] cells = splitDelim(line, delim);
                if (topRows.size() < PREVIEW_ROWS) {
                    topRows.add(cells);
                } else {
                    if (tailWindow.size() == PREVIEW_ROWS) {
                        tailWindow.pollFirst();
                    }
                    tailWindow.addLast(cells);
                }
                totalDataRows++;
            }
        } catch (IOException ex) {
            JOptionPane.showMessageDialog(this,
                    "Could not read file: " + ex.getMessage(),
                    "Read error", JOptionPane.ERROR_MESSAGE);
            return;
        }

        boolean truncated = !tailWindow.isEmpty();
        int totalDisplay = topRows.size() + (truncated ? 1 + tailWindow.size() : 0);
        String[][] data = new String[totalDisplay][header.length];
        int dataRow = 0;
        for (String[] row : topRows) {
            for (int j = 0; j < header.length; j++) {
                data[dataRow][j] = j < row.length ? row[j] : "";
            }
            dataRow++;
        }
        if (truncated) {
            for (int j = 0; j < header.length; j++) {
                data[dataRow][j] = "…";
            }
            dataRow++;
            for (String[] row : tailWindow) {
                for (int j = 0; j < header.length; j++) {
                    data[dataRow][j] = j < row.length ? row[j] : "";
                }
                dataRow++;
            }
        }

        JTable table = new JTable(data, header) {
            @Override public boolean isCellEditable(int r, int c) { return false; }
        };
        table.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        table.setFillsViewportHeight(true);
        table.getTableHeader().setReorderingAllowed(false);
        autoSizeTableColumns(table);

        JScrollPane sp = new JScrollPane(table);
        sp.setPreferredSize(new Dimension(900, 480));
        sp.setBorder(BorderFactory.createEmptyBorder());

        String shown = truncated
                ? "showing first " + topRows.size() + " + last " + tailWindow.size()
                : "showing all " + topRows.size();
        String title = "Preview — " + file.getName()
                + " (" + shown + " of " + totalDataRows + " data rows)";
        JOptionPane.showMessageDialog(this, sp, title, JOptionPane.PLAIN_MESSAGE);
    }

    /** Auto-size each column to the widest of its header text and visible cells, capped. */
    private static void autoSizeTableColumns(JTable table) {
        final int minWidth = 60;
        final int maxWidth = 280;
        FontMetrics fm = table.getFontMetrics(table.getFont());
        FontMetrics hfm = table.getTableHeader().getFontMetrics(table.getTableHeader().getFont());
        for (int c = 0; c < table.getColumnCount(); c++) {
            int width = hfm.stringWidth(String.valueOf(table.getColumnName(c))) + 24;
            for (int r = 0; r < table.getRowCount(); r++) {
                Object val = table.getValueAt(r, c);
                if (val != null) {
                    width = Math.max(width, fm.stringWidth(val.toString()) + 16);
                    if (width >= maxWidth) break;
                }
            }
            table.getColumnModel().getColumn(c)
                    .setPreferredWidth(Math.max(minWidth, Math.min(maxWidth, width)));
        }
    }

    /** Minimal CSV/TSV splitter — handles double-quoted fields. */
    private static String[] splitDelim(String line, char delim) {
        java.util.List<String> out = new ArrayList<>();
        StringBuilder cur = new StringBuilder();
        boolean inQuotes = false;
        for (int i = 0; i < line.length(); i++) {
            char c = line.charAt(i);
            if (c == '"') {
                if (inQuotes && i + 1 < line.length() && line.charAt(i + 1) == '"') {
                    cur.append('"');
                    i++;
                } else {
                    inQuotes = !inQuotes;
                }
            } else if (c == delim && !inQuotes) {
                out.add(cur.toString());
                cur.setLength(0);
            } else {
                cur.append(c);
            }
        }
        out.add(cur.toString());
        return out.toArray(new String[0]);
    }

    private JButton createBrowseButton(JTextField targetField, String extensions) {
        JButton button = new JButton("Browse");
        styleButton(button);
        button.addActionListener(e -> {
            JFileChooser chooser = new JFileChooser();
            String lastDir = prefs.get(PREF_LAST_DIR, System.getProperty("user.home"));
            chooser.setCurrentDirectory(new File(lastDir));

            String[] exts = extensions.split(",");
            FileNameExtensionFilter filter = new FileNameExtensionFilter(
                    extensions.toUpperCase() + " Files", exts);
            chooser.setFileFilter(filter);

            if (chooser.showOpenDialog(this) == JFileChooser.APPROVE_OPTION) {
                File file = chooser.getSelectedFile();
                targetField.setText(file.getAbsolutePath());
                prefs.put(PREF_LAST_DIR, file.getParent());
            }
        });
        return button;
    }

    /**
     * Multi-file Browse + Folder buttons that drive a backing list of paths.
     * <ul>
     *   <li><b>Browse</b> opens a multi-select file chooser; the picked
     *       paths fill {@code fileList} and the field switches into a
     *       hyperlink "(N files selected)" mode (single-click opens a
     *       dialog to view / edit / delete entries).</li>
     *   <li><b>Folder</b> picks a directory; the field shows the directory
     *       path verbatim (the underlying CLI knows to expand it).</li>
     *   <li>One file selected → field shows that single path, editable.</li>
     * </ul>
     */
    private JPanel createMultiFileButtonsPanel(JTextField targetField,
                                               java.util.List<String> fileList,
                                               String[] extensions) {
        JPanel buttons = new JPanel(new GridLayout(1, 0, 5, 0));
        buttons.setOpaque(false);

        JButton browse = new JButton("Browse");
        styleButton(browse);
        browse.setToolTipText("Select one or more " + String.join("/", extensions)
                + " files");
        browse.addActionListener(e -> {
            JFileChooser chooser = new JFileChooser();
            chooser.setMultiSelectionEnabled(true);
            chooser.setFileFilter(new FileNameExtensionFilter(
                    String.join("/", extensions).toUpperCase(Locale.ROOT) + " Files",
                    extensions));
            String lastDir = prefs.get(PREF_LAST_DIR, System.getProperty("user.home"));
            chooser.setCurrentDirectory(new File(lastDir));
            if (chooser.showOpenDialog(this) != JFileChooser.APPROVE_OPTION) {
                return;
            }
            File[] files = chooser.getSelectedFiles();
            if (files == null || files.length == 0) {
                return;
            }
            fileList.clear();
            for (File f : files) {
                fileList.add(f.getAbsolutePath());
            }
            updateFileFieldState(targetField, fileList);
            if (files[0].getParent() != null) {
                prefs.put(PREF_LAST_DIR, files[0].getParent());
            }
        });
        buttons.add(browse);

        JButton folder = new JButton("Folder");
        styleButton(folder);
        folder.setToolTipText("Select a folder containing "
                + String.join("/", extensions) + " files");
        folder.addActionListener(e -> {
            JFileChooser chooser = new JFileChooser();
            chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
            String lastDir = prefs.get(PREF_LAST_DIR, System.getProperty("user.home"));
            chooser.setCurrentDirectory(new File(lastDir));
            if (chooser.showOpenDialog(this) != JFileChooser.APPROVE_OPTION) {
                return;
            }
            File dir = chooser.getSelectedFile();
            // Expand to all matching files in the folder so the user can see
            // exactly what gets used (and edit the list via the hyperlink).
            fileList.clear();
            File[] children = dir.listFiles();
            if (children != null) {
                for (File f : children) {
                    if (!f.isFile()) continue;
                    String name = f.getName().toLowerCase(Locale.ROOT);
                    for (String ext : extensions) {
                        if (name.endsWith("." + ext.toLowerCase(Locale.ROOT))) {
                            fileList.add(f.getAbsolutePath());
                            break;
                        }
                    }
                }
            }
            if (fileList.isEmpty()) {
                // No matching files — just show the folder path so the user
                // sees their choice was registered.
                targetField.setText(dir.getAbsolutePath());
            } else {
                updateFileFieldState(targetField, fileList);
            }
            prefs.put(PREF_LAST_DIR, dir.getAbsolutePath());
        });
        buttons.add(folder);

        setupMultiFileFieldInteraction(targetField, fileList, extensions);
        return buttons;
    }

    /**
     * Switch a text field between "single path / editable" and a
     * "(N files selected)" hyperlink summary based on the size of {@code files}.
     */
    private void updateFileFieldState(JTextField field, java.util.List<String> files) {
        java.util.Map<java.awt.font.TextAttribute, Object> attrs =
                new java.util.HashMap<>(field.getFont().getAttributes());
        if (files != null && files.size() > 1) {
            field.setEditable(false);
            field.setText("(" + files.size() + " files selected)");
            boolean dark = FlatLaf.isLafDark();
            field.setForeground(dark ? new Color(100, 180, 255) : Color.BLUE);
            field.setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
            attrs.put(java.awt.font.TextAttribute.UNDERLINE,
                    java.awt.font.TextAttribute.UNDERLINE_ON);
        } else {
            field.setForeground(UIManager.getColor("TextField.foreground"));
            field.setCursor(Cursor.getPredefinedCursor(Cursor.TEXT_CURSOR));
            field.setEditable(true);
            attrs.put(java.awt.font.TextAttribute.UNDERLINE, -1);
            if (files != null && files.size() == 1) {
                field.setText(files.get(0));
            }
        }
        field.setFont(field.getFont().deriveFont(attrs));
    }

    /**
     * Wire up:
     * <ul>
     *   <li>Single-click on a "(N files selected)" field → open the editor
     *       dialog so the user can view / edit / remove paths.</li>
     *   <li>Manual edits in single-file mode → desync the backing list so
     *       the user-typed value wins on the next Run.</li>
     * </ul>
     */
    private void setupMultiFileFieldInteraction(JTextField field,
                                                java.util.List<String> fileList,
                                                String[] extensions) {
        field.getDocument().addDocumentListener(new javax.swing.event.DocumentListener() {
            private void update() {
                if (!field.isEditable()) {
                    return;
                }
                if (fileList.isEmpty()) {
                    return;
                }
                String text = field.getText();
                if (fileList.size() > 1) {
                    fileList.clear();
                } else if (!text.equals(fileList.get(0))) {
                    fileList.clear();
                }
            }
            @Override public void insertUpdate(javax.swing.event.DocumentEvent e) { update(); }
            @Override public void removeUpdate(javax.swing.event.DocumentEvent e) { update(); }
            @Override public void changedUpdate(javax.swing.event.DocumentEvent e) { update(); }
        });
        field.addMouseListener(new java.awt.event.MouseAdapter() {
            @Override
            public void mouseClicked(java.awt.event.MouseEvent e) {
                if (!field.isEditable() && e.getClickCount() == 1) {
                    showMultiFileEditDialog(field, fileList, extensions);
                }
            }
        });
    }

    /**
     * Modal dialog that shows one path per line and lets the user edit / add
     * / remove entries. OK rebuilds {@code fileList} from the lines and
     * re-renders the field.
     */
    private void showMultiFileEditDialog(JTextField field,
                                         java.util.List<String> fileList,
                                         String[] extensions) {
        JDialog dialog = new JDialog(this,
                isRunning ? "View File List" : "Edit File List", true);
        dialog.setLayout(new BorderLayout(0, 8));
        JPanel pad = new JPanel(new BorderLayout());
        pad.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));
        JTextArea area = new JTextArea();
        area.setFont(new Font(Font.MONOSPACED, Font.PLAIN, 12));
        area.setEditable(!isRunning);
        for (String p : fileList) {
            area.append(p + "\n");
        }
        JScrollPane sp = new JScrollPane(area);
        sp.setPreferredSize(new Dimension(720, 360));
        pad.add(sp, BorderLayout.CENTER);
        dialog.add(pad, BorderLayout.CENTER);

        JPanel actions = new JPanel(new FlowLayout(FlowLayout.RIGHT, 8, 8));
        JButton ok = new JButton(isRunning ? "Close" : "OK");
        styleButton(ok);
        JButton cancel = new JButton("Cancel");
        styleButton(cancel);
        cancel.setVisible(!isRunning);
        actions.add(ok);
        actions.add(cancel);
        dialog.add(actions, BorderLayout.SOUTH);

        ok.addActionListener(ev -> {
            if (isRunning) {
                dialog.dispose();
                return;
            }
            java.util.List<String> newPaths = new java.util.ArrayList<>();
            for (String line : area.getText().split("\\R")) {
                String p = line.trim();
                if (!p.isEmpty()) {
                    newPaths.add(p);
                }
            }
            fileList.clear();
            fileList.addAll(newPaths);
            if (newPaths.isEmpty()) {
                field.setText("");
            }
            updateFileFieldState(field, fileList);
            dialog.dispose();
        });
        cancel.addActionListener(ev -> dialog.dispose());

        dialog.pack();
        dialog.setLocationRelativeTo(this);
        dialog.setVisible(true);
    }

    private JButton createMultiBrowseButton(JTextField targetField, String extensions) {
        JButton button = new JButton("Browse");
        styleButton(button);
        button.addActionListener(e -> {
            JFileChooser chooser = new JFileChooser();
            chooser.setMultiSelectionEnabled(true);
            String lastDir = prefs.get(PREF_LAST_DIR, System.getProperty("user.home"));
            chooser.setCurrentDirectory(new File(lastDir));

            String[] exts = extensions.split(",");
            FileNameExtensionFilter filter = new FileNameExtensionFilter(
                    extensions.toUpperCase() + " Files", exts);
            chooser.setFileFilter(filter);

            if (chooser.showOpenDialog(this) == JFileChooser.APPROVE_OPTION) {
                File[] files = chooser.getSelectedFiles();
                if (files != null && files.length > 0) {
                    StringBuilder joined = new StringBuilder();
                    for (File file : files) {
                        if (joined.length() > 0) {
                            joined.append(',');
                        }
                        joined.append(file.getAbsolutePath());
                    }
                    targetField.setText(joined.toString());
                    prefs.put(PREF_LAST_DIR, files[0].getParent());
                }
            }
        });
        return button;
    }

    private JButton createSaveButton(JTextField targetField) {
        JButton button = new JButton("Save As");
        styleButton(button);
        button.addActionListener(e -> {
            JFileChooser chooser = new JFileChooser() {
                @Override
                public void approveSelection() {
                    File f = getSelectedFile();
                    if (f != null && f.exists() && getDialogType() == SAVE_DIALOG) {
                        int answer = JOptionPane.showConfirmDialog(this,
                                "The file \"" + f.getName() + "\" already exists.\nDo you want to overwrite it?",
                                "Confirm Overwrite", JOptionPane.YES_NO_OPTION,
                                JOptionPane.WARNING_MESSAGE);
                        if (answer != JOptionPane.YES_OPTION) {
                            return;
                        }
                    }
                    super.approveSelection();
                }
            };
            String lastDir = prefs.get(PREF_LAST_DIR, System.getProperty("user.home"));
            chooser.setCurrentDirectory(new File(lastDir));

            if (chooser.showSaveDialog(this) == JFileChooser.APPROVE_OPTION) {
                File file = chooser.getSelectedFile();
                targetField.setText(file.getAbsolutePath());
                if (file.getParent() != null) {
                    prefs.put(PREF_LAST_DIR, file.getParent());
                }
            }
        });
        return button;
    }

    private JButton createPrimaryButton(String text, Color color) {
        JButton button = new JButton(text);
        button.setFont(button.getFont().deriveFont(Font.BOLD, 14f));
        button.setFocusPainted(false);
        button.setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
        // Tightened from 30→18 horizontal so the six footer buttons (incl.
        // Help) fit on one row at the default 700px window width without
        // wrapping. FlatLaf's roundRect rendering already includes its own
        // internal padding, so 18px reads as comfortable.
        button.setMargin(new Insets(12, 18, 12, 18));
        // FlatLaf-native styling: client properties paint a rounded button with
        // the brand color and white text instead of fighting the LAF with
        // setOpaque/setBorderPainted/setBackground (which left the button looking
        // flat and out of place after theme switches).
        button.putClientProperty("JButton.buttonType", "roundRect");
        button.putClientProperty("JButton.background", color);
        button.putClientProperty("JButton.foreground", Color.WHITE);
        return button;
    }

    private JButton createSecondaryButton(String text) {
        JButton button = new JButton(text);
        button.setFont(button.getFont().deriveFont(Font.PLAIN, 13f));
        button.setFocusPainted(false);
        button.setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
        // Tightened from 20→14 horizontal for the same reason as the primary
        // buttons (see createPrimaryButton) — keeps the six footer buttons
        // on one row at 700px window width.
        button.setMargin(new Insets(10, 14, 10, 14));
        button.putClientProperty("JButton.buttonType", "roundRect");
        return button;
    }

    // ==================== THEME HELPERS ====================

    private void toggleDarkMode(boolean isDark) {
        try {
            if (isDark)
                FlatDarkLaf.setup();
            else
                FlatLightLaf.setup();
            customizeUIDefaults();
            prefs.putBoolean(PREF_DARK_MODE, isDark);
            SwingUtilities.updateComponentTreeUI(this);
            updateHeaderForegrounds();
        } catch (Exception e) {
            System.err.println("Theme switch failed: " + e.getMessage());
        }
    }

    private static Color adjust(Color c, int delta) {
        return new Color(
                Math.max(0, Math.min(255, c.getRed() + delta)),
                Math.max(0, Math.min(255, c.getGreen() + delta)),
                Math.max(0, Math.min(255, c.getBlue() + delta)));
    }

    private static Color withAlpha(Color c, int a) {
        return new Color(c.getRed(), c.getGreen(), c.getBlue(), Math.max(0, Math.min(255, a)));
    }

    private static Color pickOnColor(Color bg) {
        double r = bg.getRed() / 255.0;
        double g = bg.getGreen() / 255.0;
        double b = bg.getBlue() / 255.0;
        r = (r <= 0.03928) ? (r / 12.92) : Math.pow((r + 0.055) / 1.055, 2.4);
        g = (g <= 0.03928) ? (g / 12.92) : Math.pow((g + 0.055) / 1.055, 2.4);
        b = (b <= 0.03928) ? (b / 12.92) : Math.pow((b + 0.055) / 1.055, 2.4);
        double L = 0.2126 * r + 0.7152 * g + 0.0722 * b;
        return (L > 0.70) ? new Color(20, 20, 20) : Color.WHITE;
    }

    // ==================== MAIN ====================

    public static void main(String[] args) {
        System.setProperty("awt.useSystemAAFontSettings", "on");
        System.setProperty("swing.aatext", "true");

        // Install FlatLaf BEFORE constructing the JFrame. JFrame.frameInit()
        // checks LookAndFeel.getSupportsWindowDecorations() at construction
        // time — if the LAF is still default Metal at that moment (as happens
        // when setup is done inside the constructor), the FlatLaf-painted
        // title bar never gets applied and the JFrame stays OS-decorated.
        // This is also what causes the title bar to revert to native chrome
        // after a theme toggle.
        try {
            boolean dark = prefs.getBoolean(PREF_DARK_MODE, false);
            if (dark) {
                FlatDarkLaf.setup();
            } else {
                FlatLightLaf.setup();
            }
            customizeUIDefaults();
        } catch (Exception e) {
            System.err.println("Theme setup failed: " + e.getMessage());
        }
        JFrame.setDefaultLookAndFeelDecorated(true);
        JDialog.setDefaultLookAndFeelDecorated(true);

        SwingUtilities.invokeLater(() -> {
            FDRBenchGUI gui = new FDRBenchGUI();
            gui.setVisible(true);
        });
    }

    /**
     * FlowLayout that correctly reports its preferred and minimum sizes when
     * components wrap to additional rows. Plain {@link FlowLayout} always
     * reports a single-row height even when it lays out multiple rows, which
     * causes the wrapped row(s) to be clipped under {@link BorderLayout}.
     * This is the well-known Rob Camick {@code WrapLayout} idiom.
     */
    private static class WrapLayout extends FlowLayout {
        WrapLayout(int align, int hgap, int vgap) {
            super(align, hgap, vgap);
        }

        @Override
        public Dimension preferredLayoutSize(Container target) {
            return layoutSize(target, true);
        }

        @Override
        public Dimension minimumLayoutSize(Container target) {
            Dimension minimum = layoutSize(target, false);
            minimum.width -= (getHgap() + 1);
            return minimum;
        }

        private Dimension layoutSize(Container target, boolean preferred) {
            synchronized (target.getTreeLock()) {
                int targetWidth = target.getSize().width;
                Container container = target;
                while (container.getSize().width == 0 && container.getParent() != null) {
                    container = container.getParent();
                }
                targetWidth = container.getSize().width;
                if (targetWidth == 0) {
                    targetWidth = Integer.MAX_VALUE;
                }

                int hgap = getHgap();
                int vgap = getVgap();
                Insets insets = target.getInsets();
                int horizontalInsetsAndGap = insets.left + insets.right + (hgap * 2);
                int maxWidth = targetWidth - horizontalInsetsAndGap;

                Dimension dim = new Dimension(0, 0);
                int rowWidth = 0;
                int rowHeight = 0;

                int n = target.getComponentCount();
                for (int i = 0; i < n; i++) {
                    Component m = target.getComponent(i);
                    if (m.isVisible()) {
                        Dimension d = preferred ? m.getPreferredSize() : m.getMinimumSize();
                        if (rowWidth + d.width > maxWidth) {
                            addRow(dim, rowWidth, rowHeight);
                            rowWidth = 0;
                            rowHeight = 0;
                        }
                        if (rowWidth != 0) {
                            rowWidth += hgap;
                        }
                        rowWidth += d.width;
                        rowHeight = Math.max(rowHeight, d.height);
                    }
                }
                addRow(dim, rowWidth, rowHeight);
                dim.width += horizontalInsetsAndGap;
                dim.height += insets.top + insets.bottom + vgap * 2;

                return dim;
            }
        }

        private void addRow(Dimension dim, int rowWidth, int rowHeight) {
            dim.width = Math.max(dim.width, rowWidth);
            if (dim.height > 0) {
                dim.height += getVgap();
            }
            dim.height += rowHeight;
        }
    }
}
