package main.java.gui;

import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.annotations.XYLineAnnotation;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.NumberTickUnit;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.ValueMarker;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.chart.title.LegendTitle;
import org.jfree.chart.ui.RectangleEdge;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import de.rototor.pdfbox.graphics2d.PdfBoxGraphics2D;
import org.apache.pdfbox.pdmodel.PDDocument;
import org.apache.pdfbox.pdmodel.PDPage;
import org.apache.pdfbox.pdmodel.PDPageContentStream;
import org.apache.pdfbox.pdmodel.common.PDRectangle;
import org.apache.pdfbox.pdmodel.graphics.form.PDFormXObject;

import javax.imageio.IIOImage;
import javax.imageio.ImageIO;
import javax.imageio.ImageTypeSpecifier;
import javax.imageio.ImageWriteParam;
import javax.imageio.ImageWriter;
import javax.imageio.metadata.IIOMetadata;
import javax.imageio.metadata.IIOMetadataNode;
import javax.imageio.stream.ImageOutputStream;
import javax.swing.*;
import javax.swing.filechooser.FileNameExtensionFilter;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Locale;
import java.util.Map;

/**
 * Plot tab for the FDP-Estimation workflow. Reads the FDP CSV produced by
 * FDRBench and renders an FDR-vs-FDP line chart that mirrors the R
 * {@code plot_fdp_fdr()} function (combined / paired / lower-bound methods,
 * y=x reference, vertical line at FDR=1%, percent-formatted axes).
 */
public class FdpPlotPanel extends JPanel {

    // Color palette mirroring the R script's color_mapping.
    private static final Color COLOR_COMBINED  = new Color(0xF8766D);
    private static final Color COLOR_PAIRED    = new Color(0x7CAE00);
    private static final Color COLOR_LOWER     = new Color(0x00BFC4);
    private static final Color COLOR_REFERENCE = new Color(0xAAAAAA);
    private static final Color COLOR_FDR_LINE  = new Color(0x1F77FF);

    // Default chart font. Java's Font constructor silently falls back to the
    // platform default sans-serif if Arial isn't installed, so this is safe
    // to use on any OS — Windows ships Arial, macOS has it bundled, and most
    // Linux distros include it via msttcorefonts or a Liberation alias.
    private static final String CHART_FONT = "Arial";
    private static final Font AXIS_LABEL_FONT = new Font(CHART_FONT, Font.PLAIN, 13);
    private static final Font TICK_FONT       = new Font(CHART_FONT, Font.PLAIN, 12);
    private static final Font LEGEND_FONT     = new Font(CHART_FONT, Font.PLAIN, 12);
    private static final Font ANNOTATION_FONT = new Font(CHART_FONT, Font.PLAIN, 11);

    private final JLabel placeholder;
    private final JButton saveButton;
    private final JTextField fileField;
    private ChartPanel chartPanel;
    private JComponent chartContainer;
    private JFreeChart chart;
    private File loadedCsv;

    // Plot-side controls (right-hand panel). Re-created on every CSV load, so
    // they're nullable — null whenever no chart is showing.
    private JRadioButton xAutoBtn;
    private JRadioButton xCustomBtn;
    private JSpinner     xCustomValue;
    private JRadioButton yAutoBtn;
    private JRadioButton yCustomBtn;
    private JSpinner     yCustomValue;
    private Runnable     applyXRange;
    private Runnable     applyYRange;

    // Discovery annotation: text content (fixed for a given CSV) and the
    // XYTextAnnotation handles currently attached to the chart. We rebuild
    // and reposition the latter every time an axis range changes so the
    // text block stays inside the visible region.
    private java.util.List<AnnotationLine> annotationLines;
    private java.util.List<XYTextAnnotation> textAnnotations;

    // Live visibility of each method, mirrored from the checkboxes so the
    // discovery annotation can hide / show its per-method lines in sync.
    private boolean showCombined = true;
    private boolean showPaired   = true;
    private boolean showLower    = true;

    // Re-entry guard for the AxisChangeListener: setTickUnit /
    // setNumberFormatOverride inside the listener fire fresh AxisChangeEvents
    // that would otherwise loop back into the listener.
    private boolean axisUpdating;

    // Whether the current FDP CSV actually contains each upper-bound method.
    // Needed alongside showCombined/Paired to decide if exactly one upper-bound
    // series is currently on the chart — which triggers an automatic rename
    // to "Upper bound" in both the legend and the discovery annotation.
    private boolean dataHasCombined;
    private boolean dataHasPaired;

    /** One line of the discovery-text block, tagged so we can filter by method. */
    private static class AnnotationLine {
        static final String TOTAL = "total";
        static final String COMBINED = "combined";
        static final String PAIRED   = "paired";
        static final String LOWER    = "lower";
        final String text;
        final String method;
        AnnotationLine(String text, String method) {
            this.text = text;
            this.method = method;
        }
    }

    public FdpPlotPanel() {
        super(new BorderLayout(0, 8));
        setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));

        placeholder = new JLabel(
                "<html><div style='text-align:center;padding:32px;'>"
                + "Run an FDP estimation to populate the plot,<br/>"
                + "or browse to an existing FDP CSV above."
                + "</div></html>",
                SwingConstants.CENTER);
        placeholder.setForeground(UIManager.getColor("Label.disabledForeground"));
        add(placeholder, BorderLayout.CENTER);

        // ---- Top row: label / field / Browse / View / Save Image, all
        // arranged with GridBagLayout so the field stays at its preferred
        // height (COMPONENT_HEIGHT) instead of being vertically stretched by
        // BorderLayout.CENTER to match the buttons' natural height — which
        // is what made the FDP file box look taller than the Workflow tab's
        // input boxes.
        JPanel topBar = new JPanel(new GridBagLayout());
        topBar.setOpaque(false);
        GridBagConstraints g = new GridBagConstraints();
        g.gridy = 0;
        g.insets = new Insets(0, 0, 0, 5);
        g.anchor = GridBagConstraints.CENTER;

        JLabel fileLabel = new JLabel("FDP file:");
        g.gridx = 0;
        g.weightx = 0;
        g.fill = GridBagConstraints.NONE;
        g.insets = new Insets(0, 0, 0, 8);
        topBar.add(fileLabel, g);

        fileField = new JTextField();
        // Width is preferred-only; height is pinned to COMPONENT_HEIGHT so it
        // matches every other input box in the GUI.
        fileField.setPreferredSize(new Dimension(200, FDRBenchGUI.COMPONENT_HEIGHT));
        fileField.putClientProperty("JTextField.placeholderText",
                "Path to an FDP CSV (auto-filled after a successful FDP run)");
        fileField.setToolTipText("FDP estimation CSV used to render the plot");
        fileField.addActionListener(e -> loadFromTypedPath());
        g.gridx = 1;
        g.weightx = 1;
        g.fill = GridBagConstraints.HORIZONTAL;
        g.insets = new Insets(0, 0, 0, 5);
        topBar.add(fileField, g);

        JButton browseButton = stylePlotButton(new JButton("Browse"));
        browseButton.setToolTipText("Choose an FDP CSV to plot");
        browseButton.addActionListener(e -> browseForFdpFile());
        g.gridx = 2;
        g.weightx = 0;
        g.fill = GridBagConstraints.NONE;
        topBar.add(browseButton, g);

        JButton viewButton = stylePlotButton(new JButton("View"));
        viewButton.setToolTipText("Show the first " + PREVIEW_ROWS
                + " rows of the file in a table");
        viewButton.addActionListener(e -> showFilePreview());
        g.gridx = 3;
        topBar.add(viewButton, g);

        saveButton = stylePlotButton(new JButton("Save Image…"));
        saveButton.setEnabled(false);
        saveButton.addActionListener(e -> saveImage());
        g.gridx = 4;
        g.insets = new Insets(0, 8, 0, 0);
        topBar.add(saveButton, g);

        add(topBar, BorderLayout.NORTH);
    }

    private static JButton stylePlotButton(JButton b) {
        b.setFocusPainted(false);
        b.setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
        b.setMargin(new Insets(6, 12, 6, 12));
        b.putClientProperty("JButton.buttonType", "roundRect");
        return b;
    }

    private void loadFromTypedPath() {
        String path = fileField.getText().trim();
        if (path.isEmpty()) {
            return;
        }
        File f = new File(path);
        if (!f.exists()) {
            JOptionPane.showMessageDialog(this,
                    "File does not exist:\n" + path,
                    "File not found", JOptionPane.ERROR_MESSAGE);
            return;
        }
        loadFromCsv(f);
    }

    private void browseForFdpFile() {
        JFileChooser chooser = new JFileChooser();
        chooser.setFileFilter(new FileNameExtensionFilter("FDP CSV/TSV", "csv", "tsv", "txt"));
        if (loadedCsv != null && loadedCsv.getParentFile() != null) {
            chooser.setCurrentDirectory(loadedCsv.getParentFile());
        } else if (!fileField.getText().trim().isEmpty()) {
            File parent = new File(fileField.getText().trim()).getParentFile();
            if (parent != null && parent.exists()) {
                chooser.setCurrentDirectory(parent);
            }
        }
        if (chooser.showOpenDialog(this) != JFileChooser.APPROVE_OPTION) {
            return;
        }
        File f = chooser.getSelectedFile();
        fileField.setText(f.getAbsolutePath());
        loadFromCsv(f);
    }

    private static final int PREVIEW_ROWS = 100;

    private void showFilePreview() {
        String path = fileField.getText().trim();
        if (path.isEmpty()) {
            JOptionPane.showMessageDialog(this,
                    "Please select an FDP CSV first.",
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
        // Sliding window of the last PREVIEW_ROWS lines — keeps memory bounded
        // even on very large files since we never hold more than 2N rows.
        java.util.Deque<String[]> tailWindow = new java.util.ArrayDeque<>();
        long totalDataRows = 0;
        try (BufferedReader r = new BufferedReader(
                new java.io.InputStreamReader(new java.io.FileInputStream(file),
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
            header = splitCsv(line, delim);
            while ((line = r.readLine()) != null) {
                String[] cells = splitCsv(line, delim);
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

        // If the whole file fits in 2N rows, no need for the separator — every
        // row that exists is in topRows already (tailWindow stays empty).
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

    private static void autoSizeTableColumns(JTable table) {
        final int minWidth = 60;
        final int maxWidth = 280;
        FontMetrics fm  = table.getFontMetrics(table.getFont());
        FontMetrics hfm = table.getTableHeader()
                .getFontMetrics(table.getTableHeader().getFont());
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

    /**
     * Replace the placeholder with a chart built from the given CSV. If the
     * file is missing or empty, leaves the placeholder in place.
     */
    public void loadFromCsv(File csv) {
        if (csv == null || !csv.exists()) {
            return;
        }
        try {
            FdpData data = parseCsv(csv);
            chart = buildChart(data);
            // Mirror axis-range changes from any source (mouse zoom, panel
            // controls, double-click reset) into our tick / format / annotation
            // refresh pipeline.
            final XYPlot plotForListeners = chart.getXYPlot();
            final NumberAxis xAxisForListener = (NumberAxis) plotForListeners.getDomainAxis();
            final NumberAxis yAxisForListener = (NumberAxis) plotForListeners.getRangeAxis();
            xAxisForListener.addChangeListener(e -> onAxisChanged(xAxisForListener));
            yAxisForListener.addChangeListener(e -> onAxisChanged(yAxisForListener));
            // Compute the discovery-annotation text once; positions get
            // re-derived from the current axis bounds whenever those change,
            // and per-method lines are hidden / shown in sync with the
            // controls-panel checkboxes (all default to checked here).
            annotationLines = buildDiscoveryAnnotationLines(data);
            textAnnotations = null;
            showCombined = true;
            showPaired   = true;
            showLower    = true;
            dataHasCombined = data.hasCombined();
            dataHasPaired   = data.hasPaired();
            refreshLegendLabels();
            positionDiscoveryAnnotation();
            // Original axis bounds set by buildChart — captured here so a
            // double-click can restore them exactly. ChartPanel.restoreAutoBounds()
            // would instead auto-fit to the dataset, which collapses the x-axis
            // to the q-value range and pushes the y=x diagonal, the FDR=1%
            // marker, and the annotation block off-screen.
            final double initialMax = Math.max(data.maxFdp, 0.05);
            ChartPanel newPanel = new ChartPanel(chart, true, true, true, true, true);
            newPanel.setPreferredSize(new Dimension(640, 480));
            newPanel.setMouseWheelEnabled(true);

            JPanel controls = buildPlotControls(data, initialMax);

            newPanel.addMouseListener(new java.awt.event.MouseAdapter() {
                @Override
                public void mouseClicked(java.awt.event.MouseEvent e) {
                    if (e.getClickCount() == 2
                            && SwingUtilities.isLeftMouseButton(e)) {
                        // Flip both axes back to Auto so the controls panel
                        // and the chart stay in sync after the reset.
                        if (xAutoBtn != null) xAutoBtn.setSelected(true);
                        if (yAutoBtn != null) yAutoBtn.setSelected(true);
                        if (applyXRange != null) applyXRange.run();
                        if (applyYRange != null) applyYRange.run();
                    }
                }
            });

            JPanel container = new JPanel(new BorderLayout(0, 0));
            container.setOpaque(false);
            container.add(newPanel, BorderLayout.CENTER);
            container.add(controls, BorderLayout.EAST);

            if (chartContainer != null) {
                remove(chartContainer);
            } else {
                remove(placeholder);
            }
            chartPanel = newPanel;
            chartContainer = container;
            add(chartContainer, BorderLayout.CENTER);
            saveButton.setEnabled(true);
            loadedCsv = csv;
            // Keep the file field in sync so users see which CSV is being
            // plotted, regardless of whether it was set via Browse or
            // auto-loaded after an FDP run.
            if (fileField != null && !csv.getAbsolutePath().equals(fileField.getText())) {
                fileField.setText(csv.getAbsolutePath());
            }
            revalidate();
            repaint();
        } catch (IOException ex) {
            JOptionPane.showMessageDialog(this,
                    "Could not parse FDP CSV: " + ex.getMessage(),
                    "Plot error", JOptionPane.ERROR_MESSAGE);
        }
    }

    /**
     * Build the right-hand controls panel that drives axis ranges and per-method
     * visibility on the currently-loaded chart.
     */
    private JPanel buildPlotControls(FdpData data, double initialMax) {
        final XYPlot plot = chart.getXYPlot();
        final NumberAxis xAxis = (NumberAxis) plot.getDomainAxis();
        final NumberAxis yAxis = (NumberAxis) plot.getRangeAxis();
        final XYItemRenderer renderer = plot.getRenderer();
        final XYSeriesCollection ds = (XYSeriesCollection) plot.getDataset();
        final int idxCombined = ds.indexOf("Combined method");
        final int idxPaired   = ds.indexOf("Paired method");
        final int idxLower    = ds.indexOf("Lower bound");

        JPanel p = new JPanel(new GridBagLayout());
        p.setBorder(BorderFactory.createCompoundBorder(
                BorderFactory.createMatteBorder(0, 1, 0, 0, new Color(0xDDDDDD)),
                BorderFactory.createEmptyBorder(10, 12, 10, 12)));
        p.setPreferredSize(new Dimension(180, 0));
        GridBagConstraints g = new GridBagConstraints();
        g.gridx = 0;
        g.gridy = 0;
        g.weightx = 1;
        g.fill = GridBagConstraints.HORIZONTAL;
        g.anchor = GridBagConstraints.NORTHWEST;
        g.insets = new Insets(0, 0, 6, 0);

        JLabel title = new JLabel("Plot settings");
        title.setFont(title.getFont().deriveFont(Font.BOLD));
        p.add(title, g);
        g.gridy++;

        // ---- X axis range ------------------------------------------------
        p.add(new JLabel("X axis range"), g);
        g.gridy++;
        xAutoBtn   = new JRadioButton("Auto", true);
        xCustomBtn = new JRadioButton("Custom:");
        ButtonGroup xGroup = new ButtonGroup();
        xGroup.add(xAutoBtn);
        xGroup.add(xCustomBtn);
        xCustomValue = new JSpinner(new SpinnerNumberModel(
                Math.max(0.001, Math.min(1.0, initialMax)), 0.001, 1.0, 0.001));
        xCustomValue.setEditor(new JSpinner.NumberEditor(xCustomValue, "0.000"));
        xCustomValue.setEnabled(false);

        g.insets = new Insets(0, 8, 2, 0);
        p.add(xAutoBtn, g);
        g.gridy++;
        JPanel xCustomRow = new JPanel(new BorderLayout(4, 0));
        xCustomRow.setOpaque(false);
        xCustomRow.add(xCustomBtn,   BorderLayout.WEST);
        xCustomRow.add(xCustomValue, BorderLayout.CENTER);
        p.add(xCustomRow, g);
        g.gridy++;
        g.insets = new Insets(8, 0, 6, 0);

        applyXRange = () -> {
            boolean auto = xAutoBtn.isSelected();
            xCustomValue.setEnabled(!auto);
            double v = auto ? initialMax : ((Number) xCustomValue.getValue()).doubleValue();
            double tick = niceTick(v);
            xAxis.setRange(0, v);
            xAxis.setTickUnit(new NumberTickUnit(tick));
            applyPercentFormat(xAxis, tick);
            positionDiscoveryAnnotation();
        };
        xAutoBtn.addActionListener(e -> applyXRange.run());
        xCustomBtn.addActionListener(e -> applyXRange.run());
        xCustomValue.addChangeListener(e -> {
            if (xCustomBtn.isSelected()) applyXRange.run();
        });

        // ---- Y axis range ------------------------------------------------
        p.add(new JLabel("Y axis range"), g);
        g.gridy++;
        yAutoBtn   = new JRadioButton("Auto", true);
        yCustomBtn = new JRadioButton("Custom:");
        ButtonGroup yGroup = new ButtonGroup();
        yGroup.add(yAutoBtn);
        yGroup.add(yCustomBtn);
        yCustomValue = new JSpinner(new SpinnerNumberModel(
                Math.max(0.001, Math.min(1.0, initialMax)), 0.001, 1.0, 0.001));
        yCustomValue.setEditor(new JSpinner.NumberEditor(yCustomValue, "0.000"));
        yCustomValue.setEnabled(false);

        g.insets = new Insets(0, 8, 2, 0);
        p.add(yAutoBtn, g);
        g.gridy++;
        JPanel yCustomRow = new JPanel(new BorderLayout(4, 0));
        yCustomRow.setOpaque(false);
        yCustomRow.add(yCustomBtn,   BorderLayout.WEST);
        yCustomRow.add(yCustomValue, BorderLayout.CENTER);
        p.add(yCustomRow, g);
        g.gridy++;
        g.insets = new Insets(8, 0, 6, 0);

        applyYRange = () -> {
            boolean auto = yAutoBtn.isSelected();
            yCustomValue.setEnabled(!auto);
            double v = auto ? initialMax : ((Number) yCustomValue.getValue()).doubleValue();
            double tick = niceTick(v);
            yAxis.setRange(0, v);
            yAxis.setTickUnit(new NumberTickUnit(tick));
            applyPercentFormat(yAxis, tick);
            positionDiscoveryAnnotation();
        };
        yAutoBtn.addActionListener(e -> applyYRange.run());
        yCustomBtn.addActionListener(e -> applyYRange.run());
        yCustomValue.addChangeListener(e -> {
            if (yCustomBtn.isSelected()) applyYRange.run();
        });

        // ---- Upper bound (combined / paired) ----------------------------
        if (data.hasCombined() || data.hasPaired()) {
            p.add(new JLabel("Upper bound"), g);
            g.gridy++;
            g.insets = new Insets(0, 8, 2, 0);
            if (data.hasCombined() && idxCombined >= 0) {
                JCheckBox cb = new JCheckBox("Combined", true);
                cb.addActionListener(e -> {
                    renderer.setSeriesVisible(idxCombined, cb.isSelected());
                    renderer.setSeriesVisibleInLegend(idxCombined, cb.isSelected());
                    showCombined = cb.isSelected();
                    refreshLegendLabels();
                    positionDiscoveryAnnotation();
                });
                p.add(cb, g);
                g.gridy++;
            }
            if (data.hasPaired() && idxPaired >= 0) {
                JCheckBox cb = new JCheckBox("Paired", true);
                cb.addActionListener(e -> {
                    renderer.setSeriesVisible(idxPaired, cb.isSelected());
                    renderer.setSeriesVisibleInLegend(idxPaired, cb.isSelected());
                    showPaired = cb.isSelected();
                    refreshLegendLabels();
                    positionDiscoveryAnnotation();
                });
                p.add(cb, g);
                g.gridy++;
            }

            g.insets = new Insets(8, 0, 6, 0);
        }

        // ---- Lower bound -------------------------------------------------
        if (data.hasLower() && idxLower >= 0) {
            JCheckBox cb = new JCheckBox("Lower bound", true);
            // Render as "Lower bound [box]" — the label sits left of the
            // indicator so this top-level toggle reads differently from the
            // indented Upper-bound sub-checkboxes above.
            cb.setHorizontalTextPosition(SwingConstants.LEFT);
            cb.addActionListener(e -> {
                renderer.setSeriesVisible(idxLower, cb.isSelected());
                renderer.setSeriesVisibleInLegend(idxLower, cb.isSelected());
                showLower = cb.isSelected();
                positionDiscoveryAnnotation();
            });
            p.add(cb, g);
            g.gridy++;
        }

        // Push remaining vertical slack to the bottom so controls hug the top.
        g.weighty = 1;
        g.fill = GridBagConstraints.BOTH;
        p.add(new JPanel() {{ setOpaque(false); }}, g);

        return p;
    }

    private void saveImage() {
        if (chart == null) {
            return;
        }
        ExportOptions opts = askExportOptions();
        if (opts == null) {
            return;
        }
        String ext = opts.format == ExportFormat.PDF ? "pdf" : "png";
        String displayName = opts.format == ExportFormat.PDF ? "PDF Document" : "PNG Image";
        JFileChooser chooser = new JFileChooser();
        chooser.setFileFilter(new FileNameExtensionFilter(displayName, ext));
        if (loadedCsv != null && loadedCsv.getParentFile() != null) {
            chooser.setCurrentDirectory(loadedCsv.getParentFile());
            String base = loadedCsv.getName().replaceFirst("\\.[^.]+$", "");
            chooser.setSelectedFile(new File(loadedCsv.getParentFile(), base + "." + ext));
        }
        if (chooser.showSaveDialog(this) != JFileChooser.APPROVE_OPTION) {
            return;
        }
        File out = chooser.getSelectedFile();
        if (!out.getName().toLowerCase(Locale.ROOT).endsWith("." + ext)) {
            out = new File(out.getAbsolutePath() + "." + ext);
        }
        try {
            if (opts.format == ExportFormat.PDF) {
                writePdf(chart, out, opts.widthIn, opts.heightIn);
            } else {
                BufferedImage img = renderChart(chart, opts);
                writePngWithDpi(img, out, opts.dpi);
            }
        } catch (IOException ex) {
            JOptionPane.showMessageDialog(this,
                    "Could not save " + ext.toUpperCase(Locale.ROOT) + ": " + ex.getMessage(),
                    "Save error", JOptionPane.ERROR_MESSAGE);
        }
    }

    /**
     * Write the chart as a vector PDF via PdfBoxGraphics2D — strokes, text,
     * and shapes stay vector so the file scales without pixelation. The page
     * is sized in points (1 in = 72 pt).
     */
    private static void writePdf(JFreeChart chart, File out,
                                 double widthIn, double heightIn) throws IOException {
        float widthPt  = (float) (widthIn  * 72.0);
        float heightPt = (float) (heightIn * 72.0);
        try (PDDocument doc = new PDDocument()) {
            PDPage page = new PDPage(new PDRectangle(widthPt, heightPt));
            doc.addPage(page);
            PdfBoxGraphics2D pdfG2 = new PdfBoxGraphics2D(doc, widthPt, heightPt);
            try {
                pdfG2.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
                        RenderingHints.VALUE_ANTIALIAS_ON);
                pdfG2.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING,
                        RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
                chart.draw(pdfG2,
                        new java.awt.geom.Rectangle2D.Double(0, 0, widthPt, heightPt));
            } finally {
                pdfG2.dispose();
            }
            PDFormXObject xform = pdfG2.getXFormObject();
            try (PDPageContentStream cs = new PDPageContentStream(doc, page)) {
                cs.drawForm(xform);
            }
            doc.save(out);
        }
    }

    /**
     * Render the chart into a high-resolution PNG-sized buffer.
     *
     * <p>JFreeChart sizes fonts and strokes in absolute pixel units, so simply
     * asking it to draw into a 3000×3000 image makes a 12pt label still 12
     * pixels tall — it looks tiny relative to the bigger canvas, not sharper.
     * To get genuine high-resolution output we lay out the chart at a fixed
     * <em>logical</em> canvas (the figure's physical size at a 96-DPI baseline)
     * and apply a {@link Graphics2D#scale} so every drawing op is multiplied up
     * to the requested pixel dimensions. Result: a chart that looks identical
     * at any DPI, just with proportionally more pixels and crisper text.
     */
    private static BufferedImage renderChart(JFreeChart chart, ExportOptions opts) {
        int pxW = opts.width;
        int pxH = opts.height;
        // Logical canvas at the screen-DPI baseline. Layout, fonts, and
        // strokes are sized against this, not against the final pixel buffer.
        double logicalW = opts.widthIn  * SCREEN_DPI;
        double logicalH = opts.heightIn * SCREEN_DPI;
        double sx = pxW / logicalW;
        double sy = pxH / logicalH;

        BufferedImage img = new BufferedImage(pxW, pxH, BufferedImage.TYPE_INT_ARGB);
        Graphics2D g2 = img.createGraphics();
        try {
            g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
                    RenderingHints.VALUE_ANTIALIAS_ON);
            g2.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING,
                    RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
            g2.setRenderingHint(RenderingHints.KEY_FRACTIONALMETRICS,
                    RenderingHints.VALUE_FRACTIONALMETRICS_ON);
            g2.setRenderingHint(RenderingHints.KEY_STROKE_CONTROL,
                    RenderingHints.VALUE_STROKE_PURE);
            g2.scale(sx, sy);
            chart.draw(g2, new java.awt.geom.Rectangle2D.Double(0, 0, logicalW, logicalH));
        } finally {
            g2.dispose();
        }
        return img;
    }

    /** Logical DPI used for the chart's pre-scale layout. */
    private static final int SCREEN_DPI = 96;

    private enum ExportFormat { PNG, PDF }

    /**
     * Pixel dimensions, physical inches, DPI, and target file format for the
     * export. For PNG the (inches × DPI) pair drives the pixel count and the
     * pHYs metadata; for PDF the page is sized in inches directly and DPI is
     * irrelevant (vector output).
     */
    private static class ExportOptions {
        final int width;
        final int height;
        final int dpi;
        final double widthIn;
        final double heightIn;
        final ExportFormat format;
        ExportOptions(double widthIn, double heightIn, int dpi, ExportFormat format) {
            this.widthIn = widthIn;
            this.heightIn = heightIn;
            this.dpi = dpi;
            this.format = format;
            this.width  = Math.max(64, (int) Math.round(widthIn  * dpi));
            this.height = Math.max(64, (int) Math.round(heightIn * dpi));
        }
    }

    private ExportOptions askExportOptions() {
        // Width/height are in inches. For PNG: (inches × DPI) → pixels; for
        // PDF the page is sized in inches directly and DPI is hidden because
        // vector output is resolution-independent.
        JComboBox<String> formatCombo = new JComboBox<>(new String[] { "PNG", "PDF" });
        JSpinner widthSpinner  = new JSpinner(new SpinnerNumberModel(4.0, 0.5, 24.0, 0.1));
        JSpinner heightSpinner = new JSpinner(new SpinnerNumberModel(4.0, 0.5, 24.0, 0.1));
        JSpinner dpiSpinner    = new JSpinner(new SpinnerNumberModel(300,  72, 1200, 1));
        widthSpinner .setEditor(new JSpinner.NumberEditor(widthSpinner,  "0.00"));
        heightSpinner.setEditor(new JSpinner.NumberEditor(heightSpinner, "0.00"));

        JLabel dpiLabel = new JLabel("DPI:");
        JLabel pixelInfo = new JLabel(" ");
        Runnable refreshInfo = () -> {
            boolean isPdf = "PDF".equals(formatCombo.getSelectedItem());
            dpiLabel.setVisible(!isPdf);
            dpiSpinner.setVisible(!isPdf);
            double wIn = ((Number) widthSpinner.getValue()).doubleValue();
            double hIn = ((Number) heightSpinner.getValue()).doubleValue();
            if (isPdf) {
                pixelInfo.setText(String.format(Locale.US,
                        "Output: %.2f × %.2f in PDF (vector)", wIn, hIn));
            } else {
                int dpi  = ((Number) dpiSpinner.getValue()).intValue();
                int wpx = (int) Math.round(wIn * dpi);
                int hpx = (int) Math.round(hIn * dpi);
                pixelInfo.setText("Output: " + wpx + " × " + hpx + " px");
            }
        };
        // When the user toggles between PNG and PDF, snap the size spinners to
        // the format's default (PNG → 4 × 4 in, PDF → 4.5 × 4.5 in). Picked as
        // sensible defaults for proteomics figures: 4 in for raster screen
        // captures, 4.5 in for vector single-column publication figures.
        formatCombo.addActionListener(e -> {
            boolean isPdf = "PDF".equals(formatCombo.getSelectedItem());
            double sizeIn = isPdf ? 4.5 : 4.0;
            widthSpinner.setValue(sizeIn);
            heightSpinner.setValue(sizeIn);
            refreshInfo.run();
        });
        widthSpinner .addChangeListener(e -> refreshInfo.run());
        heightSpinner.addChangeListener(e -> refreshInfo.run());
        dpiSpinner   .addChangeListener(e -> refreshInfo.run());
        refreshInfo.run();

        JPanel form = new JPanel(new GridBagLayout());
        GridBagConstraints g = new GridBagConstraints();
        g.insets = new Insets(4, 6, 4, 6);
        g.fill = GridBagConstraints.HORIZONTAL;
        g.gridx = 0; g.gridy = 0; g.weightx = 0;
        form.add(new JLabel("Format:"), g);
        g.gridx = 1; g.weightx = 1;
        form.add(formatCombo, g);
        g.gridx = 0; g.gridy = 1; g.weightx = 0;
        form.add(new JLabel("Width (in):"), g);
        g.gridx = 1; g.weightx = 1;
        form.add(widthSpinner, g);
        g.gridx = 0; g.gridy = 2; g.weightx = 0;
        form.add(new JLabel("Height (in):"), g);
        g.gridx = 1; g.weightx = 1;
        form.add(heightSpinner, g);
        g.gridx = 0; g.gridy = 3; g.weightx = 0;
        form.add(dpiLabel, g);
        g.gridx = 1; g.weightx = 1;
        form.add(dpiSpinner, g);
        g.gridx = 0; g.gridy = 4; g.gridwidth = 2; g.weightx = 1;
        form.add(pixelInfo, g);

        int result = JOptionPane.showConfirmDialog(this, form,
                "Image export options",
                JOptionPane.OK_CANCEL_OPTION, JOptionPane.PLAIN_MESSAGE);
        if (result != JOptionPane.OK_OPTION) {
            return null;
        }
        double wIn = ((Number) widthSpinner.getValue()).doubleValue();
        double hIn = ((Number) heightSpinner.getValue()).doubleValue();
        int dpi  = ((Number) dpiSpinner.getValue()).intValue();
        ExportFormat fmt = "PDF".equals(formatCombo.getSelectedItem())
                ? ExportFormat.PDF : ExportFormat.PNG;
        return new ExportOptions(wIn, hIn, dpi, fmt);
    }

    /**
     * Write a PNG with the pHYs chunk populated so the file reports the
     * intended print resolution (pixels per meter). Most PNG writers leave
     * this empty, which makes Word/Photoshop default to 96 DPI regardless of
     * how big the image is.
     */
    private static void writePngWithDpi(BufferedImage img, File out, int dpi) throws IOException {
        Iterator<ImageWriter> writers = ImageIO.getImageWritersByFormatName("png");
        if (!writers.hasNext()) {
            throw new IOException("No PNG ImageWriter available");
        }
        ImageWriter writer = writers.next();
        try (ImageOutputStream ios = ImageIO.createImageOutputStream(out)) {
            writer.setOutput(ios);
            ImageWriteParam params = writer.getDefaultWriteParam();
            ImageTypeSpecifier typeSpecifier = ImageTypeSpecifier.createFromBufferedImageType(img.getType());
            IIOMetadata metadata = writer.getDefaultImageMetadata(typeSpecifier, params);
            int pixelsPerMeter = (int) Math.round(dpi * 39.3700787);
            String formatName = "javax_imageio_png_1.0";
            IIOMetadataNode pHYs = new IIOMetadataNode("pHYs");
            pHYs.setAttribute("pixelsPerUnitXAxis", Integer.toString(pixelsPerMeter));
            pHYs.setAttribute("pixelsPerUnitYAxis", Integer.toString(pixelsPerMeter));
            pHYs.setAttribute("unitSpecifier", "meter");
            IIOMetadataNode root = new IIOMetadataNode(formatName);
            root.appendChild(pHYs);
            metadata.mergeTree(formatName, root);
            writer.write(metadata, new IIOImage(img, null, metadata), params);
        } finally {
            writer.dispose();
        }
    }

    // ------------------------------------------------------------- CSV parsing

    private static class FdpData {
        final List<double[]> rows = new ArrayList<>();
        int qIdx = -1, combinedIdx = -1, pairedIdx = -1, lowerIdx = -1;
        double maxFdp = 0.0;

        boolean hasCombined() { return combinedIdx >= 0; }
        boolean hasPaired()   { return pairedIdx   >= 0; }
        boolean hasLower()    { return lowerIdx    >= 0; }
    }

    private static FdpData parseCsv(File csv) throws IOException {
        FdpData data = new FdpData();
        try (BufferedReader r = new BufferedReader(new FileReader(csv))) {
            String header = r.readLine();
            if (header == null) {
                throw new IOException("Empty CSV");
            }
            char delim = header.indexOf('\t') >= 0 ? '\t' : ',';
            String[] cols = splitCsv(header, delim);
            Map<String, Integer> idx = new HashMap<>();
            for (int i = 0; i < cols.length; i++) {
                idx.put(cols[i].trim().toLowerCase(Locale.ROOT), i);
            }
            data.qIdx        = idx.getOrDefault("q_value", -1);
            data.combinedIdx = idx.getOrDefault("combined_fdp", -1);
            data.pairedIdx   = idx.getOrDefault("paired_fdp", -1);
            data.lowerIdx    = idx.getOrDefault("lower_bound_fdp", -1);
            if (data.qIdx < 0) {
                throw new IOException("CSV is missing required q_value column");
            }
            if (!data.hasCombined() && !data.hasPaired() && !data.hasLower()) {
                throw new IOException("CSV has no *_fdp columns to plot");
            }

            String line;
            while ((line = r.readLine()) != null) {
                if (line.isEmpty()) continue;
                String[] tok = splitCsv(line, delim);
                double[] row = new double[tok.length];
                for (int i = 0; i < tok.length; i++) {
                    row[i] = parseDoubleSafe(tok[i]);
                }
                data.rows.add(row);
                data.maxFdp = Math.max(data.maxFdp, valueAt(row, data.qIdx));
                if (data.hasCombined()) data.maxFdp = Math.max(data.maxFdp, valueAt(row, data.combinedIdx));
                if (data.hasPaired())   data.maxFdp = Math.max(data.maxFdp, valueAt(row, data.pairedIdx));
                if (data.hasLower())    data.maxFdp = Math.max(data.maxFdp, valueAt(row, data.lowerIdx));
            }
        }
        if (data.rows.isEmpty()) {
            throw new IOException("CSV has no data rows");
        }
        if (data.maxFdp <= 0) {
            data.maxFdp = 1.0;
        }
        return data;
    }

    private static double valueAt(double[] row, int idx) {
        return (idx >= 0 && idx < row.length) ? row[idx] : Double.NaN;
    }

    private static double parseDoubleSafe(String s) {
        if (s == null) return Double.NaN;
        String t = s.trim();
        if (t.isEmpty() || "NA".equalsIgnoreCase(t) || "NaN".equalsIgnoreCase(t)) {
            return Double.NaN;
        }
        try {
            return Double.parseDouble(t);
        } catch (NumberFormatException e) {
            return Double.NaN;
        }
    }

    /** Minimal CSV/TSV splitter — handles double-quoted fields. */
    private static String[] splitCsv(String line, char delim) {
        List<String> out = new ArrayList<>();
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

    // -------------------------------------------------------- Chart construction

    private static JFreeChart buildChart(FdpData data) {
        XYSeriesCollection ds = new XYSeriesCollection();
        XYSeries combined = data.hasCombined() ? new XYSeries("Combined method") : null;
        XYSeries paired   = data.hasPaired()   ? new XYSeries("Paired method")   : null;
        XYSeries lower    = data.hasLower()    ? new XYSeries("Lower bound")     : null;

        for (double[] row : data.rows) {
            double q = valueAt(row, data.qIdx);
            if (Double.isNaN(q)) continue;
            if (combined != null) {
                double v = valueAt(row, data.combinedIdx);
                if (!Double.isNaN(v)) combined.add(q, v, false);
            }
            if (paired != null) {
                double v = valueAt(row, data.pairedIdx);
                if (!Double.isNaN(v)) paired.add(q, v, false);
            }
            if (lower != null) {
                double v = valueAt(row, data.lowerIdx);
                if (!Double.isNaN(v)) lower.add(q, v, false);
            }
        }

        // Series order matches the R legend: Combined, Paired, Lower.
        if (combined != null) ds.addSeries(combined);
        if (paired   != null) ds.addSeries(paired);
        if (lower    != null) ds.addSeries(lower);

        JFreeChart chart = org.jfree.chart.ChartFactory.createXYLineChart(
                null,
                "FDR threshold",
                "Estimated FDP",
                ds,
                PlotOrientation.VERTICAL,
                true,
                true,
                false);

        chart.setBackgroundPaint(Color.WHITE);

        XYPlot plot = chart.getXYPlot();
        plot.setBackgroundPaint(Color.WHITE);
        plot.setDomainGridlinePaint(new Color(0xE6E6E6));
        plot.setRangeGridlinePaint(new Color(0xE6E6E6));
        plot.setOutlinePaint(new Color(0xCCCCCC));

        XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer(true, false);
        BasicStroke lineStroke = new BasicStroke(1.6f);
        int seriesIdx = 0;
        if (combined != null) {
            renderer.setSeriesPaint(seriesIdx, COLOR_COMBINED);
            renderer.setSeriesStroke(seriesIdx, lineStroke);
            seriesIdx++;
        }
        if (paired != null) {
            renderer.setSeriesPaint(seriesIdx, COLOR_PAIRED);
            renderer.setSeriesStroke(seriesIdx, lineStroke);
            seriesIdx++;
        }
        if (lower != null) {
            renderer.setSeriesPaint(seriesIdx, COLOR_LOWER);
            renderer.setSeriesStroke(seriesIdx, lineStroke);
        }
        plot.setRenderer(renderer);

        // Cap axes at the largest observed FDP (or 5% if smaller, so the plot
        // doesn't collapse on perfectly-controlled inputs).
        double max = Math.max(data.maxFdp, 0.05);
        NumberAxis xAxis = (NumberAxis) plot.getDomainAxis();
        NumberAxis yAxis = (NumberAxis) plot.getRangeAxis();
        double initialTick = niceTick(max);
        xAxis.setRange(0, max);
        yAxis.setRange(0, max);
        xAxis.setTickUnit(new NumberTickUnit(initialTick));
        yAxis.setTickUnit(new NumberTickUnit(initialTick));
        applyPercentFormat(xAxis, initialTick);
        applyPercentFormat(yAxis, initialTick);
        xAxis.setLabelFont(AXIS_LABEL_FONT);
        yAxis.setLabelFont(AXIS_LABEL_FONT);
        xAxis.setTickLabelFont(TICK_FONT);
        yAxis.setTickLabelFont(TICK_FONT);

        // Diagonal y = x reference (gray) and vertical FDR=1% guide (blue, dashed).
        // Extend the diagonal to the full unit square so it stays visible even
        // when the user zooms the axes out past the initial data-derived max.
        plot.addAnnotation(new XYLineAnnotation(0, 0, 1.0, 1.0,
                new BasicStroke(1.0f), COLOR_REFERENCE));
        ValueMarker fdrMarker = new ValueMarker(0.01);
        fdrMarker.setPaint(COLOR_FDR_LINE);
        fdrMarker.setStroke(new BasicStroke(1.0f, BasicStroke.CAP_BUTT,
                BasicStroke.JOIN_MITER, 1.0f, new float[]{4f, 4f}, 0f));
        plot.addDomainMarker(fdrMarker);

        // Discovery text annotations are added by the caller (loadFromCsv) so
        // they can be repositioned later when the user changes the axis range.

        // Move the legend inside the plot at roughly the lower-right corner —
        // close to the R script's legend.position = c(0.7, 0.16).
        LegendTitle legend = chart.getLegend();
        if (legend != null) {
            legend.setPosition(RectangleEdge.BOTTOM);
            legend.setBackgroundPaint(new Color(255, 255, 255, 200));
            legend.setBorder(0, 0, 0, 0);
            legend.setItemFont(LEGEND_FONT);
        }

        return chart;
    }

    private static double niceTick(double range) {
        // Extra-fine entries cover mouse-rubber-band zooms into very tight
        // windows; without them a tick stride of 0.005 sits outside any
        // sub-0.005 range, the axis renders no ticks, and the labels vanish.
        if (range <= 0.0001) return 0.00002;
        if (range <= 0.0002) return 0.00005;
        if (range <= 0.0005) return 0.0001;
        if (range <= 0.001)  return 0.0002;
        if (range <= 0.002)  return 0.0005;
        if (range <= 0.005)  return 0.001;
        if (range <= 0.01)   return 0.002;
        if (range <= 0.02)   return 0.005;
        if (range <= 0.05)   return 0.01;
        if (range <= 0.10)   return 0.02;
        if (range <= 0.25)   return 0.05;
        if (range <= 0.50)   return 0.10;
        return 0.20;
    }

    /**
     * Fraction digits the percent formatter must show to render the given
     * tick unambiguously. A 0.005 tick (0.5%) would otherwise round to "0%".
     */
    private static int fractionDigitsForTick(double tick) {
        double pctTick = tick * 100.0;
        if (pctTick >= 1.0)     return 0;
        if (pctTick >= 0.1)     return 1;
        if (pctTick >= 0.01)    return 2;
        if (pctTick >= 0.001)   return 3;
        if (pctTick >= 0.0001)  return 4;
        return 5;
    }

    /**
     * Set the axis's percent formatter to use just enough fraction digits to
     * keep every tick label distinct from the next. Re-applied whenever the
     * range changes so zoomed-in axes (e.g. X max = 1%) don't collapse their
     * mid-axis tick to "0%".
     */
    private static void applyPercentFormat(NumberAxis axis, double tick) {
        NumberFormat fmt = NumberFormat.getPercentInstance(Locale.US);
        fmt.setMaximumFractionDigits(fractionDigitsForTick(tick));
        axis.setNumberFormatOverride(fmt);
    }

    /**
     * Called whenever an axis range changes — from the controls-panel
     * radios/spinners, the double-click reset, or JFreeChart's built-in
     * mouse-rubber-band zoom. Recomputes the tick stride and formatter
     * precision for the new range, and reflows the discovery annotation
     * so it stays anchored to the visible top-left.
     */
    private void onAxisChanged(NumberAxis axis) {
        if (axisUpdating) return;
        double span = axis.getUpperBound() - axis.getLowerBound();
        if (span <= 0) return;
        axisUpdating = true;
        try {
            double tick = niceTick(span);
            if (axis.getTickUnit() == null
                    || Math.abs(axis.getTickUnit().getSize() - tick) > 1e-12) {
                axis.setTickUnit(new NumberTickUnit(tick));
            }
            applyPercentFormat(axis, tick);
            positionDiscoveryAnnotation();
        } finally {
            axisUpdating = false;
        }
    }

    /**
     * Compute the per-method FDP text block reported in the chart annotation.
     * The lines are returned in display order, tagged by method so callers
     * can filter them when the matching series is toggled off in the panel.
     */
    private static List<AnnotationLine> buildDiscoveryAnnotationLines(FdpData data) {
        int discoveries = 0;
        double combinedAt001 = Double.NaN;
        double pairedAt001   = Double.NaN;
        double lowerAt001    = Double.NaN;
        double maxQAt001     = -1.0;
        for (double[] row : data.rows) {
            double q = valueAt(row, data.qIdx);
            if (Double.isNaN(q) || q > 0.01) continue;
            discoveries++;
            if (q > maxQAt001) {
                maxQAt001 = q;
                if (data.hasCombined()) combinedAt001 = valueAt(row, data.combinedIdx);
                if (data.hasPaired())   pairedAt001   = valueAt(row, data.pairedIdx);
                if (data.hasLower())    lowerAt001    = valueAt(row, data.lowerIdx);
            }
        }
        List<AnnotationLine> lines = new ArrayList<>();
        if (discoveries == 0) {
            return lines;
        }
        lines.add(new AnnotationLine("Total discoveries: " + discoveries,
                AnnotationLine.TOTAL));
        if (data.hasCombined() && !Double.isNaN(combinedAt001)) {
            lines.add(new AnnotationLine("Combined method: " + formatPct(combinedAt001),
                    AnnotationLine.COMBINED));
        }
        if (data.hasPaired() && !Double.isNaN(pairedAt001)) {
            lines.add(new AnnotationLine("Paired method: " + formatPct(pairedAt001),
                    AnnotationLine.PAIRED));
        }
        if (data.hasLower() && !Double.isNaN(lowerAt001)) {
            lines.add(new AnnotationLine("Lower bound: " + formatPct(lowerAt001),
                    AnnotationLine.LOWER));
        }
        return lines;
    }

    private boolean lineVisible(AnnotationLine ln) {
        switch (ln.method) {
            case AnnotationLine.COMBINED: return showCombined;
            case AnnotationLine.PAIRED:   return showPaired;
            case AnnotationLine.LOWER:    return showLower;
            default:                      return true;
        }
    }

    /**
     * True when exactly one upper-bound method is currently rendered on the
     * chart — either the data only contains one of them, or both are present
     * but the user has unchecked the other in the controls panel.
     */
    private boolean exactlyOneUpperBoundVisible() {
        boolean combinedVisible = dataHasCombined && showCombined;
        boolean pairedVisible   = dataHasPaired   && showPaired;
        return combinedVisible ^ pairedVisible;
    }

    /**
     * Install (or re-install) a legend label generator that swaps the visible
     * upper-bound series' name to "Upper bound" automatically whenever exactly
     * one upper-bound series is on the chart. Re-setting the generator fires
     * a RendererChangeEvent so the legend re-renders with fresh labels.
     */
    private void refreshLegendLabels() {
        if (chart == null) return;
        XYPlot plot = chart.getXYPlot();
        XYItemRenderer renderer = plot.getRenderer();
        renderer.setLegendItemLabelGenerator(
                new org.jfree.chart.labels.StandardXYSeriesLabelGenerator() {
                    @Override
                    public String generateLabel(org.jfree.data.xy.XYDataset ds, int series) {
                        String key = String.valueOf(ds.getSeriesKey(series));
                        if (exactlyOneUpperBoundVisible()
                                && ("Combined method".equals(key)
                                        || "Paired method".equals(key))) {
                            return "Upper bound";
                        }
                        return key;
                    }
                });
    }

    /**
     * Remove any previously-placed discovery annotations and re-add them at
     * positions derived from the current axis bounds, so the block always
     * lands near the top-left of the visible plot area. Lines whose method
     * is toggled off in the controls panel are skipped.
     */
    private void positionDiscoveryAnnotation() {
        if (chart == null || annotationLines == null || annotationLines.isEmpty()) {
            return;
        }
        XYPlot plot = chart.getXYPlot();
        if (textAnnotations != null) {
            for (XYTextAnnotation a : textAnnotations) {
                plot.removeAnnotation(a);
            }
        }
        textAnnotations = new ArrayList<>();
        double xLo = plot.getDomainAxis().getLowerBound();
        double xHi = plot.getDomainAxis().getUpperBound();
        double yLo = plot.getRangeAxis().getLowerBound();
        double yHi = plot.getRangeAxis().getUpperBound();
        double xSpan = xHi - xLo;
        double ySpan = yHi - yLo;
        // X anchor: park just past the FDR=1% guide (0.0105) when that point
        // sits in the LEFT HALF of the visible x-range — that leaves enough
        // room for the text to render without clipping on the right edge.
        // Otherwise (1% guide near the right side, or 1% out of view entirely)
        // fall back to 5% in from the visible left edge.
        boolean useGuideAnchor = xLo <= 0.0105 && 0.0105 <= xHi
                && (0.0105 - xLo) <= 0.5 * xSpan;
        double x = useGuideAnchor ? 0.0105 : xLo + 0.05 * xSpan;
        // Y anchor scales with the visible y-span so the block hugs the
        // top of the visible region regardless of where it starts.
        double y = yLo + 0.92 * ySpan;
        double lineStep = 0.06 * ySpan;
        // JFreeChart can't render multi-line in a single XYTextAnnotation, so
        // stack one annotation per visible line. Hidden lines are skipped so
        // the block stays compact when methods are toggled off.
        int row = 0;
        boolean relabel = exactlyOneUpperBoundVisible();
        for (AnnotationLine ln : annotationLines) {
            if (!lineVisible(ln)) continue;
            String displayText = ln.text;
            if (relabel && AnnotationLine.COMBINED.equals(ln.method)) {
                displayText = displayText.replaceFirst("Combined method:", "Upper bound:");
            } else if (relabel && AnnotationLine.PAIRED.equals(ln.method)) {
                displayText = displayText.replaceFirst("Paired method:", "Upper bound:");
            }
            XYTextAnnotation ann = new XYTextAnnotation(
                    displayText, x, y - row * lineStep);
            ann.setTextAnchor(org.jfree.chart.ui.TextAnchor.TOP_LEFT);
            ann.setFont(ANNOTATION_FONT);
            ann.setPaint(Color.BLACK);
            plot.addAnnotation(ann);
            textAnnotations.add(ann);
            row++;
        }
    }

    private static String formatPct(double v) {
        return String.format(Locale.US, "%.2f%%", v * 100.0);
    }
}
