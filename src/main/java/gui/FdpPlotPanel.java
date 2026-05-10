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
    private JFreeChart chart;
    private File loadedCsv;

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
            ChartPanel newPanel = new ChartPanel(chart, true, true, true, true, true);
            newPanel.setPreferredSize(new Dimension(640, 480));
            newPanel.setMouseWheelEnabled(true);
            // Double-click anywhere on the chart to restore the original
            // auto-range on both axes — quick reset after the user drags or
            // scroll-zooms into a region of the plot.
            newPanel.addMouseListener(new java.awt.event.MouseAdapter() {
                @Override
                public void mouseClicked(java.awt.event.MouseEvent e) {
                    if (e.getClickCount() == 2
                            && SwingUtilities.isLeftMouseButton(e)) {
                        newPanel.restoreAutoBounds();
                    }
                }
            });

            if (chartPanel != null) {
                remove(chartPanel);
            } else {
                remove(placeholder);
            }
            chartPanel = newPanel;
            add(chartPanel, BorderLayout.CENTER);
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
        NumberFormat pct = NumberFormat.getPercentInstance(Locale.US);
        pct.setMaximumFractionDigits(0);
        xAxis.setNumberFormatOverride(pct);
        yAxis.setNumberFormatOverride(pct);
        xAxis.setRange(0, max);
        yAxis.setRange(0, max);
        xAxis.setTickUnit(new NumberTickUnit(niceTick(max)));
        yAxis.setTickUnit(new NumberTickUnit(niceTick(max)));
        xAxis.setLabelFont(AXIS_LABEL_FONT);
        yAxis.setLabelFont(AXIS_LABEL_FONT);
        xAxis.setTickLabelFont(TICK_FONT);
        yAxis.setTickLabelFont(TICK_FONT);

        // Diagonal y = x reference (gray) and vertical FDR=1% guide (blue, dashed).
        plot.addAnnotation(new XYLineAnnotation(0, 0, max, max,
                new BasicStroke(1.0f), COLOR_REFERENCE));
        ValueMarker fdrMarker = new ValueMarker(0.01);
        fdrMarker.setPaint(COLOR_FDR_LINE);
        fdrMarker.setStroke(new BasicStroke(1.0f, BasicStroke.CAP_BUTT,
                BasicStroke.JOIN_MITER, 1.0f, new float[]{4f, 4f}, 0f));
        plot.addDomainMarker(fdrMarker);

        // Annotation block: total discoveries at q ≤ 0.01 plus per-method FDP.
        addDiscoveryAnnotation(plot, data, max);

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
        if (range <= 0.02) return 0.005;
        if (range <= 0.05) return 0.01;
        if (range <= 0.10) return 0.02;
        if (range <= 0.25) return 0.05;
        if (range <= 0.50) return 0.10;
        return 0.20;
    }

    private static void addDiscoveryAnnotation(XYPlot plot, FdpData data, double max) {
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
        if (discoveries == 0) {
            return;
        }
        StringBuilder text = new StringBuilder();
        text.append("Total discoveries: ").append(discoveries);
        if (data.hasCombined() && !Double.isNaN(combinedAt001)) {
            text.append("\nCombined method: ").append(formatPct(combinedAt001));
        }
        if (data.hasPaired() && !Double.isNaN(pairedAt001)) {
            text.append("\nPaired method: ").append(formatPct(pairedAt001));
        }
        if (data.hasLower() && !Double.isNaN(lowerAt001)) {
            text.append("\nLower bound: ").append(formatPct(lowerAt001));
        }

        double x = (Math.abs(max - 0.01) <= 0.02) ? max * 0.10 : 0.0105;
        double y = max * 0.92;
        // JFreeChart can't render multi-line in a single XYTextAnnotation, so
        // stack one annotation per line.
        String[] lines = text.toString().split("\n");
        double lineStep = max * 0.06;
        for (int i = 0; i < lines.length; i++) {
            XYTextAnnotation ann = new XYTextAnnotation(lines[i], x, y - i * lineStep);
            ann.setTextAnchor(org.jfree.chart.ui.TextAnchor.TOP_LEFT);
            ann.setFont(ANNOTATION_FONT);
            ann.setPaint(Color.BLACK);
            plot.addAnnotation(ann);
        }
    }

    private static String formatPct(double v) {
        return String.format(Locale.US, "%.2f%%", v * 100.0);
    }
}
