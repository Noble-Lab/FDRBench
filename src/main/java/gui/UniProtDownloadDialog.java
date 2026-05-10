package main.java.gui;

import javax.swing.*;
import javax.swing.border.TitledBorder;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.io.*;
import java.net.HttpURLConnection;
import java.net.URL;
import java.nio.charset.StandardCharsets;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;

/**
 * Dialog for downloading protein databases from UniProt.
 * Adapted from Carafe's UniProtDownloadDialog to fit FDRBench's
 * "output is a file path" convention.
 */
public class UniProtDownloadDialog extends JDialog {

    // Organism options with UniProt proteome IDs
    private static final String[][] ORGANISMS = {
            { "Homo sapiens (Human)", "UP000005640" },
            { "Mus musculus (Mouse)", "UP000000589" },
            { "Saccharomyces cerevisiae (Yeast)", "UP000002311" },
            { "Escherichia coli (K-12)", "UP000000625" },
            { "Arabidopsis thaliana", "UP000006548" }
    };

    private static final String CONTAMINANTS_RESOURCE = "/0602_Universal_Contaminants_full_tags.fasta";

    private final JTextField targetField;
    private final String defaultDownloadDir;
    private final boolean contaminantsBundled;

    private ButtonGroup organismGroup;
    private JRadioButton[] organismButtons;
    private JRadioButton otherButton;
    private JTextField otherProteomeField;

    private JCheckBox reviewedCheckbox;
    private JCheckBox isoformsCheckbox;
    private JCheckBox contaminantsCheckbox;

    private JTextField spikeInField;
    private JLabel spikeInLabel;
    private JButton spikeInBrowseButton;
    private JTextField downloadDirField;
    private JButton downloadDirBrowseButton;
    private JButton downloadButton;
    private JButton cancelButton;
    private JProgressBar progressBar;
    private JLabel statusLabel;

    private SwingWorker<Integer, String> downloadWorker;

    /**
     * @param owner              parent frame
     * @param targetField        the protein-database text field that will receive the
     *                           downloaded file path on success
     * @param defaultDownloadDir initial download directory (may be empty)
     */
    public UniProtDownloadDialog(Frame owner, JTextField targetField, String defaultDownloadDir) {
        super(owner, "Download Protein Database", true);
        this.targetField = targetField;
        this.defaultDownloadDir = defaultDownloadDir == null ? "" : defaultDownloadDir;
        this.contaminantsBundled = getClass().getResource(CONTAMINANTS_RESOURCE) != null;
        getRootPane().putClientProperty("JRootPane.titleBarShowIcon", false);
        initComponents();
        pack();
        setLocationRelativeTo(owner);
    }

    private void initComponents() {
        setLayout(new BorderLayout(10, 10));
        ((JPanel) getContentPane()).setBorder(BorderFactory.createEmptyBorder(15, 15, 15, 15));

        JPanel contentPanel = new JPanel(new GridBagLayout());
        GridBagConstraints gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.weightx = 1.0;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 0, 10, 0);

        // Organism selection panel
        JPanel organismPanel = new JPanel();
        organismPanel.setLayout(new BoxLayout(organismPanel, BoxLayout.Y_AXIS));
        organismPanel.setBorder(createSectionBorder("Select organism / Input proteome ID"));

        organismGroup = new ButtonGroup();
        organismButtons = new JRadioButton[ORGANISMS.length];

        for (int i = 0; i < ORGANISMS.length; i++) {
            organismButtons[i] = new JRadioButton(ORGANISMS[i][0] + " - " + ORGANISMS[i][1]);
            organismButtons[i].setActionCommand(ORGANISMS[i][1]);
            organismButtons[i].setAlignmentX(Component.LEFT_ALIGNMENT);
            organismGroup.add(organismButtons[i]);
            organismPanel.add(organismButtons[i]);
            if (i == 0)
                organismButtons[i].setSelected(true);
        }

        // "Other" option with proteome-ID text field
        JPanel otherPanel = new JPanel(new FlowLayout(FlowLayout.LEFT, 0, 0));
        otherPanel.setAlignmentX(Component.LEFT_ALIGNMENT);
        otherButton = new JRadioButton("Other:");
        otherButton.setActionCommand("OTHER");
        organismGroup.add(otherButton);
        otherProteomeField = new JTextField(20);
        otherProteomeField.setToolTipText("Enter UniProt proteome ID (e.g., UP000000XXX)");
        otherProteomeField.putClientProperty("JTextField.placeholderText", "e.g., UP000005640");
        otherProteomeField.setEnabled(false);
        otherButton.addActionListener(e -> otherProteomeField.setEnabled(otherButton.isSelected()));
        for (JRadioButton btn : organismButtons) {
            btn.addActionListener(e -> otherProteomeField.setEnabled(false));
        }
        otherPanel.add(otherButton);
        otherPanel.add(otherProteomeField);
        organismPanel.add(otherPanel);

        contentPanel.add(organismPanel, gbc);

        // Options panel
        gbc.gridy++;
        JPanel optionsPanel = new JPanel();
        optionsPanel.setLayout(new BoxLayout(optionsPanel, BoxLayout.Y_AXIS));
        optionsPanel.setBorder(createSectionBorder("Options"));

        reviewedCheckbox = new JCheckBox("Reviewed sequences only (Swiss-Prot)", true);
        reviewedCheckbox.setToolTipText("Download only manually reviewed entries from Swiss-Prot");

        isoformsCheckbox = new JCheckBox("Include isoforms", false);
        isoformsCheckbox.setToolTipText("Include alternative protein isoforms");

        contaminantsCheckbox = new JCheckBox("Add common contaminants (cRAP)", contaminantsBundled);
        if (contaminantsBundled) {
            contaminantsCheckbox.setToolTipText("Append common laboratory contaminants to the database");
        } else {
            contaminantsCheckbox.setEnabled(false);
            contaminantsCheckbox.setToolTipText(
                    "Bundled contaminants FASTA is not available in this build of FDRBench");
        }

        optionsPanel.add(reviewedCheckbox);
        optionsPanel.add(isoformsCheckbox);
        optionsPanel.add(contaminantsCheckbox);

        contentPanel.add(optionsPanel, gbc);

        // Spike-in sequences panel
        gbc.gridy++;
        JPanel spikeInPanel = new JPanel(new BorderLayout(5, 0));
        spikeInPanel.setBorder(createSectionBorder("Spike-in sequences"));

        JPanel spikeInRow = new JPanel(new BorderLayout(5, 0));
        spikeInLabel = new JLabel("FASTA file path");
        spikeInField = new JTextField();
        spikeInField.setToolTipText("Optional: Path to spike-in FASTA file (e.g., iRT standards)");
        spikeInField.putClientProperty("JTextField.placeholderText", "Optional");
        spikeInBrowseButton = new JButton("Browse");
        spikeInBrowseButton.addActionListener(e -> {
            JFileChooser chooser = new JFileChooser();
            chooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
            chooser.setFileFilter(new javax.swing.filechooser.FileNameExtensionFilter(
                    "FASTA Files", "fasta", "fa"));
            if (!spikeInField.getText().isEmpty()) {
                File parent = new File(spikeInField.getText()).getParentFile();
                if (parent != null) {
                    chooser.setCurrentDirectory(parent);
                }
            } else if (!downloadDirField.getText().isEmpty()) {
                chooser.setCurrentDirectory(new File(downloadDirField.getText()));
            }
            if (chooser.showOpenDialog(this) == JFileChooser.APPROVE_OPTION) {
                spikeInField.setText(chooser.getSelectedFile().getAbsolutePath());
            }
        });
        spikeInRow.add(spikeInLabel, BorderLayout.WEST);
        spikeInRow.add(spikeInField, BorderLayout.CENTER);
        spikeInRow.add(spikeInBrowseButton, BorderLayout.EAST);
        spikeInPanel.add(spikeInRow, BorderLayout.CENTER);

        contentPanel.add(spikeInPanel, gbc);

        // Download directory panel
        gbc.gridy++;
        JPanel dirPanel = new JPanel(new BorderLayout(5, 0));
        dirPanel.setBorder(createSectionBorder("Download to"));

        downloadDirField = new JTextField(defaultDownloadDir);

        downloadDirBrowseButton = new JButton("Browse");
        downloadDirBrowseButton.setToolTipText(
                "Choose the directory where the downloaded FASTA file will be saved.");
        downloadDirBrowseButton.addActionListener(e -> {
            JFileChooser chooser = new JFileChooser();
            chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
            if (!downloadDirField.getText().isEmpty()) {
                File current = new File(downloadDirField.getText());
                if (current.exists()) {
                    chooser.setCurrentDirectory(current);
                } else if (current.getParentFile() != null && current.getParentFile().exists()) {
                    chooser.setCurrentDirectory(current.getParentFile());
                }
            }
            if (chooser.showOpenDialog(this) == JFileChooser.APPROVE_OPTION) {
                downloadDirField.setText(chooser.getSelectedFile().getAbsolutePath());
            }
        });

        dirPanel.add(downloadDirField, BorderLayout.CENTER);
        dirPanel.add(downloadDirBrowseButton, BorderLayout.EAST);

        contentPanel.add(dirPanel, gbc);

        // Progress panel
        gbc.gridy++;
        gbc.insets = new Insets(0, 0, 0, 0);
        JPanel progressPanel = new JPanel(new BorderLayout(5, 5));
        progressBar = new JProgressBar();
        progressBar.setIndeterminate(false);
        progressBar.setStringPainted(true);
        progressBar.setString("Ready");
        progressBar.setBorder(BorderFactory.createEmptyBorder(6, 0, 6, 0));
        statusLabel = new JLabel(" ");
        progressPanel.add(progressBar, BorderLayout.CENTER);
        progressPanel.add(statusLabel, BorderLayout.SOUTH);

        contentPanel.add(progressPanel, gbc);

        add(contentPanel, BorderLayout.CENTER);

        // Button panel
        JPanel buttonPanel = new JPanel(new FlowLayout(FlowLayout.CENTER, 10, 0));
        downloadButton = new JButton("Download");
        downloadButton.putClientProperty("JButton.buttonType", "default");
        downloadButton.addActionListener(this::onDownload);
        cancelButton = new JButton("Cancel");
        cancelButton.addActionListener(e -> onCancel());
        buttonPanel.add(cancelButton);
        buttonPanel.add(downloadButton);

        add(buttonPanel, BorderLayout.SOUTH);

        getRootPane().setDefaultButton(downloadButton);

        setMinimumSize(new Dimension(500, 420));
    }

    private void onDownload(ActionEvent e) {
        String downloadDir = downloadDirField.getText().trim();
        if (downloadDir.isEmpty()) {
            JOptionPane.showMessageDialog(this,
                    "Please specify a download directory.",
                    "Download Directory Required", JOptionPane.WARNING_MESSAGE);
            return;
        }

        File dir = new File(downloadDir);
        if (!dir.exists()) {
            int result = JOptionPane.showConfirmDialog(this,
                    "Directory does not exist. Create it?",
                    "Create Directory", JOptionPane.YES_NO_OPTION);
            if (result != JOptionPane.YES_OPTION) {
                return;
            }
            if (!dir.mkdirs()) {
                JOptionPane.showMessageDialog(this,
                        "Failed to create directory.",
                        "Error", JOptionPane.ERROR_MESSAGE);
                return;
            }
        }

        String proteomeId;
        String organismName;
        if (otherButton.isSelected()) {
            proteomeId = otherProteomeField.getText().trim();
            if (proteomeId.isEmpty()) {
                JOptionPane.showMessageDialog(this,
                        "Please enter a proteome ID.",
                        "Proteome ID Required", JOptionPane.WARNING_MESSAGE);
                return;
            }
            organismName = "custom";
        } else {
            proteomeId = organismGroup.getSelection().getActionCommand();
            organismName = "unknown";
            for (String[] entry : ORGANISMS) {
                if (entry[1].equals(proteomeId)) {
                    String[] parts = entry[0].split(" \\(");
                    if (parts.length > 1) {
                        organismName = parts[1].replace(")", "").toLowerCase();
                    } else {
                        organismName = parts[0].toLowerCase();
                    }
                    break;
                }
            }
        }

        String url = buildUniProtUrl(proteomeId);
        String filename = generateFilename(organismName);
        File outputFile = new File(dir, filename);

        setControlsEnabled(false);
        progressBar.setIndeterminate(true);
        progressBar.setString("Connecting...");
        statusLabel.setText("Connecting to UniProt...");

        downloadWorker = new SwingWorker<Integer, String>() {
            @Override
            protected Integer doInBackground() throws Exception {
                publish("Downloading from UniProt...");

                URL urlObj = new URL(url);
                HttpURLConnection conn = (HttpURLConnection) urlObj.openConnection();
                conn.setRequestProperty("Accept", "text/plain");
                conn.setConnectTimeout(30000);
                conn.setReadTimeout(300000);

                int responseCode = conn.getResponseCode();
                if (responseCode != 200) {
                    throw new IOException("UniProt returned error: " + responseCode);
                }

                int sequenceCount = 0;
                try (BufferedReader reader = new BufferedReader(
                        new InputStreamReader(conn.getInputStream(), StandardCharsets.UTF_8));
                     BufferedWriter writer = new BufferedWriter(
                             new OutputStreamWriter(new FileOutputStream(outputFile), StandardCharsets.UTF_8))) {

                    String line;
                    while ((line = reader.readLine()) != null) {
                        if (isCancelled())
                            break;
                        writer.write(line);
                        writer.newLine();
                        if (line.startsWith(">")) {
                            sequenceCount++;
                            if (sequenceCount % 1000 == 0) {
                                publish("Downloaded " + sequenceCount + " sequences...");
                            }
                        }
                    }
                    publish("Downloaded " + sequenceCount + " sequences from UniProt.");

                    if (contaminantsCheckbox.isSelected() && contaminantsBundled) {
                        publish("Appending contaminants...");
                        appendContaminants(writer);
                    }

                    String spikeInPath = spikeInField.getText().trim();
                    if (!spikeInPath.isEmpty()) {
                        File spikeInFile = new File(spikeInPath);
                        if (spikeInFile.exists()) {
                            publish("Appending spike-in sequences...");
                            appendSpikeIn(writer, spikeInFile);
                        }
                    }
                }

                return sequenceCount;
            }

            @Override
            protected void process(java.util.List<String> chunks) {
                if (!chunks.isEmpty()) {
                    statusLabel.setText(chunks.get(chunks.size() - 1));
                }
            }

            @Override
            protected void done() {
                progressBar.setIndeterminate(false);
                setControlsEnabled(true);

                if (isCancelled()) {
                    statusLabel.setText("Download cancelled.");
                    progressBar.setString("Cancelled");
                    if (outputFile.exists()) {
                        outputFile.delete();
                    }
                    return;
                }

                try {
                    Integer count = get();
                    statusLabel.setText("Download complete! (" + count + " proteins)");
                    progressBar.setString("Complete (" + count + " proteins)");
                    progressBar.setValue(100);

                    targetField.setText(outputFile.getAbsolutePath());

                    JOptionPane.showMessageDialog(UniProtDownloadDialog.this,
                            "Downloaded " + count + " proteins to:\n" + outputFile.getAbsolutePath(),
                            "Download Complete", JOptionPane.INFORMATION_MESSAGE);

                    dispose();
                } catch (Exception ex) {
                    statusLabel.setText("Download failed.");
                    progressBar.setString("Failed");
                    JOptionPane.showMessageDialog(UniProtDownloadDialog.this,
                            "Download failed: " + ex.getMessage(),
                            "Error", JOptionPane.ERROR_MESSAGE);
                    if (outputFile.exists()) {
                        outputFile.delete();
                    }
                }
            }
        };

        downloadWorker.execute();
    }

    private String buildUniProtUrl(String proteomeId) {
        StringBuilder url = new StringBuilder();
        url.append("https://rest.uniprot.org/uniprotkb/stream?");
        url.append("query=proteome:").append(proteomeId);

        if (reviewedCheckbox.isSelected()) {
            url.append("+AND+reviewed:true");
        }

        url.append("&format=fasta");

        if (isoformsCheckbox.isSelected()) {
            url.append("&includeIsoform=true");
        }

        return url.toString();
    }

    private String generateFilename(String organismName) {
        StringBuilder name = new StringBuilder();
        name.append(organismName.replaceAll("[^a-zA-Z0-9]", "_"));
        name.append("_");
        name.append(reviewedCheckbox.isSelected() ? "reviewed" : "all");
        name.append("_");

        LocalDateTime now = LocalDateTime.now();
        DateTimeFormatter formatter = DateTimeFormatter.ofPattern("yyyy-MM-dd_HH-mm-ss");
        name.append(now.format(formatter));

        if (contaminantsCheckbox.isSelected()) {
            name.append("_contaminants");
        }

        if (isoformsCheckbox.isSelected()) {
            name.append("_isoforms");
        }

        name.append(".fasta");
        return name.toString();
    }

    private void appendContaminants(BufferedWriter writer) throws IOException {
        try (InputStream is = getClass().getResourceAsStream(CONTAMINANTS_RESOURCE)) {
            if (is == null) {
                return;
            }
            try (BufferedReader reader = new BufferedReader(
                    new InputStreamReader(is, StandardCharsets.UTF_8))) {
                writer.newLine();
                String line;
                while ((line = reader.readLine()) != null) {
                    writer.write(line);
                    writer.newLine();
                }
            }
        }
    }

    private void appendSpikeIn(BufferedWriter writer, File spikeInFile) throws IOException {
        try (BufferedReader reader = new BufferedReader(
                new InputStreamReader(new FileInputStream(spikeInFile), StandardCharsets.UTF_8))) {
            writer.newLine();
            writer.write("# Spike-in sequences");
            writer.newLine();
            String line;
            while ((line = reader.readLine()) != null) {
                writer.write(line);
                writer.newLine();
            }
        }
    }

    private void setControlsEnabled(boolean enabled) {
        downloadButton.setEnabled(enabled);
        for (JRadioButton btn : organismButtons) {
            btn.setEnabled(enabled);
        }
        otherButton.setEnabled(enabled);
        otherProteomeField.setEnabled(enabled && otherButton.isSelected());
        reviewedCheckbox.setEnabled(enabled);
        isoformsCheckbox.setEnabled(enabled);
        contaminantsCheckbox.setEnabled(enabled && contaminantsBundled);
        spikeInLabel.setEnabled(enabled);
        spikeInField.setEnabled(enabled);
        spikeInBrowseButton.setEnabled(enabled);
        downloadDirField.setEnabled(enabled);
        downloadDirBrowseButton.setEnabled(enabled);
    }

    private void onCancel() {
        if (downloadWorker != null && !downloadWorker.isDone()) {
            downloadWorker.cancel(true);
        }
        dispose();
    }

    public void showDialog() {
        setVisible(true);
    }

    private TitledBorder createSectionBorder(String title) {
        Color line = UIManager.getColor("Component.borderColor");
        if (line == null) {
            line = Color.LIGHT_GRAY;
        }
        return BorderFactory.createTitledBorder(
                BorderFactory.createMatteBorder(1, 0, 0, 0, line),
                title,
                TitledBorder.CENTER,
                TitledBorder.TOP);
    }
}
