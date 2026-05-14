package main.java.gui;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Cursor;
import java.awt.Desktop;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.Insets;
import java.net.URI;

/**
 * Slim notification bar that sits above the main header. Hidden by default;
 * {@link #show(UpdateChecker.ReleaseInfo, String)} reveals it with text and
 * action buttons; {@link #hideBanner()} tucks it away again.
 *
 * <p>Buttons: <b>View release</b> opens the GitHub release page in the system
 * browser; <b>Skip this version</b> persists the tag so the auto-check stops
 * pestering for that release; <b>Dismiss</b> just hides the banner for the
 * current session.
 */
public class UpdateBanner extends JPanel {

    // Muted amber tone so the banner reads as informational, not alarming.
    private static final Color BG_COLOR     = new Color(0xFFF4CE);
    private static final Color BORDER_COLOR = new Color(0xE0C36A);
    private static final Color TEXT_COLOR   = new Color(0x4E3A00);

    private final JLabel messageLabel;
    private final JButton viewButton;
    private final JButton skipButton;
    private final JButton dismissButton;

    private UpdateChecker.ReleaseInfo release;

    public UpdateBanner() {
        super(new BorderLayout(8, 0));
        setOpaque(true);
        setBackground(BG_COLOR);
        setBorder(BorderFactory.createCompoundBorder(
                BorderFactory.createMatteBorder(0, 0, 1, 0, BORDER_COLOR),
                BorderFactory.createEmptyBorder(6, 12, 6, 8)));

        messageLabel = new JLabel(" ", SwingConstants.LEFT);
        messageLabel.setForeground(TEXT_COLOR);
        messageLabel.setFont(messageLabel.getFont().deriveFont(Font.PLAIN));
        add(messageLabel, BorderLayout.CENTER);

        JPanel buttons = new JPanel(new FlowLayout(FlowLayout.RIGHT, 6, 0));
        buttons.setOpaque(false);
        viewButton    = makeLinkButton("View release");
        skipButton    = makeLinkButton("Skip this version");
        dismissButton = makeLinkButton("Dismiss");
        viewButton.addActionListener(e -> openReleaseInBrowser());
        skipButton.addActionListener(e -> skipThisVersion());
        dismissButton.addActionListener(e -> hideBanner());
        buttons.add(viewButton);
        buttons.add(skipButton);
        buttons.add(dismissButton);
        add(buttons, BorderLayout.EAST);

        setVisible(false);
    }

    private static JButton makeLinkButton(String text) {
        JButton b = new JButton(text);
        b.setFocusPainted(false);
        b.setBorderPainted(false);
        b.setContentAreaFilled(false);
        b.setOpaque(false);
        b.setMargin(new Insets(2, 6, 2, 6));
        b.setForeground(new Color(0x1A5FB4));
        b.setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
        Font f = b.getFont();
        b.setFont(f.deriveFont(f.getStyle() | Font.BOLD, f.getSize2D()));
        return b;
    }

    /**
     * Show the banner advertising {@code release}. {@code currentVersion} is
     * included in the message so users can see what they're upgrading from.
     */
    public void show(UpdateChecker.ReleaseInfo release, String currentVersion) {
        if (release == null) return;
        this.release = release;
        messageLabel.setText("<html><b>FDRBench " + escape(release.tagName)
                + "</b> is available. You're running " + escape(currentVersion) + ".</html>");
        setVisible(true);
        revalidate();
        repaint();
    }

    public void hideBanner() {
        setVisible(false);
        revalidate();
        repaint();
    }

    private void openReleaseInBrowser() {
        if (release == null) return;
        try {
            if (Desktop.isDesktopSupported()
                    && Desktop.getDesktop().isSupported(Desktop.Action.BROWSE)) {
                Desktop.getDesktop().browse(URI.create(release.htmlUrl));
            } else {
                JOptionPane.showMessageDialog(this,
                        "Open this URL in your browser:\n" + release.htmlUrl,
                        "Release URL", JOptionPane.INFORMATION_MESSAGE);
            }
        } catch (Exception ex) {
            JOptionPane.showMessageDialog(this,
                    "Could not open browser: " + ex.getMessage()
                            + "\n\nURL: " + release.htmlUrl,
                    "Open release", JOptionPane.ERROR_MESSAGE);
        }
    }

    private void skipThisVersion() {
        if (release != null) {
            UpdateChecker.skipVersion(release.tagName);
        }
        hideBanner();
    }

    private static String escape(String s) {
        if (s == null) return "";
        return s.replace("&", "&amp;")
                .replace("<", "&lt;")
                .replace(">", "&gt;");
    }

    /**
     * Convenience for callers running off the EDT: marshal the show call onto
     * the EDT and skip silently when no release was found or the user opted to
     * skip this tag.
     */
    public void showOnEdt(UpdateChecker.ReleaseInfo release, String currentVersion) {
        if (release == null) return;
        if (UpdateChecker.isVersionSkipped(release.tagName)) return;
        SwingUtilities.invokeLater(() -> show(release, currentVersion));
    }
}
