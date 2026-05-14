package main.java.gui;

import java.io.IOException;
import java.io.InputStream;
import java.net.URI;
import java.net.http.HttpClient;
import java.net.http.HttpRequest;
import java.net.http.HttpResponse;
import java.time.Duration;
import java.util.Optional;
import java.util.Properties;
import java.util.prefs.Preferences;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Checks GitHub for newer FDRBench releases and remembers the user's
 * preferences (auto-check on/off, last-check timestamp, skipped version).
 *
 * <p>Uses only the JDK 11+ HttpClient and a small regex parser so we don't
 * need to pull in a JSON dependency just for a single endpoint that returns
 * a handful of fields.
 */
public final class UpdateChecker {

    private static final String RELEASES_API_URL =
            "https://api.github.com/repos/Noble-Lab/FDRBench/releases/latest";
    private static final Duration CONNECT_TIMEOUT = Duration.ofSeconds(5);
    private static final Duration REQUEST_TIMEOUT = Duration.ofSeconds(10);

    /** Skip the startup auto-check if a successful check ran within this window. */
    private static final long AUTO_CHECK_THROTTLE_MS = 12L * 60L * 60L * 1000L; // 12h

    // Preferences keys (java.util.prefs.Preferences). All under the same node so
    // the user can wipe them at once if they ever export/import preferences.
    private static final String PREFS_NODE         = "main/java/gui/UpdateChecker";
    private static final String KEY_AUTO_CHECK     = "autoCheck";
    private static final String KEY_SKIPPED_TAG    = "skippedTag";
    private static final String KEY_LAST_CHECK_MS  = "lastCheckMs";

    private UpdateChecker() {}

    // ---------------------------------------------------------------- versions

    /**
     * Current FDRBench version, read from the filtered {@code project.properties}
     * resource (Maven substitutes {@code ${project.version}} at build time).
     * Falls back to {@code "0.0.0"} if the resource can't be loaded.
     */
    public static String currentVersion() {
        try (InputStream in = UpdateChecker.class
                .getResourceAsStream("/project.properties")) {
            if (in == null) return "0.0.0";
            Properties p = new Properties();
            p.load(in);
            String v = p.getProperty("version", "0.0.0");
            return v.trim();
        } catch (IOException e) {
            return "0.0.0";
        }
    }

    /**
     * Compare two semver-ish strings ("1.2.3", "v1.2.3", "1.2"). Returns
     * {@code true} iff {@code latest} is strictly newer than {@code current}.
     * Non-numeric trailing segments (e.g. {@code "-beta"}) are stripped and
     * treated as equal numerics — good enough for FDRBench's plain X.Y.Z line.
     */
    public static boolean isNewer(String latest, String current) {
        int[] a = parseSemver(latest);
        int[] b = parseSemver(current);
        for (int i = 0; i < Math.max(a.length, b.length); i++) {
            int ai = i < a.length ? a[i] : 0;
            int bi = i < b.length ? b[i] : 0;
            if (ai != bi) return ai > bi;
        }
        return false;
    }

    private static int[] parseSemver(String v) {
        if (v == null) return new int[0];
        String s = v.trim();
        if (s.startsWith("v") || s.startsWith("V")) s = s.substring(1);
        // Drop any pre-release / build suffix so we only compare numerics.
        int dash = s.indexOf('-');
        if (dash >= 0) s = s.substring(0, dash);
        String[] parts = s.split("\\.");
        int[] out = new int[parts.length];
        for (int i = 0; i < parts.length; i++) {
            try {
                out[i] = Integer.parseInt(parts[i].replaceAll("\\D.*$", ""));
            } catch (NumberFormatException e) {
                out[i] = 0;
            }
        }
        return out;
    }

    // ---------------------------------------------------------------- HTTP

    public static final class ReleaseInfo {
        public final String tagName;
        public final String name;
        public final String htmlUrl;
        public ReleaseInfo(String tagName, String name, String htmlUrl) {
            this.tagName = tagName;
            this.name    = name;
            this.htmlUrl = htmlUrl;
        }
    }

    /**
     * Fetch the latest release from GitHub. Throws on network / HTTP errors;
     * returns {@code Optional.empty()} only if the response is missing
     * required fields.
     */
    public static Optional<ReleaseInfo> fetchLatestRelease()
            throws IOException, InterruptedException {
        HttpClient client = HttpClient.newBuilder()
                .connectTimeout(CONNECT_TIMEOUT)
                .followRedirects(HttpClient.Redirect.NORMAL)
                .build();
        HttpRequest req = HttpRequest.newBuilder()
                .uri(URI.create(RELEASES_API_URL))
                .timeout(REQUEST_TIMEOUT)
                .header("Accept", "application/vnd.github+json")
                .header("User-Agent", "FDRBench/" + currentVersion())
                .GET()
                .build();
        HttpResponse<String> resp = client.send(req,
                HttpResponse.BodyHandlers.ofString());
        if (resp.statusCode() != 200) {
            throw new IOException("GitHub responded HTTP " + resp.statusCode());
        }
        String body = resp.body();
        String tag  = extractStringField(body, "tag_name");
        String name = extractStringField(body, "name");
        String url  = extractStringField(body, "html_url");
        if (tag == null || tag.isEmpty()) {
            return Optional.empty();
        }
        return Optional.of(new ReleaseInfo(tag, name == null ? tag : name,
                url == null ? "https://github.com/Noble-Lab/FDRBench/releases/latest" : url));
    }

    /**
     * Extract the first occurrence of {@code "key": "value"} from a JSON
     * payload. Sufficient for the few flat top-level fields we need from
     * GitHub's releases response; intentionally not a general JSON parser.
     */
    private static String extractStringField(String json, String key) {
        Pattern pat = Pattern.compile(
                "\"" + Pattern.quote(key) + "\"\\s*:\\s*\"((?:\\\\.|[^\"\\\\])*)\"");
        Matcher m = pat.matcher(json);
        if (!m.find()) return null;
        // Unescape the minimal set of JSON escapes we might encounter.
        return m.group(1)
                .replace("\\\"", "\"")
                .replace("\\\\", "\\")
                .replace("\\/", "/")
                .replace("\\n", "\n")
                .replace("\\r", "\r")
                .replace("\\t", "\t");
    }

    // ---------------------------------------------------------------- prefs

    private static Preferences prefs() {
        return Preferences.userRoot().node(PREFS_NODE);
    }

    public static boolean isAutoCheckEnabled() {
        return prefs().getBoolean(KEY_AUTO_CHECK, true);
    }

    public static void setAutoCheckEnabled(boolean enabled) {
        prefs().putBoolean(KEY_AUTO_CHECK, enabled);
    }

    public static boolean isVersionSkipped(String tag) {
        if (tag == null) return false;
        return tag.equals(prefs().get(KEY_SKIPPED_TAG, ""));
    }

    public static void skipVersion(String tag) {
        if (tag == null) return;
        prefs().put(KEY_SKIPPED_TAG, tag);
    }

    public static long lastCheckMillis() {
        return prefs().getLong(KEY_LAST_CHECK_MS, 0L);
    }

    public static void recordCheckAttempt() {
        prefs().putLong(KEY_LAST_CHECK_MS, System.currentTimeMillis());
    }

    /**
     * True when an auto-check should run on startup: the user hasn't disabled
     * it AND the throttle window (12h) has elapsed since the last attempt.
     */
    public static boolean shouldAutoCheckOnStartup() {
        if (!isAutoCheckEnabled()) return false;
        long elapsed = System.currentTimeMillis() - lastCheckMillis();
        return elapsed >= AUTO_CHECK_THROTTLE_MS;
    }
}
