#!/usr/bin/env python3
"""Build a self-contained, standalone version of the cminpack home page.

The home page (index.html) normally references an external stylesheet
(fd.css) and a handful of images (the background stripe, the list bullets and
the W3C validity badges at the bottom of the page). This script produces a
single self-contained HTML file by:

  * inlining fd.css into a <style> element (replacing the <link> tag), and
  * embedding every referenced image as a base64 data: URI, both in the HTML
    (<img src="images/...">) and inside the CSS (url(images/...)).

The result has no external dependencies and can be dropped anywhere (archived,
emailed, or served without the images/ directory).

Usage:
    python3 build_standalone.py [INPUT_HTML] [OUTPUT_HTML]

Defaults: INPUT_HTML=index.html, OUTPUT_HTML=index-standalone.html
(paths are resolved relative to this script's directory).
"""
import base64
import os
import re
import sys

MIME = {
    ".png": "image/png",
    ".gif": "image/gif",
    ".jpg": "image/jpeg",
    ".jpeg": "image/jpeg",
    ".svg": "image/svg+xml",
}


def data_uri(path):
    """Return a base64 data: URI for the file at path."""
    ext = os.path.splitext(path)[1].lower()
    mime = MIME.get(ext, "application/octet-stream")
    with open(path, "rb") as f:
        b64 = base64.b64encode(f.read()).decode("ascii")
    return "data:%s;base64,%s" % (mime, b64)


def inline_css_images(css, base_dir):
    """Replace url(images/foo.png) in CSS with url(data:...)."""
    def repl(m):
        ref = m.group(1).strip("'\"")
        asset = os.path.join(base_dir, ref)
        if os.path.isfile(asset):
            return "url(%s)" % data_uri(asset)
        sys.stderr.write("warning: CSS asset not found: %s\n" % asset)
        return m.group(0)
    return re.sub(r"url\(([^)]+)\)", repl, css)


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    inp = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "index.html")
    out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(here, "index-standalone.html")
    base_dir = os.path.dirname(os.path.abspath(inp))

    with open(inp, encoding="utf-8") as f:
        html = f.read()

    # 1) inline the stylesheet(s) referenced by <link rel="stylesheet" href="...">
    def link_repl(m):
        href = m.group(1)
        css_path = os.path.join(base_dir, href)
        if not os.path.isfile(css_path):
            sys.stderr.write("warning: stylesheet not found: %s\n" % css_path)
            return m.group(0)
        with open(css_path, encoding="utf-8") as f:
            css = f.read()
        css = inline_css_images(css, base_dir)
        return '<style type="text/css">\n%s\n</style>' % css

    html = re.sub(
        r'<link\b[^>]*\brel="stylesheet"[^>]*\bhref="([^"]+)"[^>]*/?>',
        link_repl, html)

    # 2) embed <img src="images/..."> (and any other local image src) as data URIs
    def img_repl(m):
        pre, src, post = m.group(1), m.group(2), m.group(3)
        if re.match(r"[a-zA-Z][a-zA-Z0-9+.-]*:", src) or src.startswith("data:"):
            return m.group(0)  # already absolute or a data URI
        asset = os.path.join(base_dir, src)
        if os.path.isfile(asset):
            return '%s"%s"%s' % (pre, data_uri(asset), post)
        sys.stderr.write("warning: image not found: %s\n" % asset)
        return m.group(0)

    html = re.sub(r'(<img\b[^>]*\bsrc=)"([^"]+)"([^>]*>)', img_repl, html)

    with open(out, "w", encoding="utf-8") as f:
        f.write(html)
    sys.stderr.write("wrote %s (%d bytes)\n" % (out, len(html.encode("utf-8"))))


if __name__ == "__main__":
    main()
