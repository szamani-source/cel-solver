#!/usr/bin/env python3
"""Check the HTML user manual in docs/manual (series M, decision 0116).

Offline checks, run by tests/test_manual.sh in CI:
  1. NAV in assets/manual.js and the .html files in docs/manual agree.
  2. Every page has the common head (charset, viewport, title, manual.css,
     manual.js), exactly one <h1> and one <main>, and balanced tags.
  3. Every local link and image resolves, including #anchors on other pages.
  4. Every citation (a.cite -> references.html#key) has an entry there, and
     every entry in references.html is cited somewhere.
  5. Every h2 has an id (the side menu builds the on-page contents from them).
  6. Equation delimiters are balanced, and, when node is installed, every
     equation is typeset by the vendored KaTeX in strict mode (no errors).
  7. Page-specific content rules (REQUIRED below): each work package adds the
     terms its page must cover.
  8. Once case-setup.html exists: every input keyword that the C++ code reads
     from a dictionary appears on it inside <code>.

Usage: python3 tools/check_manual.py [--online]
  --online also checks that every DOI in references.html resolves (needs network).
"""

import html.parser
import json
import os
import re
import shutil
import subprocess
import sys
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MAN = os.path.join(ROOT, "docs", "manual")

VOID = {"area", "base", "br", "col", "embed", "hr", "img", "input", "link", "meta",
        "source", "track", "wbr"}

# Terms each page must mention (case-insensitive, in its visible text or code).
# A package adds its page here when it writes the page.
REQUIRED = {
    "index.html": ["OpenFOAM v2512", "GPL-3.0", "Eulerian", "Lagrangian", "Mohr-Coulomb",
                   "Drucker-Prager"],
    "references.html": [],
}

errors = []


def err(msg):
    errors.append(msg)


class Page(html.parser.HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.stack = []
        self.ids = set()
        self.links = []          # (tag, attr value)
        self.cites = []          # href of a.cite
        self.h1 = 0
        self.main = 0
        self.h2_without_id = []
        self.title = ""
        self._in_title = False
        self._in_h2 = None
        self.heads = {"charset": False, "viewport": False, "css": False, "js": False}
        self.text = []
        self.code_text = []
        self._code_depth = 0
        self.problems = []

    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        if "id" in a:
            if a["id"] in self.ids:
                self.problems.append("duplicate id '%s'" % a["id"])
            self.ids.add(a["id"])
        if tag == "a" and "href" in a:
            self.links.append(a["href"])
            if "cite" in (a.get("class") or "").split():
                self.cites.append(a["href"])
        if tag in ("img", "script", "source") and "src" in a:
            self.links.append(a["src"])
        if tag == "link" and "href" in a:
            self.links.append(a["href"])
            if a["href"] == "assets/manual.css":
                self.heads["css"] = True
        if tag == "script" and a.get("src") == "assets/manual.js":
            self.heads["js"] = True
        if tag == "meta":
            if "charset" in a:
                self.heads["charset"] = True
            if a.get("name") == "viewport":
                self.heads["viewport"] = True
        if tag == "h1":
            self.h1 += 1
        if tag == "main":
            self.main += 1
        if tag == "h2":
            self._in_h2 = "id" in a
            if "id" not in a:
                self.h2_without_id.append(self.getpos()[0])
        if tag == "title":
            self._in_title = True
        if tag in ("code", "pre"):
            self._code_depth += 1
        if tag not in VOID:
            self.stack.append((tag, self.getpos()[0]))

    def handle_startendtag(self, tag, attrs):
        self.handle_starttag(tag, attrs)
        if tag not in VOID and self.stack and self.stack[-1][0] == tag:
            self.stack.pop()

    def handle_endtag(self, tag):
        if tag in VOID:
            return
        if tag == "title":
            self._in_title = False
        if tag in ("code", "pre"):
            self._code_depth = max(0, self._code_depth - 1)
        if not self.stack:
            self.problems.append("line %d: </%s> without an open tag" % (self.getpos()[0], tag))
            return
        if self.stack[-1][0] != tag:
            # Optional end tags (p, li, td, ...) close implicitly; find the match.
            names = [t for t, _ in self.stack]
            if tag in names:
                while self.stack and self.stack[-1][0] != tag:
                    t, line = self.stack.pop()
                    if t not in ("p", "li", "td", "th", "tr", "dt", "dd", "option", "thead",
                                 "tbody"):
                        self.problems.append("line %d: <%s> not closed before </%s>"
                                             % (line, t, tag))
                self.stack.pop()
            else:
                self.problems.append("line %d: </%s> does not match <%s>"
                                     % (self.getpos()[0], tag, self.stack[-1][0]))
            return
        self.stack.pop()

    def handle_data(self, data):
        if self._in_title:
            self.title += data
        self.text.append(data)
        if self._code_depth:
            self.code_text.append(data)


def read_nav():
    src = open(os.path.join(MAN, "assets", "manual.js"), encoding="utf-8").read()
    m = re.search(r"var NAV = (\[.*?\n\]);", src, re.S)
    if not m:
        err("manual.js: NAV list not found")
        return []
    pages = re.findall(r'\["([^"]+\.html)",\s*"([^"]+)"\]', m.group(1))
    return [p for p, _ in pages]


def parse_all():
    pages = {}
    for name in sorted(os.listdir(MAN)):
        if not name.endswith(".html"):
            continue
        p = Page()
        p.feed(open(os.path.join(MAN, name), encoding="utf-8").read())
        p.close()
        pages[name] = p
    return pages


def check_structure(pages, nav):
    for n in nav:
        if n not in pages:
            err("NAV lists %s but the file does not exist" % n)
    for n in pages:
        if n not in nav:
            err("%s is not in NAV (assets/manual.js)" % n)
    if len(set(nav)) != len(nav):
        err("NAV lists a page twice")
    for n, p in pages.items():
        for k, ok in p.heads.items():
            if not ok:
                err("%s: missing head item '%s'" % (n, k))
        if not p.title.strip().endswith("celFoam manual"):
            err("%s: <title> should end with 'celFoam manual'" % n)
        if p.h1 != 1:
            err("%s: %d <h1> (want 1)" % (n, p.h1))
        if p.main != 1:
            err("%s: %d <main> (want 1)" % (n, p.main))
        for line in p.h2_without_id:
            err("%s line %d: <h2> without id" % (n, line))
        for prob in p.problems:
            err("%s: %s" % (n, prob))
        left = [t for t, _ in p.stack if t not in ("p", "li", "td", "th", "tr", "html",
                                                   "body", "head")]
        if left:
            err("%s: unclosed tags at end: %s" % (n, left))


def check_links(pages):
    for n, p in pages.items():
        for href in p.links:
            if re.match(r"^(https?:|mailto:)", href):
                continue
            if href.startswith("#"):
                if href[1:] and href[1:] not in p.ids:
                    err("%s: anchor %s not found" % (n, href))
                continue
            path, _, anchor = href.partition("#")
            target = os.path.normpath(os.path.join(MAN, path))
            if not os.path.exists(target):
                err("%s: link to missing file %s" % (n, href))
                continue
            if anchor:
                tn = os.path.relpath(target, MAN)
                if tn in pages and anchor not in pages[tn].ids:
                    err("%s: %s has no anchor #%s" % (n, tn, anchor))


def check_citations(pages):
    if "references.html" not in pages:
        err("references.html missing")
        return
    refs = pages["references.html"]
    ref_keys = {i for i in refs.ids if i.startswith("ref-")}
    cited = set()
    for n, p in pages.items():
        for href in p.cites:
            path, _, key = href.partition("#")
            if path != "references.html" or not key.startswith("ref-"):
                err("%s: citation link %s should be references.html#ref-<key>" % (n, href))
                continue
            cited.add(key)
            if key not in ref_keys:
                err("%s: cites %s, which is not in references.html" % (n, key))
    for k in sorted(ref_keys - cited):
        err("references.html: %s is never cited" % k)


MATH_RE = re.compile(r"\$\$(.+?)\$\$|\\\((.+?)\\\)", re.S)


def visible_source(src):
    """Page source without <pre>, <code>, <script> and <style> blocks."""
    return re.sub(r"<(pre|code|script|style)\b.*?</\1>", " ", src, flags=re.S)


def check_math():
    formulas = []
    for name in sorted(os.listdir(MAN)):
        if not name.endswith(".html"):
            continue
        src = visible_source(open(os.path.join(MAN, name), encoding="utf-8").read())
        if src.count("$$") % 2:
            err("%s: odd number of $$" % name)
        if src.count("\\(") != src.count("\\)"):
            err("%s: unbalanced \\( \\)" % name)
        for m in MATH_RE.finditer(src):
            tex = m.group(1) if m.group(1) is not None else m.group(2)
            tex = (tex.replace("&lt;", "<").replace("&gt;", ">").replace("&amp;", "&"))
            formulas.append((name, tex, m.group(1) is not None))
    node = shutil.which("node")
    if not node:
        print("note: node not found, equations not typeset-checked")
        return len(formulas)
    js = r"""
const katex = require(process.argv[1]);
const items = JSON.parse(require('fs').readFileSync(0, 'utf8'));
let bad = 0;
for (const [page, tex, display] of items) {
  try { katex.renderToString(tex, {displayMode: display, throwOnError: true, strict: 'error'}); }
  catch (e) { bad++; console.log(page + ': ' + e.message.split('\n')[0] + '  in: ' + tex.slice(0, 80)); }
}
process.exit(bad ? 1 : 0);
"""
    r = subprocess.run([node, "-e", js, os.path.join(MAN, "assets", "katex", "katex.min.js")],
                       input=json.dumps(formulas), capture_output=True, text=True)
    if r.returncode:
        for line in r.stdout.strip().splitlines() + r.stderr.strip().splitlines()[-3:]:
            err("equation: " + line)
    return len(formulas)


def check_required(pages):
    for n, terms in REQUIRED.items():
        if n not in pages:
            err("%s required but missing" % n)
            continue
        text = " ".join(pages[n].text).lower()
        text = re.sub(r"\s+", " ", text)
        for t in terms:
            if t.lower() not in text:
                err("%s: must mention '%s'" % (n, t))


def code_keywords():
    """Dictionary keywords the C++ code reads (lookup, get, getOrDefault, ...)."""
    pat = re.compile(r"\b(?:lookupOrDefault|getOrDefault|getOrAdd|lookup|get|found|subDict"
                     r"|optionalSubDict|isDict|readIfPresent|readEntry|getCheck"
                     r"|getCheckOrDefault)\s*(?:<[^>()]*>)?\s*\(\s*\"([A-Za-z_][A-Za-z0-9_]*)\"")
    keys = set()
    for base in ("src", "applications"):
        for d, _, files in os.walk(os.path.join(ROOT, base)):
            for f in files:
                if f.endswith((".C", ".H")):
                    keys.update(pat.findall(open(os.path.join(d, f), encoding="utf-8",
                                                 errors="replace").read()))
    return keys


def check_keywords(pages):
    if "case-setup.html" not in pages:
        return 0
    code = " ".join(pages["case-setup.html"].code_text)
    words = set(re.findall(r"[A-Za-z_][A-Za-z0-9_]*", code))
    keys = code_keywords()
    for k in sorted(keys - words):
        err("case-setup.html: keyword '%s' (read by the code) is not documented" % k)
    return len(keys)


def check_dois_online():
    src = open(os.path.join(MAN, "references.html"), encoding="utf-8").read()
    dois = sorted(set(re.findall(r'href="https://doi\.org/([^"]+)"', src)))
    for d in dois:
        url = "https://doi.org/api/handles/" + d
        try:
            with urllib.request.urlopen(url, timeout=20) as r:
                data = json.load(r)
            if data.get("responseCode") != 1:
                err("DOI does not resolve: " + d)
        except Exception as e:  # noqa: BLE001
            err("DOI check failed for %s: %s" % (d, e))
    return len(dois)


def main():
    nav = read_nav()
    pages = parse_all()
    check_structure(pages, nav)
    check_links(pages)
    check_citations(pages)
    nmath = check_math()
    check_required(pages)
    nkeys = check_keywords(pages)
    ndoi = check_dois_online() if "--online" in sys.argv else 0
    print("pages %d, equations %d, keywords %d, DOIs checked %d"
          % (len(pages), nmath, nkeys, ndoi))
    if errors:
        for e in errors:
            print("FAIL " + e)
        print("%d problem(s)" % len(errors))
        return 1
    print("manual OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
