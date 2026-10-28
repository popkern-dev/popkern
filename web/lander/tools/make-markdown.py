#!/usr/bin/env python3
"""Write the site's agent-readable copies: one Markdown file per page, a
copy of the scripting guide, and /llms.txt, the guide for agents that
links them.

    python3 tools/make-markdown.py <site-dir> <scripting-guide.md> [<base-url>]

Each page in PAGES gets a .md beside it (index.html -> index.md,
blog/announcing-popkern.html -> blog/announcing-popkern.md) holding its
<main> as Markdown, with every link made absolute against <base-url>
(default https://popkern.dev) so a copy fetched on its own still resolves.
The pages name their copy with <link rel="alternate" type="text/markdown">,
and a hidden block at the top of each body points to llms.txt.

The scripting guide (docs/reference/scripting-guide.md, the source popkern
live's guide.html is rendered from) is copied to /scripting-guide.md with
its repo-relative links pointed at GitHub. llms.txt is tools/llms.txt.in
filled in: the homepage tagline, the guide's boss-plate script as the
example, and the list of page copies.

Run by tools/deploy-site-routes.sh on the staged tree, never on the source:
the copies are derived from the HTML and the guide, so they cannot drift
from them. Exits non-zero if a page is missing, has no <main>, or converts
to nothing; if the guide's example is gone; or if a shell command in
llms.txt is not also on the install page (the install page is the one
that is reviewed). Standard library only.
"""

import html
import re
import sys
from html.parser import HTMLParser
from pathlib import Path
from string import Template
from urllib.parse import urljoin

# (page, section in llms.txt). Order is the order llms.txt lists them.
PAGES = [
    ("index.html", "Site"),
    ("install/index.html", "Site"),
    ("support/index.html", "Site"),
    ("security/index.html", "Site"),
    ("blog/index.html", "Blog"),
    ("blog/announcing-popkern.html", "Blog"),
]

REPO_BLOB = "https://github.com/popkern-dev/popkern/blob/main/"
# The guide's worked example, quoted into llms.txt.
EXAMPLE_SECTION = "## 3. The boss plate, line by line"
# Lines of llms.txt that must match the install page word for word.
INSTALL_COMMAND = re.compile(r"(curl |irm |claude mcp add )")

SKIP = {"svg", "button", "script", "style", "template", "nav"}
BLOCK = {"p", "div", "section", "article", "header", "figure", "figcaption",
         "ul", "ol", "li", "pre", "blockquote", "hr", "h1", "h2", "h3", "h4",
         "h5", "h6", "table", "tr"}


def md_path(page):
    return page[:-len(".html")] + ".md"


class MainToMarkdown(HTMLParser):
    """Converts the content of <main> to Markdown. Handles the elements the
    site's pages use: headings, paragraphs, lists, links, emphasis, inline
    and block code, figures with images, and rules."""

    def __init__(self, page_url):
        super().__init__(convert_charrefs=True)
        self.page_url = page_url
        self.in_main = False
        self.skip_depth = 0
        self.out = []          # finished blocks
        self.line = []         # inline text of the block being built
        self.lists = []        # stack of ["ul"|"ol", counter]
        self.pre = None        # list of chunks while inside <pre>
        self.links = []        # stack of hrefs for open <a>
        self.prefix = ""
        self.bold_closed = False

    # -- helpers --
    def flush(self):
        text = re.sub(r"[ \t\r\n]+", " ", "".join(self.line)).strip()
        self.line = []
        card = next((h for h, start in self.links if start is None), None)
        if text and card and self.prefix.startswith("#"):
            text = f"[{text}]({card})"
        if text:
            self.out.append(self.prefix + text)
        self.prefix = ""

    def emit(self, s):
        self.line.append(s)

    # -- parser hooks --
    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        if tag == "main":
            self.in_main = True
            return
        if not self.in_main:
            return
        if self.skip_depth or tag in SKIP or "code-head" in a.get("class", ""):
            self.skip_depth += tag not in {"br", "img", "hr", "use", "path"}
            return
        if self.pre is not None:
            return
        if tag in BLOCK:
            if self.links and self.links[-1][1] is not None:
                # A card: the link wraps whole blocks. Markdown links are
                # inline, so the link moves onto the headings inside it.
                href, start = self.links[-1]
                self.line[start] = ""
                self.links[-1] = (href, None)
            self.flush()
        if tag in {"h1", "h2", "h3", "h4", "h5", "h6"}:
            self.prefix = "#" * int(tag[1]) + " "
        elif tag in {"ul", "ol"}:
            self.lists.append([tag, 0])
        elif tag == "li" and self.lists:
            kind = self.lists[-1]
            kind[1] += 1
            indent = "  " * (len(self.lists) - 1)
            self.prefix = indent + ("- " if kind[0] == "ul" else f"{kind[1]}. ")
        elif tag == "pre":
            self.pre = []
        elif tag == "hr":
            self.out.append("---")
        elif tag == "br":
            self.emit(" ")
        elif tag in {"b", "strong"}:
            self.emit("**")
        elif tag in {"i", "em"}:
            self.emit("*")
        elif tag == "code":
            self.emit("`")
        elif tag == "a":
            self.links.append((urljoin(self.page_url, a.get("href", "")),
                               len(self.line)))
            self.emit("[")
        elif tag == "img":
            src = urljoin(self.page_url, a.get("src", "").split("?")[0])
            self.flush()
            self.out.append(f"![{a.get('alt', '').strip()}]({src})")

    def handle_endtag(self, tag):
        if tag == "main":
            self.flush()
            self.in_main = False
            return
        if not self.in_main:
            return
        if self.skip_depth:
            if tag not in {"br", "img", "hr", "use", "path"}:
                self.skip_depth -= 1
            return
        if tag == "pre" and self.pre is not None:
            body = "".join(self.pre).strip("\n")
            self.pre = None
            self.out.append("```\n" + body + "\n```")
            return
        if self.pre is not None:
            return
        if tag in {"b", "strong"}:
            self.emit("**")
            # A bolded lead-in runs straight into its sentence in the markup
            # (the CSS separates them); Markdown needs the space.
            self.bold_closed = True
        elif tag in {"i", "em"}:
            self.emit("*")
        elif tag == "code":
            self.emit("`")
        elif tag == "a" and self.links:
            href, start = self.links.pop()
            if start is None:
                self.flush()
                return
            label = re.sub(r"\s+", " ", "".join(self.line[start + 1:])).strip()
            self.line[start:] = [f"[{label}]({href})"]
        elif tag in {"ul", "ol"}:
            self.flush()
            if self.lists:
                self.lists.pop()
        elif tag in BLOCK:
            self.flush()

    def handle_data(self, data):
        if not self.in_main or self.skip_depth:
            return
        if self.pre is not None:
            self.pre.append(data)
            return
        if self.bold_closed and data.strip():
            if data[0].isalnum():
                self.emit(" ")
            self.bold_closed = False
        self.emit(data)

    def markdown(self):
        blocks = []
        for b in self.out:
            # List items stay tight; everything else is its own paragraph.
            tight = blocks and re.match(r"\s*(- |\d+\. )", b) and \
                re.match(r"\s*(- |\d+\. )", blocks[-1].split("\n")[-1])
            if tight:
                blocks[-1] += "\n" + b
            else:
                blocks.append(b)
        text = "\n\n".join(blocks)
        return re.sub(r"\*\*\s*\*\*", "", text).strip() + "\n"


def head_meta(source):
    title = re.search(r"<title>(.*?)</title>", source, re.S)
    desc = re.search(r'<meta name="description" content="([^"]*)"', source)
    if not title or not desc:
        return None, None
    return html.unescape(title.group(1).strip()), html.unescape(desc.group(1))


def main(argv):
    if len(argv) not in (3, 4):
        print(__doc__, file=sys.stderr)
        return 2
    site, guide_path = Path(argv[1]), Path(argv[2])
    base = (argv[3] if len(argv) == 4 else "https://popkern.dev").rstrip("/") + "/"
    sections = {}
    for page, section in PAGES:
        path = site / page
        if not path.is_file():
            print(f"make-markdown: {page} is missing", file=sys.stderr)
            return 1
        source = path.read_text(encoding="utf-8")
        title, desc = head_meta(source)
        if title is None:
            print(f"make-markdown: {page} has no <title> or description",
                  file=sys.stderr)
            return 1
        conv = MainToMarkdown(urljoin(base, page))
        conv.feed(source)
        body = conv.markdown()
        if len(body.split()) < 20:
            print(f"make-markdown: {page} converted to almost nothing "
                  f"({len(body.split())} words): no <main>?", file=sys.stderr)
            return 1
        canonical = urljoin(base, page).removesuffix("index.html")
        (site / md_path(page)).write_text(
            f"<!-- Markdown copy of {canonical} -->\n\n{body}", encoding="utf-8")
        sections.setdefault(section, []).append(
            f"- [{title}]({urljoin(base, md_path(page))}): {desc}")

    _, tagline = head_meta((site / "index.html").read_text("utf-8"))
    guide = guide_path.read_text(encoding="utf-8")
    # The guide's links are relative to docs/reference/ in the repo.
    guide_copy = re.sub(
        r"\]\((?!https?:|#)([^)]+)\)",
        lambda m: "](" + urljoin(REPO_BLOB + "docs/reference/", m.group(1)) + ")",
        guide)
    (site / "scripting-guide.md").write_text(guide_copy, encoding="utf-8")

    example = re.search(re.escape(EXAMPLE_SECTION) + r"\s+```rhai\n(.*?)\n```",
                        guide, re.S)
    if not example:
        print(f"make-markdown: {guide_path} has no rhai block under "
              f"'{EXAMPLE_SECTION}'", file=sys.stderr)
        return 1
    page_list = []
    for section, entries in sections.items():
        page_list += [f"## {section} pages (Markdown copies)", "", *entries, ""]
    template = Template(
        (Path(__file__).parent / "llms.txt.in").read_text(encoding="utf-8"))
    llms = template.substitute(base=base, tagline=tagline,
                               example=example.group(1),
                               pages="\n".join(page_list).rstrip())
    install_page = html.unescape(
        (site / "install/index.html").read_text(encoding="utf-8"))
    for command in llms.splitlines():
        if INSTALL_COMMAND.match(command) and command not in install_page:
            print(f"make-markdown: llms.txt says `{command}` but the install "
                  "page does not", file=sys.stderr)
            return 1
    (site / "llms.txt").write_text(llms, encoding="utf-8")
    print(f"make-markdown: {len(PAGES)} pages, the scripting guide and "
          f"llms.txt in {site}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
