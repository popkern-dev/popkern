#!/usr/bin/env python3
"""Stamp a content hash onto every local asset reference in the HTML.

Without a build step there is nothing to fingerprint filenames, so a CDN or a
browser can serve a stale stylesheet long after the file changed. Cloudflare
caches by extension and does not honour a later Cache-Control on an object it
already holds, so the only reliable cache-bust is a new URL.

Run after changing anything in assets/. Idempotent: an existing ?v= is replaced.
"""
import hashlib, pathlib, re, sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
PAGES = ["index.html", "404.html", "blog/index.html",
         "blog/announcing-popkern.html", "playground/index.html",
         "install/index.html", "support/index.html"]
ASSET_RE = re.compile(r"((?:href|src)=\")(/?(?:\.\./)*assets/[^\"?]+)(?:\?v=[0-9a-f]+)?(\")")

def digest(ref: str, page: pathlib.Path) -> str | None:
    if ref.startswith("/"):
        target = (ROOT / ref.lstrip("/")).resolve()
    else:
        target = (page.parent / ref).resolve()
    if not target.is_file():
        print(f"  !! missing {target}", file=sys.stderr)
        return None
    return hashlib.sha256(target.read_bytes()).hexdigest()[:8]

changed = 0
for name in PAGES:
    page = ROOT / name
    text = page.read_text(encoding="utf-8")

    def repl(m):
        global changed
        h = digest(m.group(2), page)
        if h is None:
            return m.group(0)
        changed += 1
        return f"{m.group(1)}{m.group(2)}?v={h}{m.group(3)}"

    new = ASSET_RE.sub(repl, text)
    if new != text:
        page.write_text(new, encoding="utf-8")
    print(f"  {name}: {len(ASSET_RE.findall(text))} asset refs stamped")

print(f"stamped {changed} references")
