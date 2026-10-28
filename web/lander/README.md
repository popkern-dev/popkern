# popkern.dev landing site

The marketing site for [popkern](https://github.com/popkern-dev/popkern), an
open boundary-representation CAD kernel written in Rust.

Plain static HTML, CSS, and one small JS file. **No build step, no dependencies,
no framework.** What is in this folder is exactly what gets served.

It lives in the popkern repository as `web/lander/` (since 2026-10-14; its
history before that is the `popkernlander` repo's, merged in). The site it
publishes is this folder at the root with popkern live (`web/live/`) under
`live/`, which is where every "Try it now" link points.

```
index.html                  the landing page
install/index.html          install instructions (install.sh / install.ps1)
install.sh                  Linux and macOS installer (curl ... | sh)
install.ps1                 Windows installer (irm ... | iex)
support/index.html          ways to support popkern
blog/index.html             post index
blog/announcing-popkern.html  launch post
blog/feed.xml               RSS (hand-maintained)
playground/index.html       redirect to live/ (the playground shipped as popkern live)
404.html
assets/css/site.css         the whole design system, one file
assets/js/site.js           copy buttons + footer year, nothing else
assets/img/                 gallery cutouts, OG card, favicon
_headers                    Cloudflare Pages cache/security headers
llms.txt, *.md              NOT in this folder: the deploy writes them (see below)
.nojekyll                   stops GitHub Pages running Jekyll over the tree
tools/                      logo, social-card and asset-stamp generators (not deployed)
```

## Before announcing

1. **Social links are set** (2026-10-24): GitHub, Discord, Reddit, Bluesky
   and X, in every page's nav (GitHub, Discord) and footer (all five), the
   homepage's Community column and "Come and break it" card, and the help
   lists on the support page and launch post. Any link still spelled
   `#TODO-…` renders with a small `?` marker (the `a[href^="#TODO"]` rule in
   `assets/css/site.css`).
2. **Set the domain.** `popkern.dev` is assumed in the canonical tag, the
   `og:url`/`og:image` tags, `sitemap.xml`, `robots.txt` and `blog/feed.xml`.
   Change it everywhere if the real domain differs:
   `grep -rl popkern.dev . --exclude-dir=.git`
3. **Check the launch post.** `blog/announcing-popkern.html` is a full draft,
   not a stub. Read it before it goes out. The date (28 Oct 2026) appears in
   the post, the blog index, and the feed.
4. **Fill in the sponsor links** on `support/index.html`
   (`#TODO-sponsor-monthly`, `#TODO-sponsor-once`).
5. **Ship `install.sh` and `install.ps1`.** The homepage, install page,
   support page and launch post link to them. Both are in this folder and
   `tools/deploy-pages.sh` already serves them (as `text/plain`, which `irm`
   needs to hand `iex` a string). They install the archives
   `tools/dist-linux.sh`, `tools/dist-windows.sh` and `tools/dist-macos.sh` publish to
   `https://dl.popkern.dev`. (`live/` ships with the site: see Deploying.)

## Deploying

`tools/deploy-site-routes.sh` (repo root) builds popkern live, copies this
folder to the site root and live to `live/`, and verifies the served copy. It
publishes to `~/code/.routes/static/popkern/`, which the routes server serves
at **`/routes/popkern/`** (locally on port 7681, and on the tunnel hostname).
The README, `tools/` and editor swap files stay behind.

The deploy also writes the site's agent-readable copies into the staged
tree with `tools/make-markdown.py`: a Markdown version of each page beside
it (`index.md`, `blog/announcing-popkern.md`, ...) and `/llms.txt`, an
index of them. Each page names its copy with
`<link rel="alternate" type="text/markdown">` and opens its body with a
hidden (`sr-only`, `aria-hidden`) pointer to `llms.txt`. A new page goes
in the script's `PAGES` list and gets the same two tags.

```sh
tools/deploy-site-routes.sh --no-build   # lander edits only: reuse the built live bundle
```

The public site, **popkern.dev**, is Cloudflare Pages (project `popkern`).
`tools/deploy-pages.sh` runs the routes deploy above first, then publishes
that same tree: the bakes go gzipped to the R2 bucket `popkern-site` (Pages
caps a file at 25 MiB; the arm's bake is 105 MB) and a Pages Function serves
them on the site's own origin. The script's header has the details.

```sh
tools/deploy-pages.sh --no-build         # lander edits only: reuse the built live bundle
```

> `404.html` is the one file using root-absolute paths (`/assets/…`), because
> a 404 can be served from any depth. Every other page uses relative paths and
> works at any prefix. On a project page (`<user>.github.io/<repo>/`), either
> use a custom domain or change those three paths in `404.html`.

## Local preview

Any static server works for the lander alone:

```sh
cd web/lander && python3 -m http.server 8000
```

The `live/` links 404 there; the deployed site has them.

## Where the images come from

The hero and gallery images are popkern live's own renders of its gallery
examples, so each one is a model anyone can open on the live page.
`tools/render-gallery.mjs` opens every example in `tools/gallery-jobs.json`
pre-built in headless Chrome, sets its view (Shade and Color, as on the page)
and camera, hides every card so only the transparent canvas is left, and
screenshots it at 2x. `tools/cutout.py` crops each to its subject and
downscales it (1400 px wide for the hero, 900 for the gallery) into
`assets/img/`. Then run `tools/stamp-assets.py` for fresh `?v=` hashes.

```sh
tools/build-live.sh --web && tools/bake-live.sh          # repo root
node web/lander/tools/render-gallery.mjs /tmp/renders "$POPKERN_CHROME"
python3 web/lander/tools/cutout.py /tmp/renders/arm.png web/lander/assets/img/six-axis-arm.png 1400
```

The colors mean what they mean on the page: a part's own finish (CAD), one
color per part (Parts), or zebra reflection lines across surface seams.

`assets/img/og.png` (1200x630 social card) is generated, not hand-drawn, by
`tools/make-og.py`: the six-axis arm cutout over a teal glow, the brand mark with its
amber arc, and the homepage tagline. It needs Pillow plus the Inter and
JetBrains Mono TTFs (from their GitHub releases):

```sh
uv run --with pillow python tools/make-og.py . <fonts-dir> assets/img/og.png
```

`<fonts-dir>` holds `inter/extras/ttf/Inter-{Bold,Regular}.ttf` and
`jb/fonts/ttf/JetBrainsMono-Regular.ttf`. Re-run it whenever the tagline or
palette changes.

## Cache busting

The sandbox is exposed through a Cloudflare Tunnel, and Cloudflare caches by
file extension at the edge. It will not re-read a `.css` or `.png` just because
the origin later starts sending `Cache-Control: no-cache`, so a stale
stylesheet can outlive several rounds of edits and make the live page look
unchanged.

With no build step there are no fingerprinted filenames, so after touching
anything in `assets/` run:

```sh
python3 tools/stamp-assets.py
```

It rewrites every local asset reference to `...?v=<sha256 prefix>`. A changed
file gets a new URL, which no cache anywhere has seen. The script is idempotent
and safe to re-run. If the HTML itself looks stale, load the page once with a
throwaway query (`/routes/popkern/?v=2`).

## Editing notes

- The nav and footer are duplicated in each page on purpose. With
  no build step and four pages, a copy is cheaper than a templating dependency.
  If a fifth page appears, reconsider.
- `assets/css/site.css` is token-driven. The palette lives in `:root`, with a
  full light-mode override under `prefers-color-scheme: light`: teal
  `--brand`/`--brand-ink` is the main colour, `--blue` the far end of the
  headline gradient, and amber `--amber` the logo arc and primary buttons.
  The amber logo arc is also hardcoded in three non-CSS places: the inline
  `brand-mark` SVG in each page's nav and footer, `assets/img/favicon.svg`, and
  the generated OG card.
- Fonts are Inter and JetBrains Mono from Google Fonts, with system fallbacks
  that hold the layout if the fetch fails.
