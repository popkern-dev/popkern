# popkern.dev landing site

The marketing site for [popkern](https://github.com/graemedouglas/popkern), an
open boundary-representation CAD kernel written in Rust.

Plain static HTML, CSS, and one small JS file. **No build step, no dependencies,
no framework.** What is in this repo is exactly what gets served, which is what
makes it deployable to GitHub Pages or Cloudflare Pages without a toolchain.

```
index.html                  the landing page
install/index.html          install instructions (install.sh / install.ps1)
support/index.html          ways to support popkern
blog/index.html             post index
blog/announcing-popkern.html  launch post
blog/feed.xml               RSS (hand-maintained)
playground/index.html       "coming soon" page for the wasm playground
404.html
assets/css/site.css         the whole design system, one file
assets/js/site.js           copy buttons + footer year, nothing else
assets/img/                 gallery cutouts, OG card, favicon
_headers                    Cloudflare Pages cache/security headers
.github/workflows/pages.yml GitHub Pages deploy
.nojekyll                   stops Pages running Jekyll over the tree
```

## Before announcing

1. **Fill in the social links.** Search the tree for `#TODO-`:
   `#TODO-x`, `#TODO-bluesky`, `#TODO-discord`. Unset links render with a small
   `?` marker so they are obvious in the page; that marker comes from the
   `a[href^="#TODO"]` rule in `assets/css/site.css` and disappears by itself
   once the hrefs are real.
2. **Set the domain.** `popkern.dev` is assumed in the canonical tag, the
   `og:url`/`og:image` tags, `sitemap.xml`, `robots.txt` and `blog/feed.xml`.
   Change it everywhere if the real domain differs:
   `grep -rl popkern.dev . --exclude-dir=.git`
3. **Check the launch post.** `blog/announcing-popkern.html` is a full draft,
   not a stub. Read it before it goes out. The date (30 Oct 2026) appears in
   the post, the blog index, and the feed.
4. **Fill in the sponsor links** on `support/index.html`
   (`#TODO-sponsor-monthly`, `#TODO-sponsor-once`).
5. **Ship `live/`, `install.sh` and `install.ps1`.** The homepage, install
   page, support page and launch post all link to them.

## Deploying

**Cloudflare Pages**: connect the repo; framework preset *None*; build command
empty; build output directory `/`. `_headers` is picked up automatically.

**GitHub Pages**: Settings, then Pages, then Source: **GitHub Actions**. The included
workflow uploads the repo root on every push to `main`.

> If you deploy to a *project* page (`<user>.github.io/<repo>/`) rather than a
> custom domain, note that `404.html` is the one file using root-absolute paths
> (`/assets/…`), because a 404 can be served from any depth. Every other page
> uses relative paths and works at any prefix. On a project page, either use a
> custom domain or change those three paths in `404.html`.

## Local preview

Served in this sandbox at **`/routes/popkern/`** by
`~/code/.routes/conf.d/popkern.caddy`, which points straight at this directory
and sends `Cache-Control: no-cache` so a stale stylesheet never survives a
design pass. Edits are live on the next request; a change to the snippet itself
needs `caddy reload --config ~/Caddyfile`.

(A symlink at `~/code/.routes/static/popkern` also exists and would work on its
own via the generic static fallback, but that path caches.)

Anywhere else, any static server works:

```sh
python3 -m http.server 8000
```

## Where the images come from

The gallery images are the kernel's own debug renders from
`popkern/examples/output/`, run through a small script that keys out the flat
grey background, crops to the subject, and rescales. The originals are produced
by `popkern-viz`'s headless CPU renderer as part of the example runs, so the
colours mean something: they mark face provenance, not decoration.

Two of them (`torus-boolean`, `shelled-torus`) were lilac in the source render
and got re-hued to the site blue. Only saturated pixels move, so the grey
bodies and the thin construction lines are untouched.

To refresh or swap one, re-run the cutout step against a different render. The
transform is: sample the modal corner colour, alpha out anything within ~26 of
it with a 14-wide soft ramp, crop to the alpha bbox with 8px padding, cap width
at 1100px.

`assets/img/og.png` (1200x630 social card) is generated, not hand-drawn, by
`tools/make-og.py`: the torus cutout over a teal glow, the brand mark with its
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
