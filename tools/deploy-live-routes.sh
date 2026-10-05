#!/usr/bin/env bash
# Deploy popkern live (the browser workbench) to the sandbox's /routes/ web server
# (PK-426, design §10 phase 6e).
#
#   tools/deploy-live-routes.sh              # build, bake, rsync, verify
#   tools/deploy-live-routes.sh --no-build   # redeploy the existing bundle
#   tools/deploy-live-routes.sh --no-verify  # skip the headless-Chrome pass
#
# The sandbox's Caddy (see ~/Caddyfile) serves ~/code/.routes/static/<name>/
# at https://<host>/routes/<name>/ — nginx fronts it on port 7681 and the
# Cloudflare tunnel points at that port, so the deployed page is reachable
# from any browser that can reach the tunnel hostname. The static tree
# persists across container rebuilds; no Caddy reload is needed for it.
#
# Only the page's runtime assets are deployed: the node test harness, the
# wasm-bindgen .d.ts files and npm packaging metadata stay behind. The page
# uses only relative URLs and meta-tag CSP, so it is subpath-agnostic.
#
# Daemon variant (phase 6b+): a live popkern-daemon session is wired with a
# reverse-proxy snippet instead of a static copy — drop
#   handle_path /routes/<name>/* {
#       reverse_proxy 127.0.0.1:<port>
#   }
#   redir /routes/<name> /routes/<name>/
# into ~/code/.routes/conf.d/<name>.caddy and run
# `caddy reload --config ~/Caddyfile`. (handle_path strips the prefix, so
# the daemon serves "/" and never needs to know where it is mounted; the
# redir catches the slashless spelling, which otherwise falls through to
# the static file server and 404s.)
set -euo pipefail
cd "$(dirname "$0")/.."

BUILD=1
VERIFY=1
for arg in "$@"; do
  case "$arg" in
    --no-build) BUILD=0 ;;
    --no-verify) VERIFY=0 ;;
    *) echo "unknown argument: $arg" >&2; exit 2 ;;
  esac
done

ROUTE_NAME="${POPKERN_ROUTE_NAME:-popkern-live}"
# Earlier name of the route (the page was "popkern playground" until
# 2026-10-03): it keeps answering, as a redirect to the current one.
LEGACY_ROUTE_NAMES=(popkern-playground)
ROUTES_DIR="${ROUTES_DIR:-$HOME/code/.routes}"
DEST="$ROUTES_DIR/static/$ROUTE_NAME"
LOCAL_BASE="${POPKERN_ROUTES_BASE:-http://127.0.0.1:7681}"
URL="$LOCAL_BASE/routes/$ROUTE_NAME"

if [[ ! -d "$ROUTES_DIR/static" ]]; then
  echo "no sandbox routes tree at $ROUTES_DIR/static — is this the ccvm sandbox?" >&2
  exit 1
fi

[[ $BUILD == 1 ]] && tools/build-live.sh --web
# Pre-built gallery examples (reused while the bundle and scripts are
# unchanged; a new bundle rebakes them all, about 90 s).
[[ $BUILD == 1 ]] && tools/bake-live.sh

if [[ ! -f web/live/pkg/popkern_live_bg.wasm ]]; then
  echo "no wasm bundle at web/live/pkg/ — run without --no-build" >&2
  exit 1
fi

# Versioned deploy (2026-09-24): every build lands whole in an immutable
# content-addressed directory v/<stamp>/, and the mutable surface is two
# root entry pages that carry a <base> into it. A cache layer (browser
# heuristic cache, tunnel edge) that holds one file but not its siblings
# can therefore never pair JS glue with a wasm module from a different
# build — that mix traps as `unreachable executed` at kernel init, which
# is what an in-place deploy risks on every ship. A stale root page is
# still coherent: it points at the previous version directory, which is
# kept (last 3 builds survive pruning).
#
# The stamp hashes every staged file (path and content), not only the
# wasm: a CSS/JS/HTML-only change must land in a new directory too, or
# the "immutable" tree would change under an unchanged name.

# tar pipe (the sandbox has no rsync) into a fresh tree, then swap it
# in whole — a half-copied deploy is never what Caddy serves.
mkdir -p "$DEST/v"
STAGE="$DEST/v/.stage-$$"
rm -rf "$STAGE"
mkdir -p "$STAGE"
tar -C web/live \
  --exclude './test' \
  --exclude './core.*' \
  --exclude './package.json' \
  --exclude './pkg/*.d.ts' \
  --exclude './pkg/package.json' \
  -cf - . | tar -C "$STAGE" -xf -
STAMP=$( (cd "$STAGE" && find . -type f -print0 | LC_ALL=C sort -z | xargs -0 sha256sum) \
          | sha256sum | cut -c1-12)
VDIR="$DEST/v/$STAMP"
rm -rf "$VDIR"
mv "$STAGE" "$VDIR"
# tar restores the source tree's mtimes, so the new version directory
# would look as old as web/live/ and the newest-3 pruning below
# could delete the very version the root pages are about to point at
# (it did, 2026-09-26: the root page served a 404 wasm). Stamp it now.
touch "$VDIR"

# Root entry pages: the version's own pages with a <base> injected right
# after <head>, so every relative URL — module imports, the worker, the
# examples manifest, the wasm fetch — resolves inside the immutable tree.
# Requires base-uri 'self' in the pages' CSP (index.html documents why).
for page in index.html guide.html; do
  awk -v tag="<base href=\"./v/$STAMP/\">" '{print} $0 == "<head>" {print tag}' \
    "$VDIR/$page" > "$DEST/.$page.tmp"
  mv "$DEST/.$page.tmp" "$DEST/$page"
done

# Clear root-level files from pre-versioned deploys (everything now lives
# under v/), then keep only the last 3 version directories.
find "$DEST" -mindepth 1 -maxdepth 1 ! -name v ! -name index.html ! -name guide.html \
  -exec rm -rf {} +
ls -dt "$DEST"/v/*/ 2>/dev/null | tail -n +4 | xargs -r rm -rf

# Legacy routes: entry pages that redirect here (meta refresh: the CSP
# allows no script), and no app of their own, so an old bookmark lands
# on the current build instead of a stale one.
for legacy in "${LEGACY_ROUTE_NAMES[@]}"; do
  [[ "$legacy" == "$ROUTE_NAME" ]] && continue
  LEGACY_DEST="$ROUTES_DIR/static/$legacy"
  mkdir -p "$LEGACY_DEST"
  find "$LEGACY_DEST" -mindepth 1 -maxdepth 1 -exec rm -rf {} +
  for page in index.html guide.html; do
    target="../$ROUTE_NAME/"
    [[ $page == guide.html ]] && target="../$ROUTE_NAME/guide.html"
    cat > "$LEGACY_DEST/$page" <<HTML
<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta http-equiv="Content-Security-Policy" content="default-src 'none'; base-uri 'none'; form-action 'none'">
<meta http-equiv="refresh" content="0; url=$target">
<title>popkern live</title>
</head>
<body><p>popkern playground is now <a href="$target">popkern live</a>.</p></body>
</html>
HTML
  done
done

echo "deployed $STAMP to $VDIR"
echo "served at $URL/  (and at /routes/$ROUTE_NAME/ on the tunnel hostname)"

# ---- Served-URL checks: reachable, typed, and self-contained ----
check() { # path expected-type
  local code type
  read -r code type < <(curl -s -o /dev/null -w '%{http_code} %{content_type}\n' "$URL/$1")
  if [[ "$code" != 200 || "$type" != "$2"* ]]; then
    echo "FAIL $URL/$1 -> $code $type (wanted 200 $2)" >&2
    return 1
  fi
  echo "  ok $1 ($type)"
}
check "" "text/html"
check "guide.html" "text/html"
check "v/$STAMP/pkg/popkern_live_bg.wasm" "application/wasm"
check "v/$STAMP/pkg/popkern_live.js" "text/javascript"
check "v/$STAMP/examples.json" "application/json"
check "v/$STAMP/examples/boss_plate.rhai" ""
check "v/$STAMP/examples/boss_plate_in_steps/migrations/001_plate.rhai" ""
check "v/$STAMP/guide.html" "text/html"

# The root entry pages must reference the version that was just deployed.
for page in "" guide.html; do
  if ! curl -s "$URL/$page" | grep -q "<base href=\"./v/$STAMP/\">"; then
    echo "FAIL $URL/$page does not point at v/$STAMP" >&2
    exit 1
  fi
done
echo "  ok root pages point at v/$STAMP"

# The root entry pages are the only mutable files, so the server must
# tell browsers to revalidate them. Without Cache-Control: no-cache,
# heuristic caching (10% of Last-Modified age) showed a days-stale
# example manifest after a redeploy (2026-09-06). The header lives in
# ~/Caddyfile's /routes/ handler; fail loudly if this sandbox's server
# lost it. (Everything under v/ is immutable, so caching it is safe.)
if ! curl -sI "$URL/" | grep -qi '^cache-control: *no-cache'; then
  echo "FAIL $URL/ served without Cache-Control: no-cache —" >&2
  echo "     browsers will cache a stale entry page; add the header to the" >&2
  echo "     static /routes/ handler in ~/Caddyfile and caddy reload" >&2
  exit 1
fi
echo "  ok root page (Cache-Control: no-cache)"

if [[ $VERIFY == 1 ]]; then
  # The DevTools harness against the deployed copy: real kernel run,
  # face pick, parameter edit, cancel path, request accounting.
  OUT_DIR=docs/reports/live
  mkdir -p "$OUT_DIR"
  POPKERN_PAGE_ORIGIN="$URL" node web/live/test/browser.mjs \
    "$OUT_DIR/routes-deploy.png" > "$OUT_DIR/routes-deploy.json"
  node -e '
    const out = require("./docs/reports/live/routes-deploy.json");
    const fail = (m) => { console.error("FAIL: " + m); process.exit(1); };
    if (out.error) fail(out.error);
    if (!/^served/.test(out.status_after_first_run ?? "")) fail("first run: " + out.status_after_first_run);
    if (out.bakes_present && !/^served \(pre-built\)/.test(out.status_after_first_run)) fail("the default script built live on load: its bake is missing");
    if (out.live_run && out.live_run.outcome !== "served") fail("live run of the default script: " + JSON.stringify(out.live_run));
    if (out.foreign_requests.length) fail("foreign requests: " + out.foreign_requests.join(", "));
    if (out.exceptions.length) fail("page exceptions: " + JSON.stringify(out.exceptions[0]));
    if (out.server_requests.length) fail("harness server was hit; the deploy was not what ran");
    console.log("  ok harness: " + out.status_after_first_run
      + " | kernel ready " + out.kernel_ready_ms + " ms"
      + " | first render " + out.first_render_ms + " ms"
      + " | " + out.bodies_on_screen + " body on screen");
  '
  echo "  report: $OUT_DIR/routes-deploy.json + routes-deploy.png"
fi
