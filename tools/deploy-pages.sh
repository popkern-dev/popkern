#!/usr/bin/env bash
# Deploy popkern.dev: the lander with popkern live under live/, on
# Cloudflare Pages (project `popkern`, custom domains popkern.dev and
# www.popkern.dev; PK-832).
#
#   tools/deploy-pages.sh              # build, bake, stage, publish, verify
#   tools/deploy-pages.sh --no-build   # publish the existing bundle and bakes
#   tools/deploy-pages.sh --no-verify  # skip the headless-Chrome passes
#
# What goes public is what tools/deploy-site-routes.sh (run first, with
# the same flags) has just staged and checked on the routes server: the
# lander at the root and live/ holding the current version directory
# v/<stamp>/ (older versions stay on the routes server only; Pages
# revalidates every page, so no browser holds a stale entry page).
#
# Bakes are not Pages assets. Pages caps a file at 25 MiB (the arm's bake
# is 105 MB), and Cloudflare does not compress application/octet-stream,
# so every bakes/*.bin goes gzipped (about 3:1) to the R2 bucket
# `popkern-site` under its own path, and a Pages Function on
# live/v/<stamp>/bakes/* streams it back with Content-Encoding: gzip on
# the page's own origin: the page's CSP is connect-src 'self', and
# fetch() decodes it before bakes.js sees a byte. Anything the bucket
# does not hold (bakes/index.json) falls through to the static assets.
#
# Release archives are not here: install.sh fetches them from
# https://dl.popkern.dev (the popkern-releases R2 bucket,
# tools/dist-linux.sh, tools/dist-windows.sh and tools/dist-macos.sh
# --upload).
#
# Needs wrangler, logged in to the account that owns the project and
# both buckets, and POPKERN_CHROME (a headless Chrome) unless
# --no-verify. The custom domains were attached once through the API
# (POST /accounts/<id>/pages/projects/popkern/domains); Cloudflare made
# their DNS records, the zone being in the same account.
set -euo pipefail
cd "$(dirname "$0")/.."

PROJECT=popkern
BUCKET=popkern-site
SITE_URL="${POPKERN_SITE_URL:-https://popkern.dev}"
VERIFY=1
for arg in "$@"; do
  case "$arg" in
    --no-build) ;;
    --no-verify) VERIFY=0 ;;
    *) echo "unknown argument: $arg" >&2; exit 2 ;;
  esac
done

tools/deploy-site-routes.sh "$@"

ROUTES_DIR="${ROUTES_DIR:-$HOME/code/.routes}"
SITE="$ROUTES_DIR/static/${POPKERN_ROUTE_NAME:-popkern}"
STAMP=$(sed -nE 's|^<base href="\./v/([0-9a-f]+)/">$|\1|p' "$SITE/live/index.html")
if [[ -z "$STAMP" || ! -d "$SITE/live/v/$STAMP" ]]; then
  echo "no current live version under $SITE/live" >&2
  exit 1
fi

# install.ps1 must stay ASCII (irm decodes a charset-less response as
# Latin-1) and parse; install.sh must parse.
sh -n web/lander/install.sh
if LC_ALL=C grep -q "[^[:print:][:space:]]" web/lander/install.ps1; then
  echo "install.ps1 has non-ASCII bytes" >&2; exit 1
fi
for ps in pwsh powershell; do
  if command -v "$ps" >/dev/null; then
    "$ps" -NoProfile -NonInteractive -Command '$e = $null; [void][System.Management.Automation.Language.Parser]::ParseFile($args[0], [ref]$null, [ref]$e); if ($e) { $e; exit 1 }' \
      "$(cygpath -w "$PWD/web/lander/install.ps1" 2>/dev/null || echo "$PWD/web/lander/install.ps1")"
    break
  fi
done

# ---- Stage: the routes copy, current version only, bakes set aside ----
STAGE=$(mktemp -d)
trap 'rm -rf "$STAGE"' EXIT
PUBLIC="$STAGE/public"
mkdir -p "$PUBLIC/live/v" "$STAGE/bakes"
tar -C "$SITE" --exclude './live/v' -cf - . | tar -C "$PUBLIC" -xf -
cp -R "$SITE/live/v/$STAMP" "$PUBLIC/live/v/$STAMP"
BAKES="$PUBLIC/live/v/$STAMP/bakes"
shopt -s nullglob
for bake in "$BAKES"/*.bin; do
  gzip -9 -n -c "$bake" > "$STAGE/bakes/$(basename "$bake").gz"
  sha256sum "$bake" | sed "s|  .*/|  |" >> "$STAGE/bakes/SHA256SUMS"
  rm "$bake"
done
shopt -u nullglob
too_big=$(find "$PUBLIC" -type f -size +25M)
if [[ -n "$too_big" ]]; then
  echo "files over Pages' 25 MiB limit:" >&2; echo "$too_big" >&2; exit 1
fi

mkdir -p "$STAGE/functions/live/v/[stamp]/bakes"
cat > "$STAGE/functions/live/v/[stamp]/bakes/[name].js" <<'JS'
// A bake from the popkern-site bucket, stored gzipped (tools/deploy-pages.sh).
export async function onRequest({ request, env }) {
  if (request.method !== "GET" && request.method !== "HEAD") {
    return new Response(null, { status: 405, headers: { allow: "GET, HEAD" } });
  }
  const key = new URL(request.url).pathname.replace(/^\/+/, "");
  const object = request.method === "HEAD" ? await env.BAKES.head(key) : await env.BAKES.get(key);
  if (!object) return env.ASSETS.fetch(request);
  const headers = new Headers();
  object.writeHttpMetadata(headers);
  headers.set("etag", object.httpEtag);
  headers.set("x-content-type-options", "nosniff");
  const body = request.method === "HEAD" ? null : object.body;
  return new Response(body, { headers, encodeBody: "manual" });
}
JS
cat > "$PUBLIC/_routes.json" <<JSON
{ "version": 1, "include": ["/live/v/$STAMP/bakes/*"], "exclude": ["/live/v/$STAMP/bakes/index.json"] }
JSON
cat > "$STAGE/wrangler.toml" <<TOML
name = "$PROJECT"
pages_build_output_dir = "public"
compatibility_date = "2026-10-01"

[[r2_buckets]]
binding = "BAKES"
bucket_name = "$BUCKET"
TOML

# ---- Upload the bakes (immutable: their path carries the version) ----
cd "$STAGE"
if ! wrangler r2 bucket info "$BUCKET" >/dev/null 2>&1; then
  wrangler r2 bucket create "$BUCKET"
fi
for gz in bakes/*.gz; do
  name=$(basename "$gz" .gz)
  wrangler r2 object put "$BUCKET/live/v/$STAMP/bakes/$name" --remote --file "$gz" \
    --content-type application/octet-stream --content-encoding gzip \
    --cache-control "public, max-age=31536000, immutable" >/dev/null
  echo "  ok r2 live/v/$STAMP/bakes/$name ($(stat -c%s "$gz") bytes gzipped)"
done

# ---- Publish ----
# wrangler 4.15x delegates `pages` commands to Workers static assets by
# default; --force keeps this a Pages project (*.pages.dev, _headers,
# Pages custom domains). Run from the stage: its wrangler.toml is the
# project's configuration (bindings included).
if ! wrangler pages project list 2>/dev/null | grep -qw "$PROJECT"; then
  wrangler pages project create "$PROJECT" --production-branch main --force
fi
wrangler pages deploy --project-name "$PROJECT" --branch main \
  --commit-hash "$(git -C "$OLDPWD" rev-parse HEAD)" --commit-dirty=true --force
cd "$OLDPWD"

# ---- Served-URL checks ----
fetch_ok() { # path expected-type
  local code type
  read -r code type < <(curl -sL -o /dev/null -w '%{http_code} %{content_type}\n' "$SITE_URL/$1")
  [[ "$code" == 200 && "$type" == "$2"* ]]
}
check() { # path expected-type
  for attempt in 1 2 3 4 5 6; do
    if fetch_ok "$1" "$2"; then echo "  ok $SITE_URL/$1"; return 0; fi
    sleep 5
  done
  echo "FAIL $SITE_URL/$1 (wanted 200 $2)" >&2
  return 1
}
check "" "text/html"
check "install/" "text/html"
check "blog/announcing-popkern.html" "text/html"
check "live/" "text/html"
check "live/v/$STAMP/pkg/popkern_live_bg.wasm" "application/wasm"
check "live/v/$STAMP/bakes/index.json" "application/json"
check "install.sh" "text/plain"
check "install.ps1" "text/plain"
check "llms.txt" "text/plain"
check "index.md" "text/markdown"
check "scripting-guide.md" "text/markdown"
check "blog/announcing-popkern.md" "text/markdown"
for file in install.sh install.ps1; do
  curl -fsS "$SITE_URL/$file" | cmp -s - "web/lander/$file" ||
    { echo "FAIL $SITE_URL/$file does not match web/lander/$file" >&2; exit 1; }
done
echo "  ok install.sh and install.ps1 match web/lander/"
entry=$(curl -fsS "$SITE_URL/live/")
grep -q "<base href=\"./v/$STAMP/\">" <<<"$entry" ||
  { echo "FAIL $SITE_URL/live/ does not point at v/$STAMP" >&2; exit 1; }
echo "  ok live/ points at v/$STAMP"
# Every bake, decoded as a browser would, byte for byte.
while read -r sum name; do
  got=$(curl -fsS --compressed "$SITE_URL/live/v/$STAMP/bakes/$name" | sha256sum | cut -c1-64)
  [[ "$got" == "$sum" ]] || { echo "FAIL bake $name served wrong bytes" >&2; exit 1; }
done < "$STAGE/bakes/SHA256SUMS"
echo "  ok $(wc -l < "$STAGE/bakes/SHA256SUMS") bakes served from R2, byte-identical"

if [[ $VERIFY == 1 ]]; then
  OUT_DIR=docs/reports/live
  mkdir -p "$OUT_DIR"
  POPKERN_PAGE_ORIGIN="$SITE_URL/live" node web/live/test/browser.mjs \
    "$OUT_DIR/pages-deploy.png" > "$OUT_DIR/pages-deploy.json"
  node -e '
    const out = require("./docs/reports/live/pages-deploy.json");
    const fail = (m) => { console.error("FAIL: " + m); process.exit(1); };
    if (out.error) fail(out.error);
    if (!/^served \(pre-built\)/.test(out.status_after_first_run ?? "")) fail("first run: " + out.status_after_first_run);
    if (out.foreign_requests.length) fail("foreign requests: " + out.foreign_requests.join(", "));
    if (out.exceptions.length) fail("page exceptions: " + JSON.stringify(out.exceptions[0]));
    console.log("  ok harness on " + process.argv[1] + ": " + out.status_after_first_run
      + " | first render " + Math.round(out.first_render_ms) + " ms");
  ' "$SITE_URL"
fi
echo "published v/$STAMP to $SITE_URL"
