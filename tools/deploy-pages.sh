#!/usr/bin/env bash
# Deploy the public site to Cloudflare Pages (project `popkern`).
#
#   tools/deploy-pages.sh            # the pre-launch parking page + install.sh/.ps1
#
# Until launch the site is web/parking/ (a "coming soon" page, noindex)
# with web/lander/install.sh and install.ps1 at /install.sh and
# /install.ps1, served at
# https://popkern-c2q.pages.dev (Pages' own subdomain: popkern.pages.dev
# belongs to another account; popkern.dev is not attached yet). Release
# archives are not here: install.sh fetches them from
# https://dl.popkern.dev (the popkern-releases R2 bucket,
# tools/dist-linux.sh, tools/dist-windows.sh and tools/dist-macos.sh
# --upload). The full lander + live deploy is PK-832.
#
# Needs wrangler, logged in to the account that owns the project.
set -euo pipefail
cd "$(dirname "$0")/.."

PROJECT=popkern
STAGE=$(mktemp -d)
trap 'rm -rf "$STAGE"' EXIT

cp -R web/parking/. "$STAGE/"
cp web/lander/install.sh "$STAGE/install.sh"
cp web/lander/install.ps1 "$STAGE/install.ps1"
sh -n "$STAGE/install.sh"
# Parse-check install.ps1 where a PowerShell is at hand (it must also
# stay ASCII: irm decodes a charset-less response as Latin-1).
if LC_ALL=C grep -q "[^[:print:][:space:]]" "$STAGE/install.ps1"; then
  echo "install.ps1 has non-ASCII bytes" >&2; exit 1
fi
for ps in pwsh powershell; do
  if command -v "$ps" >/dev/null; then
    "$ps" -NoProfile -NonInteractive -Command '$e = $null; [void][System.Management.Automation.Language.Parser]::ParseFile($args[0], [ref]$null, [ref]$e); if ($e) { $e; exit 1 }' \
      "$(cygpath -w "$STAGE/install.ps1" 2>/dev/null || echo "$STAGE/install.ps1")"
    break
  fi
done

# Run from the stage so a stray wrangler config in the checkout never
# applies. wrangler 4.15x delegates `pages` commands to Workers static
# assets by default; --force keeps this a Pages project (*.pages.dev,
# _headers, Pages custom domains).
cd "$STAGE"
if ! wrangler pages project list 2>/dev/null | grep -qw "$PROJECT"; then
  wrangler pages project create "$PROJECT" --production-branch main --force
fi
wrangler pages deploy . --project-name "$PROJECT" --branch main --commit-dirty=true --force
cd - >/dev/null

URL=https://popkern-c2q.pages.dev
served() { # <file>: the deployed copy matches web/lander/<file>
  curl -fsS -o "$STAGE/served" "$URL/$1" && cmp -s "$STAGE/served" "web/lander/$1"
}
for attempt in 1 2 3 4 5 6; do
  if served install.sh && served install.ps1; then
    echo "  ok $URL/install.sh and install.ps1 match web/lander/"
    curl -fsSI "$URL/install.ps1" | grep -qi '^content-type: text/plain' &&
      echo "  ok $URL/install.ps1 is served as text (irm returns it as a string)"
    curl -fsS "$URL/" | grep -q "Coming soon" && echo "  ok $URL/ is the parking page"
    exit 0
  fi
  sleep 5
done
echo "FAIL $URL/install.sh or install.ps1 does not match web/lander/" >&2
exit 1
