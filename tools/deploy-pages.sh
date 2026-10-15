#!/usr/bin/env bash
# Deploy the public site to Cloudflare Pages (project `popkern`).
#
#   tools/deploy-pages.sh            # the pre-launch parking page + install.sh
#
# Until launch the site is web/parking/ (a "coming soon" page, noindex)
# with web/lander/install.sh at /install.sh, served at
# https://popkern-c2q.pages.dev (Pages' own subdomain: popkern.pages.dev
# belongs to another account; popkern.dev is not attached yet). Release
# archives are not here: install.sh fetches them from
# https://dl.popkern.dev (the popkern-releases R2 bucket,
# tools/dist-linux.sh --upload). The full lander + live deploy is PK-832.
#
# Needs wrangler, logged in to the account that owns the project.
set -euo pipefail
cd "$(dirname "$0")/.."

PROJECT=popkern
STAGE=$(mktemp -d)
trap 'rm -rf "$STAGE"' EXIT

cp -R web/parking/. "$STAGE/"
cp web/lander/install.sh "$STAGE/install.sh"
sh -n "$STAGE/install.sh"

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
for attempt in 1 2 3 4 5 6; do
  if curl -fsS -o "$STAGE/served.sh" "$URL/install.sh" && cmp -s "$STAGE/served.sh" web/lander/install.sh; then
    echo "  ok $URL/install.sh matches web/lander/install.sh"
    curl -fsS "$URL/" | grep -q "Coming soon" && echo "  ok $URL/ is the parking page"
    exit 0
  fi
  sleep 5
done
echo "FAIL $URL/install.sh does not match web/lander/install.sh" >&2
exit 1
