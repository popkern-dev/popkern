# Release-bucket upload helpers shared by tools/dist-linux.sh,
# tools/dist-windows.sh and tools/dist-macos.sh. Source it; it runs nothing by itself.
#
# Each platform builds on its own machine and publishes into the same
# v<version>/ folder of the popkern-releases R2 bucket, served at
# https://dl.popkern.dev. Archives are written once (immutable).
# v<version>/SHA256SUMS is shared: every platform's publish merges its own
# lines into the copy already published, so the second platform to publish
# a release keeps the first one's checksums. That makes SHA256SUMS mutable
# for as long as a release is still gathering platforms, so it is cached
# for minutes, not a year.
#
# Needs wrangler, logged in to the account that owns the bucket (or
# CLOUDFLARE_API_TOKEN + CLOUDFLARE_ACCOUNT_ID for unattended runs). Run
# it from the checkout, not $HOME: wrangler keeps a per-project cache in
# ./.wrangler, and in $HOME that folder shadows its own login config.

DIST_BUCKET=${DIST_BUCKET:-popkern-releases}
DIST_BASE=${DIST_BASE:-https://dl.popkern.dev}
DIST_IMMUTABLE="public, max-age=31536000, immutable"
DIST_SHORT="public, max-age=300"

# dist_put <key> <file> <content-type> <cache-control>
dist_put() {
  wrangler r2 object put "$DIST_BUCKET/$1" --remote --file "$2" \
    --content-type "$3" --cache-control "$4" >/dev/null
  echo "uploaded $DIST_BASE/$1"
}

# dist_publish_release <out-dir> <version> <latest:0|1> <archive>...
# Uploads the archives, merges their checksums from <out-dir>/SHA256SUMS
# into the published v<version>/SHA256SUMS, and with latest=1 points
# `latest` at <version>. The merged file is left at <out-dir>/SHA256SUMS.published.
dist_publish_release() {
  local out=$1 version=$2 latest=$3; shift 3
  local file
  # wrangler reports a missing bucket as a missing key, which the merge
  # would read as a new release: confirm the bucket (and the login) first.
  if ! wrangler r2 bucket info "$DIST_BUCKET" >/dev/null 2>&1; then
    echo "cannot see bucket $DIST_BUCKET: is wrangler logged in to its account? (wrangler whoami)" >&2
    return 1
  fi
  for file in "$@"; do
    dist_put "v$version/$(basename "$file")" "$file" "$(dist_content_type "$file")" "$DIST_IMMUTABLE"
  done
  dist_merge_sums "$out" "$version"
  dist_put "v$version/SHA256SUMS" "$out/SHA256SUMS.published" "text/plain; charset=utf-8" "$DIST_SHORT"
  if [[ $latest == 1 ]]; then
    printf '%s\n' "$version" > "$out/latest"
    dist_put latest "$out/latest" "text/plain; charset=utf-8" "public, max-age=60"
  fi
}

dist_content_type() {
  case "$1" in
    *.tar.gz) echo application/gzip ;;
    *.zip) echo application/zip ;;
    *) echo application/octet-stream ;;
  esac
}

# dist_merge_sums <out-dir> <version>: the published SHA256SUMS, read
# from the bucket itself (no CDN copy in between), with every line naming
# one of this build's archives replaced by this build's, sorted by file
# name. A missing key means a new release; any other failure stops the
# publish rather than risk replacing another platform's lines with nothing.
dist_merge_sums() {
  local out=$1 version=$2
  if ! wrangler r2 object get "$DIST_BUCKET/v$version/SHA256SUMS" --remote --pipe \
      > "$out/SHA256SUMS.remote" 2> "$out/SHA256SUMS.err"; then
    if grep -q 'specified key does not exist' "$out/SHA256SUMS.err"; then
      : > "$out/SHA256SUMS.remote"
    else
      echo "cannot read $DIST_BUCKET/v$version/SHA256SUMS:" >&2
      cat "$out/SHA256SUMS.err" >&2
      return 1
    fi
  fi
  rm -f "$out/SHA256SUMS.err"
  awk 'NR == FNR { mine[$2] = 1; print; next } !($2 in mine) && NF == 2' \
    "$out/SHA256SUMS" "$out/SHA256SUMS.remote" | sort -k2 > "$out/SHA256SUMS.published"
  rm -f "$out/SHA256SUMS.remote"
}
