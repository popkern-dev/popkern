#!/usr/bin/env bash
# Regenerate the KP20 pyramid-apex showcase end to end:
#   1. run the example (STEP + OBJ exports),
#   2. mesh-finish the whole-object OBJ (weld + close sliver seams),
#   3. render the standard six views with f3d,
#   4. publish everything under a timestamped name set (defeats
#      caching between the file server and viewers) and prune older
#      stamped sets, keeping the most recent PREVIOUS set as backup.
set -euo pipefail
cd "$(dirname "$0")/.."

STAMP=$(date -u +%Y%m%d-%H%M%S)
OUT=examples/output
PREFIX="kp20_apex_${STAMP}"

cargo run --release --bin kp20_pyramid_apex_showcase

python3 tools/mesh-finish.py \
  "$OUT/kp20_pyramid_apex_whole.obj" "$OUT/${PREFIX}_whole.obj"

# Blend-detail companion (fillets + cap only), same finish.
python3 - "$OUT/kp20_pyramid_apex_whole.obj" "$OUT/.blend_detail_raw.obj" <<'EOF'
import sys
keep = ("fillet_", "paper_cap")
out, active = [], True
for line in open(sys.argv[1]):
    if line.startswith("g "):
        active = line[2:].strip().startswith(keep)
        out.append(line)
        continue
    if line.startswith("f ") and not active:
        continue
    out.append(line)
open(sys.argv[2], "w").writelines(out)
EOF
python3 tools/mesh-finish.py \
  "$OUT/.blend_detail_raw.obj" "$OUT/${PREFIX}_blend_detail.obj" || true
rm -f "$OUT/.blend_detail_raw.obj"

cp "$OUT/kp20_pyramid_apex_cap.step" "$OUT/${PREFIX}_cap.step"
cp "$OUT/kp20_pyramid_apex_cap.obj" "$OUT/${PREFIX}_cap.obj"

render() { # name, extra f3d args...
  local name=$1
  shift
  f3d --output "$OUT/${PREFIX}_${name}.png" --resolution 1024,768 \
    --up +Z -q -a -t "$@"
}
render whole --camera-direction=-1,-0.5,-0.35 "$OUT/${PREFIX}_whole.obj"
render whole_top --camera-direction=-0.05,-0.05,-1 "$OUT/${PREFIX}_whole.obj"
render bottom --camera-direction=0.3,0.3,0.85 "$OUT/${PREFIX}_whole.obj"
render closeup --camera-position=0.45,-0.45,1.65 \
  --camera-focal-point=0,0,1.05 "$OUT/${PREFIX}_whole.obj"
render blend_detail --camera-position=0.45,-0.45,1.65 \
  --camera-focal-point=0,0,1.05 "$OUT/${PREFIX}_blend_detail.obj"

# Prune: keep this set and one previous set.
stamps=$(ls "$OUT" | grep -oE '^kp20_apex_[0-9]{8}-[0-9]{6}' | sort -u | head -n -2 || true)
for old in $stamps; do
  rm -f "$OUT/${old}"_*
done

echo "published set: ${PREFIX}"
ls "$OUT" | grep "^${PREFIX}"
