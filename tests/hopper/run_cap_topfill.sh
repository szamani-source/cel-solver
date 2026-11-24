#!/bin/bash
# Couto silo filled from the top with the compaction cap (WP K.5, decision
# 0124): loose wheat (CEL_RHO, default 772 kg/m^3, the standard test weight
# in Bhadra et al. 2018) fed until the Couto mass (1503 kg) is in, cap
# lambda* 0.023 (Zeng and Wang 2019), restRelax, 2 s settle. Same as the
# decision 0115 run 2 otherwise, which is the no-cap reference.
#
#     tests/hopper/run_cap_topfill.sh          (about 8 h on 2 processes)
#
# Checkpoints every 0.05 s, mirrored to $CEL_MIRROR
# (/mnt/project-files/hopper/cap/topfill); run again to resume.
source /usr/lib/openfoam/openfoam2512/etc/bashrc
cd "$(dirname "$0")/../.." || exit 1
work=${CEL_WORK:-/tmp/cel_topfill}
mirror=${CEL_MIRROR:-/mnt/project-files/hopper/cap/topfill}
name=couto_mc_topfill
case=$work/$name
mkdir -p "$work" "$mirror"
m=$mirror/$name
[ -d "$m" ] || { [ -d "$m.old" ] && mv "$m.old" "$m"; }
if [ ! -f "$case/history.dat" ] && [ -f "$m/history.dat" ]; then
    echo "restoring $case from $mirror"
    cp -a "$m" "$work/"
fi
sync_case() {
    [ -d "$case" ] || return
    for pd in "$case" "$case"/processor*; do
        [ -d "$pd" ] || continue
        (cd "$pd" && ls -d [0-9]* 2>/dev/null | grep -vx 0 | sort -g | head -n -2 \
            | xargs -r rm -rf)
    done
    local dst=$mirror/$name
    rm -rf "$dst.new" "$dst.old"
    cp -a "$case" "$dst.new" || return
    [ -d "$dst" ] && mv "$dst" "$dst.old"
    mv "$dst.new" "$dst" && rm -rf "$dst.old"
}
( while true; do sleep 300; sync_case; done ) &
syncer=$!
export CEL_NPROCS=${CEL_NPROCS:-2} CEL_RHO=${CEL_RHO:-772} CEL_SETTLE=${CEL_SETTLE:-2.0}
export CEL_FILL_MASS=${CEL_FILL_MASS:-1503} CEL_HEAD=${CEL_HEAD:-0.60}
export CEL_CHECKPOINT=${CEL_CHECKPOINT:-0.05} CEL_RESUME=1
rules='mu2 0.8412;\nI0 0.279;\nd 0.004;\nrhoS 1380;\nshearCap yes;\nspreadLayers yes;\nflowImpedance yes;\nnormalCap yes;\nfullBlend yes;\n'
rules="${rules}restRelax { time 0.05; speed 0.01; }\ncapLambda ${CEL_CAP_LAMBDA:-0.023};\ncapP0 ${CEL_CAP_P0:-100};\n"
export CEL_HOPPER_RULES="$(printf "$rules")"
python3 tests/hopper/top_fill.py "$work" >> "$work/topfill.out" 2>&1
status=$?
kill $syncer 2>/dev/null
sync_case
echo "EXIT $status" >> "$work/topfill.out"
cp "$work/topfill.out" "$mirror/" 2>/dev/null
