#!/bin/bash
# Couto silo with and without the compaction cap (WP K.5/K.7, decision 0124):
# the run_final_4mm.sh setup (gravity-on start, outlet opens at 0.5 s, run to
# 1.5 s, T0 rules + fullBlend + nonlocal), checkpointed and mirrored.
#
#     tests/hopper/run_cap_study.sh mc|dp cap|nocap coarse|fine
#
# coarse: 12 mm cylinder rows, 12 mm hopper rows (decision 0115 mesh);
# fine: the 4 mm mesh of decision 0117. Cap values (decision 0124):
# CEL_CAP_LAMBDA (default 0.023, wheat, Zeng and Wang 2019), CEL_CAP_P0 (100).
# Work: $CEL_WORK (default /tmp/cel_cap); report $CEL_WORK/<tag>.out.
# Mirror: $CEL_MIRROR (default /mnt/project-files/hopper/cap/ckpt). Run again
# after a lost container: it resumes from the latest checkpoint.
source /usr/lib/openfoam/openfoam2512/etc/bashrc
cd "$(dirname "$0")/../.." || exit 1
model=$1 cap=$2 grid=$3
tag=${model}_${cap}_${grid}
work=${CEL_WORK:-/tmp/cel_cap}/$tag
mirror=${CEL_MIRROR:-/mnt/project-files/hopper/cap/ckpt}/$tag
name=couto_${model}_aligned_buffered
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

if [ "$grid" = fine ]; then
    export CEL_H_TOP=0.012 CEL_NCOLS=58 CEL_WALL_RATIO=0.25 \
           CEL_H_HOP=0.004 CEL_OUTLET_WALL_RATIO=1.0
else
    export CEL_H_TOP=0.012
fi
export CEL_NPROCS=${CEL_NPROCS:-1} CEL_DP_FIT=plane_strain
export CEL_CHECKPOINT=${CEL_CHECKPOINT:-0.05} CEL_RESUME=1
rules='mu2 0.8412;\nI0 0.279;\nd 0.004;\nrhoS 1380;\nshearCap yes;\nspreadLayers yes;\nflowImpedance yes;\nnormalCap yes;\nfullBlend yes;\nnonlocalA 0.48;\n'
if [ "$cap" = cap ]; then
    rules="${rules}capLambda ${CEL_CAP_LAMBDA:-0.023};\ncapP0 ${CEL_CAP_P0:-100};\n"
fi
export CEL_HOPPER_RULES="$(printf "$rules")"
python3 tests/hopper/hopper_pressure_checks.py "$work" "$model" aligned_buffered \
    >> "$work/$tag.out" 2>&1
status=$?
kill $syncer 2>/dev/null
sync_case
echo "EXIT $status" >> "$work/$tag.out"
cp "$work/$tag.out" "$mirror/" 2>/dev/null
