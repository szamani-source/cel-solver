#!/bin/bash
# Mesh refinement of the Couto silo with and without the nonlocal law
# (WP N.4, decision 0119): the decision 0117 set-up (gravity-on start,
# outlet opens at 0.5 s, run to 1.5 s, 4 MPI ranks) on three meshes of one
# blockMesh family, hopper cells 8, 5.7 and 4 mm (outlet about 6.2, 4.4 and
# 3.1 mm; 6,204, 11,250 and 20,398 cells).
#
#     tests/hopper/run_refine.sh mc|dp 8|5.7|4 nl|local
#
# About 1, 3 and 8 h on 4 cores. nl: the final validation rules (decision
# 0117, nonlocalA 0.48); local: the same without nonlocalA.
#
# Needs ./Allwmake and the test driver: wmake tests/void/voidCheck.
# Work directory: $CEL_WORK (default /tmp/cel_final); the report is
# $CEL_WORK/<model>.out (hopper_pressure_checks.py output).
#
# Checkpoints (decision 0117): a restart time every $CEL_CHECKPOINT s of
# simulated time (default 0.05). Run the script again after a lost
# container and it resumes from the latest checkpoint. While it runs, the
# case (latest two checkpoints only) is mirrored every 5 min to
# $CEL_MIRROR (default /mnt/project-files/hopper/final4mm/ckpt) and is
# restored from there when the work directory has gone.
source /usr/lib/openfoam/openfoam2512/etc/bashrc
cd "$(dirname "$0")/../.." || exit 1
model=$1
level=$2
law=$3
case "$level" in
    8)   ncols=29 htop=0.024 hhop=0.008 ;;
    5.7) ncols=41 htop=0.017 hhop=0.0057 ;;
    4)   ncols=58 htop=0.012 hhop=0.004 ;;
    *) echo "level must be 8, 5.7 or 4"; exit 2 ;;
esac
case "$law" in
    nl) nl='nonlocalA 0.48;\n' ;;
    local) nl='' ;;
    *) echo "law must be nl or local"; exit 2 ;;
esac
run=${model}_${level}mm_$law
work=${CEL_WORK:-/tmp/cel_refine}/$run
mirror=${CEL_MIRROR:-/mnt/project-files/hopper/nonlocal-refine/ckpt}/$run
case=$work/$model/couto_${model}_aligned_buffered
mkdir -p "$work/$model" "$mirror"

m=$mirror/couto_${model}_aligned_buffered
[ -d "$m" ] || { [ -d "$m.old" ] && mv "$m.old" "$m"; }
if [ ! -f "$case/history.dat" ] && [ -f "$m/history.dat" ]; then
    echo "restoring $case from $mirror"
    cp -a "$mirror/couto_${model}_aligned_buffered" "$work/$model/"
fi

# keep the latest two checkpoints, mirror the case
sync_case() {
    [ -d "$case/processor0" ] || return
    for pd in "$case"/processor*; do
        (cd "$pd" && ls -d [0-9]* 2>/dev/null | grep -vx 0 | sort -g | head -n -2 \
            | xargs -r rm -rf)
    done
    local name dst
    name=$(basename "$case")
    dst=$mirror/$name
    mkdir -p "$mirror"
    rm -rf "$dst.new" "$dst.old"
    cp -a "$case" "$dst.new" || return
    [ -d "$dst" ] && mv "$dst" "$dst.old"
    mv "$dst.new" "$dst" && rm -rf "$dst.old"
}
( while true; do sleep 300; sync_case; done ) &
syncer=$!

# mesh: the decision 0117 family (wall column a quarter of the axis column,
# outlet columns even), scaled by the level
export CEL_NPROCS=4 CEL_H_TOP=$htop CEL_NCOLS=$ncols CEL_WALL_RATIO=0.25 \
       CEL_H_HOP=$hhop CEL_OUTLET_WALL_RATIO=1.0
export CEL_DP_FIT=plane_strain
export CEL_CHECKPOINT=${CEL_CHECKPOINT:-0.05} CEL_RESUME=1
# T0 rules + fullBlend + nonlocal fluidity (decisions 0102, 0105, 0111)
export CEL_HOPPER_RULES="$(printf 'mu2 0.8412;\nI0 0.279;\nd 0.004;\nrhoS 1380;\nshearCap yes;\nspreadLayers yes;\nflowImpedance yes;\nnormalCap yes;\nfullBlend yes;\n'"$nl")"
python3 tests/hopper/hopper_pressure_checks.py "$work/$model" "$model" aligned_buffered \
    > "$work/$model.out" 2>&1
status=$?
kill $syncer 2>/dev/null
sync_case
cp "$work/$model.out" "$mirror/" 2>/dev/null
echo "EXIT $status" >> "$work/$model.out"
