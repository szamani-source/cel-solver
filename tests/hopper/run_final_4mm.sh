#!/bin/bash
# Final validation of the Couto silo on the fine 4 mm mesh (decision 0116):
# gravity-on start, outlet opens at 0.5 s, run to 1.5 s, 4 MPI ranks.
#
#     tests/hopper/run_final_4mm.sh mc|dp      (about 6 h each on 4 cores)
#
# Needs ./Allwmake and the test driver: wmake tests/void/voidCheck.
# Work directory: $CEL_WORK (default /tmp/cel_final); the report is
# $CEL_WORK/<model>.out (hopper_pressure_checks.py output).
#
# Checkpoints (decision 0116): a restart time every $CEL_CHECKPOINT s of
# simulated time (default 0.05). Run the script again after a lost
# container and it resumes from the latest checkpoint. While it runs, the
# case (latest two checkpoints only) is mirrored every 5 min to
# $CEL_MIRROR (default /mnt/project-files/hopper/final4mm/ckpt) and is
# restored from there when the work directory has gone.
source /usr/lib/openfoam/openfoam2512/etc/bashrc
cd "$(dirname "$0")/../.." || exit 1
model=$1
work=${CEL_WORK:-/tmp/cel_final}
mirror=${CEL_MIRROR:-/mnt/project-files/hopper/final4mm/ckpt}/$model
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

# mesh: 58 columns (16 mm at the axis to 4 mm at the wall), 4 mm rows in the
# hopper, 12 mm rows above, outlet columns even (3.1-3.2 mm)
export CEL_NPROCS=4 CEL_H_TOP=0.012 CEL_NCOLS=58 CEL_WALL_RATIO=0.25 \
       CEL_H_HOP=0.004 CEL_OUTLET_WALL_RATIO=1.0
export CEL_DP_FIT=plane_strain
export CEL_CHECKPOINT=${CEL_CHECKPOINT:-0.05} CEL_RESUME=1
# T0 rules + fullBlend + nonlocal fluidity (decisions 0102, 0105, 0111)
export CEL_HOPPER_RULES="$(printf 'mu2 0.8412;\nI0 0.279;\nd 0.004;\nrhoS 1380;\nshearCap yes;\nspreadLayers yes;\nflowImpedance yes;\nnormalCap yes;\nfullBlend yes;\nnonlocalA 0.48;\n')"
python3 tests/hopper/hopper_pressure_checks.py "$work/$model" "$model" aligned_buffered \
    > "$work/$model.out" 2>&1
status=$?
kill $syncer 2>/dev/null
sync_case
cp "$work/$model.out" "$mirror/" 2>/dev/null
echo "EXIT $status" >> "$work/$model.out"
