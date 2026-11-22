#!/bin/bash
# Final validation of the Couto silo on the fine 4 mm mesh (decision 0116):
# gravity-on start, outlet opens at 0.5 s, run to 1.5 s, 4 MPI ranks.
#
#     tests/hopper/run_final_4mm.sh mc|dp      (about 6 h each on 4 cores)
#
# Needs ./Allwmake and the test driver: wmake tests/void/voidCheck.
# Work directory: $CEL_WORK (default /tmp/cel_final); the report is
# $CEL_WORK/<model>.out (hopper_pressure_checks.py output).
source /usr/lib/openfoam/openfoam2512/etc/bashrc
cd "$(dirname "$0")/../.." || exit 1
work=${CEL_WORK:-/tmp/cel_final}
mkdir -p "$work"
# mesh: 58 columns (16 mm at the axis to 4 mm at the wall), 4 mm rows in the
# hopper, 12 mm rows above, outlet columns even (3.1-3.2 mm)
export CEL_NPROCS=4 CEL_H_TOP=0.012 CEL_NCOLS=58 CEL_WALL_RATIO=0.25 \
       CEL_H_HOP=0.004 CEL_OUTLET_WALL_RATIO=1.0
export CEL_DP_FIT=plane_strain
# T0 rules + fullBlend + nonlocal fluidity (decisions 0102, 0105, 0111)
export CEL_HOPPER_RULES="$(printf 'mu2 0.8412;\nI0 0.279;\nd 0.004;\nrhoS 1380;\nshearCap yes;\nspreadLayers yes;\nflowImpedance yes;\nnormalCap yes;\nfullBlend yes;\nnonlocalA 0.48;\n')"
python3 tests/hopper/hopper_pressure_checks.py "$work/$1" "$1" aligned_buffered \
    > "$work/$1.out" 2>&1
echo "EXIT $?" >> "$work/$1.out"
