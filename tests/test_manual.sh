#!/bin/bash
# Series M (user manual, decision 0116): the HTML pages in docs/manual are consistent.
# Checks navigation, links, citations, equations (typeset with the vendored KaTeX when
# node is installed), each page's required content and, once the case setup page
# exists, that every input keyword the code reads is documented. No OpenFOAM needed.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
python3 "$root/tools/check_manual.py"
