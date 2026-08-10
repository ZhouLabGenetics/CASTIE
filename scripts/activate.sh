#!/usr/bin/env bash

# Make CASTIE's source-tree command-line scripts available in `pixi shell`
# and to commands launched with `pixi run`.
export PATH="${PIXI_PROJECT_ROOT}/extdata:${PATH}"
