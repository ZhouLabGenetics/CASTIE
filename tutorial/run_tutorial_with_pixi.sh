#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${REPO_ROOT}"

pixi install
pixi run build
pixi run python tutorial/generate_simulated_data.py

CASTIE_SOURCE_ROOT="${REPO_ROOT}" \
  pixi run bash tutorial/run_tutorial_in_container.sh
