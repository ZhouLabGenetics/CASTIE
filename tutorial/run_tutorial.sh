#!/usr/bin/env bash
set -euo pipefail

TUTORIAL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
IMAGE="${CASTIE_IMAGE:-yijia0802/castie:Latest}"

python3 "${TUTORIAL_DIR}/generate_simulated_data.py"
mkdir -p "${TUTORIAL_DIR}/output"

docker run --rm --platform linux/amd64 \
  -v "${TUTORIAL_DIR}:/tutorial" \
  -w /tutorial \
  "${IMAGE}" \
  bash run_tutorial_in_container.sh

echo "Tutorial complete. Results are in ${TUTORIAL_DIR}/output"
