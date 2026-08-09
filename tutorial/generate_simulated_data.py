#!/usr/bin/env python3
"""Generate a small deterministic CASTIE tutorial dataset.

The script uses only the Python standard library and writes PLINK 1 binary
files plus a cell-level phenotype table. The data are synthetic and intended
only to demonstrate the CASTIE workflow.
"""

from __future__ import annotations

import math
import random
import shutil
from pathlib import Path


SEED = 20260809
N_DONORS = 40
CELLS_PER_DONOR = 5
N_VARIANTS = 200


def poisson(rng: random.Random, rate: float) -> int:
    threshold = math.exp(-rate)
    product = 1.0
    value = 0
    while product > threshold:
        value += 1
        product *= rng.random()
    return value - 1


def write_plink(prefix: Path, genotypes: list[list[int]]) -> None:
    with prefix.with_suffix(".fam").open("w", encoding="utf-8") as handle:
        for donor in range(1, N_DONORS + 1):
            sample = f"donor_{donor:03d}"
            handle.write(f"{sample}\t{sample}\t0\t0\t0\t-9\n")

    with prefix.with_suffix(".bim").open("w", encoding="utf-8") as handle:
        for variant in range(N_VARIANTS):
            position = 100_000 + variant * 1_000
            handle.write(f"1\trs{variant + 1}\t0\t{position}\tA\tG\n")

    # PLINK .bed SNP-major encoding: 00=A1/A1, 10=heterozygous, 11=A2/A2.
    codes = {0: 0b00, 1: 0b10, 2: 0b11}
    with prefix.with_suffix(".bed").open("wb") as handle:
        handle.write(bytes((0x6C, 0x1B, 0x01)))
        for variant_genotypes in genotypes:
            for start in range(0, N_DONORS, 4):
                packed = 0
                for offset, genotype in enumerate(variant_genotypes[start : start + 4]):
                    packed |= codes[genotype] << (2 * offset)
                handle.write(bytes((packed,)))


def main() -> None:
    rng = random.Random(SEED)
    tutorial_dir = Path(__file__).resolve().parent
    data_dir = tutorial_dir / "data"
    if data_dir.exists():
        shutil.rmtree(data_dir)
    data_dir.mkdir()

    donor_x1 = [rng.gauss(0, 1) for _ in range(N_DONORS)]
    donor_x2 = [donor % 2 for donor in range(N_DONORS)]

    genotypes: list[list[int]] = []
    for variant in range(N_VARIANTS):
        maf = 0.10 + 0.05 * (variant % 7)
        values = [int(rng.random() < maf) + int(rng.random() < maf) for _ in range(N_DONORS)]
        if len(set(values)) == 1:
            values[0] = 1
        genotypes.append(values)

    write_plink(data_dir / "genotypes", genotypes)
    shutil.copyfile(data_dir / "genotypes.bed", data_dir / "grm_variants.bed")
    shutil.copyfile(data_dir / "genotypes.bim", data_dir / "grm_variants.bim")
    shutil.copyfile(data_dir / "genotypes.fam", data_dir / "grm_variants.fam")

    with (data_dir / "phenotypes.tsv").open("w", encoding="utf-8") as handle:
        handle.write("IND_ID\tX1\tX2\tpf1\tpf2\tgene_1\n")
        for donor in range(N_DONORS):
            causal_genotype = genotypes[9][donor]
            for _ in range(CELLS_PER_DONOR):
                pf1 = rng.gauss(0, 1)
                pf2 = rng.gauss(0, 1)
                log_rate = 1.0 + 0.20 * donor_x1[donor] + 0.15 * pf1 + 0.25 * causal_genotype
                count = poisson(rng, math.exp(log_rate))
                handle.write(
                    f"donor_{donor + 1:03d}\t{donor_x1[donor]:.6f}\t{donor_x2[donor]}"
                    f"\t{pf1:.6f}\t{pf2:.6f}\t{count}\n"
                )

    with (data_dir / "gene_1_region.tsv").open("w", encoding="utf-8") as handle:
        handle.write("1\t95000\t305000\n")

    print(f"Wrote simulated data to {data_dir}")


if __name__ == "__main__":
    main()
