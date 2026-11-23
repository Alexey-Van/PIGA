#!/usr/bin/env python3
"""Repair PEPPER gVCF PL vectors that omit the trailing <*>/<*> likelihood.

For a diploid record with REF plus two ALT alleles (the second ALT being the
gVCF symbolic <*> allele), VCF Number=G requires six PL values. Some
PEPPER-Margin-DeepVariant outputs contain only the first five. This script
appends the largest observed PL value as the least-supported likelihood. It
does not alter any other field or silently repair any other malformed shape.

Input may be plain or gzip/BGZF VCF. The script writes plain VCF to stdout so
the caller can BGZF-compress it and build a fresh index.
"""

import argparse
import gzip
import math
import sys


def genotype_count(alleles: int, ploidy: int) -> int:
    """Number of unordered genotypes for allele count and ploidy."""
    return math.comb(alleles + ploidy - 1, ploidy)


def open_text(path: str):
    """Open plain or gzip-compressed VCF as text."""
    with open(path, "rb") as handle:
        compressed = handle.read(2) == b"\x1f\x8b"
    if compressed:
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "rt", encoding="utf-8")


def fix_record(line: str, repaired: list[int]) -> str:
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 10:
        return line

    chrom, pos, alt_text = fields[0], fields[1], fields[4]
    alts = alt_text.split(",")
    if not alts or alts[-1] != "<*>":
        return line

    format_keys = fields[8].split(":")
    if "PL" not in format_keys:
        return line
    pl_index = format_keys.index("PL")
    gt_index = format_keys.index("GT") if "GT" in format_keys else None
    expected = None

    for sample_index, sample_value in enumerate(fields[9:]):
        sample_parts = sample_value.split(":")
        if pl_index >= len(sample_parts) or sample_parts[pl_index] in ("", "."):
            continue
        pl_values = sample_parts[pl_index].split(",")

        if gt_index is None or gt_index >= len(sample_parts):
            raise ValueError(f"{chrom}:{pos}: PL exists but GT is unavailable")
        gt = sample_parts[gt_index]
        ploidy = len(gt.replace("|", "/").split("/")) if gt not in ("", ".") else 1
        expected = genotype_count(1 + len(alts), ploidy)

        if len(pl_values) == expected:
            continue

        # Supported producer shape: diploid REF + real ALT + <*> with the
        # final (<*>/<*>) likelihood omitted. Preserve known values and use
        # their maximum as the conservative, least-supported likelihood.
        if len(alts) == 2 and ploidy == 2 and len(pl_values) == expected - 1:
            try:
                trailing_value = max(int(value) for value in pl_values)
            except ValueError as exc:
                raise ValueError(
                    f"{chrom}:{pos} sample column {sample_index + 10}: "
                    f"cannot infer trailing PL from {sample_parts[pl_index]!r}"
                ) from exc
            sample_parts[pl_index] = ",".join((*pl_values, str(trailing_value)))
            fields[sample_index + 9] = ":".join(sample_parts)
            repaired[0] += 1
            continue

        raise ValueError(
            f"{chrom}:{pos} sample column {sample_index + 10}: "
            f"PL has {len(pl_values)} values; expected {expected}; "
            f"ALT={alt_text}, GT={gt}"
        )

    return "\t".join(fields) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gvcf", help="PEPPER gVCF, plain or gzip/BGZF compressed")
    args = parser.parse_args()

    repaired = [0]
    try:
        with open_text(args.gvcf) as source:
            for line in source:
                if line.startswith("#"):
                    sys.stdout.write(line)
                else:
                    sys.stdout.write(fix_record(line, repaired))
    except (OSError, UnicodeError, ValueError) as exc:
        print(f"gVCF PL compatibility error: {exc}", file=sys.stderr)
        return 1

    print(
        f"Normalized {repaired[0]} trailing <*> PL value(s) from {args.gvcf}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
