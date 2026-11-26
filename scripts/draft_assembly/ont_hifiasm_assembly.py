#!/usr/bin/env python3
"""Assemble and polish two haplotypes from WhatsHap-tagged ONT reads."""
import argparse
import os
import subprocess

import pysam


def par_intervals(bed_path):
    intervals = []
    with open(bed_path) as bed:
        for line in bed:
            if line.strip() and not line.startswith("#"):
                chrom, start, end = line.split()[:3]
                intervals.append((chrom, int(start), int(end)))
    return intervals


def tagged_read_names(bam_path, sex, par_bed):
    names = {1: set(), 2: set()}
    par = par_intervals(par_bed)
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            haplotype = read.get_tag("HP") if read.has_tag("HP") else None
            if haplotype is None and sex == "male":
                if read.reference_name == "chrY":
                    haplotype = 1
                elif read.reference_name == "chrX":
                    start, end = read.reference_start, read.reference_end
                    par_overlap = sum(
                        max(0, min(end, par_end) - max(start, par_start))
                        for chrom, par_start, par_end in par
                        if chrom == "chrX"
                    )
                    if par_overlap < end - start:
                        haplotype = 2
            if haplotype in names:
                names[haplotype].add(read.query_name)
    return names


def write_haplotype_fastqs(read_fastqs, read_names, output_dir):
    paths = {
        haplotype: os.path.join(output_dir, f"hap{haplotype}.fastq")
        for haplotype in (1, 2)
    }
    found = {1: set(), 2: set()}

    with open(paths[1], "w") as hap1, open(paths[2], "w") as hap2:
        outputs = {1: hap1, 2: hap2}
        for fastq_path in read_fastqs:
            with pysam.FastxFile(fastq_path) as reads:
                for read in reads:
                    haplotype = None
                    for candidate in (1, 2):
                        if read.name in read_names[candidate]:
                            haplotype = candidate
                            break
                    if haplotype is None:
                        continue
                    if read.quality is None:
                        raise ValueError(f"Input file is not FASTQ: {fastq_path}")
                    comment = f" {read.comment}" if read.comment else ""
                    outputs[haplotype].write(
                        f"@{read.name}{comment}\n{read.sequence}\n+\n{read.quality}\n"
                    )
                    found[haplotype].add(read.name)

    for haplotype in (1, 2):
        if not found[haplotype]:
            raise RuntimeError(f"No source FASTQ reads found for haplotype {haplotype}")
        missing = read_names[haplotype] - found[haplotype]
        if missing:
            raise RuntimeError(
                f"{len(missing)} haplotype {haplotype} BAM reads are absent from source FASTQ"
            )
    return paths


def run(command):
    subprocess.run(command, check=True)


def gfa_to_fasta(gfa_path, fasta_path):
    with open(gfa_path) as source, open(fasta_path, "w") as target:
        for line in source:
            fields = line.rstrip("\n").split("\t")
            if fields[0] == "S":
                if fields[2] == "*":
                    raise RuntimeError(f"Sequence missing for segment {fields[1]} in {gfa_path}")
                target.write(f">{fields[1]}\n{fields[2]}\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tagged-bam", required=True)
    parser.add_argument("--read-fastqs", nargs="+", required=True)
    parser.add_argument("--sex", choices=("female", "male"), required=True)
    parser.add_argument("--par-bed", required=True)
    parser.add_argument("--medaka-model", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--threads", required=True, type=int)
    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    read_names = tagged_read_names(args.tagged_bam, args.sex, args.par_bed)
    haplotype_fastqs = write_haplotype_fastqs(
        args.read_fastqs, read_names, args.output_dir
    )

    for haplotype in (1, 2):
        reads = haplotype_fastqs[haplotype]
        prefix = os.path.join(args.output_dir, f"{args.sample}.hap{haplotype}")
        run(["hifiasm", "--ont", "-i", "-t", str(args.threads), "-o", prefix, reads])

        draft = prefix + ".draft.fa"
        gfa_to_fasta(prefix + ".bp.p_ctg.gfa", draft)
        medaka_dir = prefix + ".medaka"
        medaka = [
            "medaka_consensus",
            "-i",
            reads,
            "-d",
            draft,
            "-o",
            medaka_dir,
            "-t",
            str(args.threads),
            "-f",
        ]
        # medaka_consensus accepts an explicit model via -m; --auto_model is
        # not a medaka_consensus option (model resolution is a separate tool).
        medaka.extend(["-m", args.medaka_model])
        run(medaka)

        output_fasta = os.path.join(
            args.output_dir, f"{args.sample}.hap{haplotype}.medaka.fasta"
        )
        with open(os.path.join(medaka_dir, "consensus.fasta")) as source, open(
            output_fasta, "w"
        ) as target:
            for line in source:
                if line.startswith(">"):
                    target.write(f">hap{haplotype}_{line[1:]}")
                else:
                    target.write(line)


if __name__ == "__main__":
    main()
