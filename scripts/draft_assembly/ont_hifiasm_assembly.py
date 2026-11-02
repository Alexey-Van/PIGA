#!/usr/bin/env python3
"""Haplotype-specific ONT draft assemblies with hifiasm and Medaka."""
import argparse
import os
import subprocess


def run(command, stdin=None, stdout=None):
    subprocess.run(command, check=True, stdin=stdin, stdout=stdout)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", required=True)
    parser.add_argument("--variants", required=True)
    parser.add_argument("--reads", required=True)
    parser.add_argument("--medaka-model", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--threads", required=True, type=int)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    tagged_bam = os.path.join(args.output_dir, f"{args.sample}.haplotagged.bam")
    raw_sam = os.path.join(args.output_dir, f"{args.sample}.ont.sam")
    with open(raw_sam, "wb") as sam_out:
        run(["minimap2", "-t", str(args.threads), "-ax", "map-ont", args.reference, args.reads],
            stdout=sam_out)
    sorted_bam = os.path.join(args.output_dir, f"{args.sample}.ont.srt.bam")
    with open(raw_sam, "rb") as sam:
        run(["samtools", "sort", "-@", str(args.threads), "-o", sorted_bam, "-"], stdin=sam)
    os.remove(raw_sam)
    run(["samtools", "index", sorted_bam])
    run(["whatshap", "haplotag", "--ignore-read-groups", "--reference", args.reference,
         "--output", tagged_bam, args.variants, sorted_bam])

    for hap in (1, 2):
        hap_reads = os.path.join(args.output_dir, f"{args.sample}.hap{hap}.fastq")
        with open(hap_reads, "wb") as fq:
            view = subprocess.Popen(["samtools", "view", "-bh", "-d", f"HP:{hap}", tagged_bam], stdout=subprocess.PIPE)
            subprocess.run(["samtools", "fastq", "-"], stdin=view.stdout, stdout=fq, check=True)
            view.stdout.close()
            if view.wait() != 0:
                raise subprocess.CalledProcessError(view.returncode, view.args)
        if os.path.getsize(hap_reads) == 0:
            raise RuntimeError(f"WhatsHap assigned no reads to haplotype {hap}; cannot assemble it")

        prefix = os.path.join(args.output_dir, f"{args.sample}.hap{hap}")
        run(["hifiasm", "--ont", "-t", str(args.threads), "-o", prefix, hap_reads])
        gfa = prefix + ".bp.p_ctg.gfa"
        draft = prefix + ".draft.fa"
        with open(gfa) as source, open(draft, "w") as target:
            for line in source:
                fields = line.rstrip("\n").split("\t")
                if fields[0] == "S":
                    target.write(f">{fields[1]}\n{fields[2]}\n")

        medaka_dir = prefix + ".medaka"
        medaka = ["medaka_consensus", "-i", hap_reads, "-d", draft, "-o", medaka_dir,
                  "-t", str(args.threads)]
        if args.medaka_model.lower() == "auto":
            medaka.extend(["--auto_model", "consensus", hap_reads])
        else:
            medaka.extend(["-m", args.medaka_model])
        run(medaka)
        consensus = os.path.join(medaka_dir, "consensus.fasta")
        destination = os.path.join(args.output_dir, f"{args.sample}.hap{hap}.fasta")
        with open(consensus) as source, open(destination, "w") as target:
            for line in source:
                if line.startswith(">"):
                    target.write(f">hap{hap}_{line[1:]}" if not line.startswith(f">hap{hap}_") else line)
                else:
                    target.write(line)


if __name__ == "__main__":
    main()
