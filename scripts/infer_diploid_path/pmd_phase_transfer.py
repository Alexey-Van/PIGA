#!/usr/bin/env python3
"""Transfer PEPPER-Margin phase sets to matching consensus genotypes."""
import argparse
import pysam


def allele_key(record):
    return (record.contig, record.pos, record.ref, tuple(record.alts or ()))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-vcf", required=True)
    parser.add_argument("--phased-vcf", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--output-vcf", required=True)
    args = parser.parse_args()

    with pysam.VariantFile(args.phased_vcf) as phased_file:
        if args.sample not in phased_file.header.samples:
            raise ValueError(f"Sample {args.sample} missing from PEPPER phased VCF")
        phases = {}
        for record in phased_file:
            call = record.samples[args.sample]
            gt = call.get("GT")
            ps = call.get("PS")
            if gt and len(gt) == 2 and call.phased and ps is not None:
                phases[allele_key(record)] = (tuple(gt), ps)

    # Open the source without requiring an index: this script streams through
    # every record, so random access is not needed and Snakemake may not have
    # created a .tbi for this intermediate VCF.
    with pysam.VariantFile(args.input_vcf, "r") as source:
        if args.sample not in source.header.samples:
            raise ValueError(f"Sample {args.sample} missing from consensus VCF")
        header = source.header.copy()
        # PS is required before mutating records. Header-format membership can
        # differ between pysam/htslib versions, so inspect the actual format
        # definitions and add it if absent.
        if not any(line.key == "FORMAT" and line.get("ID") == "PS" for line in header.records):
            header.formats.add("PS", 1, "Integer", "Phase set from PEPPER-Margin")
        with pysam.VariantFile(args.output_vcf, "wz", header=header) as out:
            matched = 0
            for record in source:
                # Records are associated with the source header. Translate
                # them to the output header so newly declared FORMAT fields
                # (notably PS) are recognized when sample values are set.
                record.translate(header)
                match = phases.get(allele_key(record))
                if match:
                    gt, ps = match
                    call = record.samples[args.sample]
                    current_gt = call.get("GT")
                    if current_gt and len(current_gt) == 2 and tuple(sorted(current_gt)) == tuple(sorted(gt)):
                        call["GT"] = gt
                        call["PS"] = ps
                        call.phased = True
                        matched += 1
                out.write(record)
    print(f"Transferred PEPPER phase for {matched} matching genotypes", flush=True)


if __name__ == "__main__":
    main()
