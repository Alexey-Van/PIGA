# Draft diploid genome assembly

This pipeline performs reference-guided draft diploid genome assembly.

<img align="middle" width="942" src="draft_assembly.jpg"/>

## draft_assembly

### Description

- Map Illumina reads to each sample's personal reference and call short-read variants.
- Combine those calls with lifted, population-phased variants, phase with ONT reads, and correct switch errors against the population phasing.
- Haplotag ONT reads with the corrected phased VCF and split the original FASTQ reads by haplotype, retaining their basecaller metadata.
- Assemble each haplotype with hifiasm in ONT mode and polish it with Medaka.
- For male samples, assign otherwise untagged chrX non-PAR reads to the X haplotype and chrY reads to the Y haplotype; PAR reads remain assigned only when supported by haplotags.

This replaces the original PacBio ZMW/subread interval assembly and Arrow polishing with ONT read phasing, hifiasm, and Medaka. Short reads still contribute to variant calling, switch-error correction, and conservative homozygous-error polishing. The final diploid-path stage uses ONT reads for structural-variant validation and the Illumina-only k-mer database for Merqury clipping.

### Requirement

- [BWA](https://github.com/lh3/bwa)
- [GATK](https://github.com/broadinstitute/gatk)
- [minimap2](https://github.com/lh3/minimap2)
- [WhatsHap](https://github.com/whatshap/whatshap)
- [Samtools](https://github.com/samtools/samtools)
- [Bcftools](https://github.com/samtools/bcftools)
- [Seqkit](https://github.com/shenwei356/seqkit)
- [hifiasm](https://github.com/chhylp123/hifiasm)
- [Medaka](https://github.com/nanoporetech/medaka)
- [Freebayes](https://github.com/freebayes/freebayes)
- [DeepVariant](https://github.com/google/deepvariant)
- pysam
- [CrossMap](https://github.com/liguowang/CrossMap)

### Configuration (`draft_assembly.yaml`)

The configuration file should contain:

`samples`: Specify the path to a text file listing all samples.
The file must be space-delimited, with:

- Column 1: Sample name
- Column 2: Sample sex

`sr_fastqs`: Specify the paths to the paired-end short-read FASTQ files.
Use `{sample}` as a wildcard; it will be automatically replaced with sample names from the `samples` file.

`lr_fastqs`: Specify the path(s) to ONT long-read FASTQ files.
Use `{sample}` as a wildcard.

`medaka_model`: Medaka consensus model matching the flow cell chemistry and basecaller. `auto` selects from read metadata; use an explicit model name if that metadata is missing.

`reference`:

- `CHM13`: Provide the path to the T2T CHM13 human reference genome (FASTA format, indexed by `samtools faidx`).

`par_region`: Path to the pseudoautosomal region (PAR) BED file for CHM13.
This is used to correctly handle sex-chromosome reconstruction.

`prefix`: Prefix used for naming output files.

`consensus_fasta` _(optional)_: Path to the sample-specific consensus FASTA file.
If omitted, the workflow will use the default:
`c5_personal_ref/consensus_fasta/{sample}/{sample}.personal_ref.fasta`

`chain` _(optional)_: Path to the chain file describing coordinate mapping between CHM13 and the consensus assembly.
If omitted, the workflow will use the default:
`c5_personal_ref/consensus_fasta/{sample}/{sample}.personal_ref.chain`

`sample_phased_vcf` _(optional)_: Path to the sample's phased VCF.
If omitted, the workflow will use the default:
`c4_phase_snv/sample_vcf/{sample}/{sample}.shapeit.vcf.gz`

### Usage

```bash
snakemake -s Snakefile --cores 64 --jobs 64 --configfile config/draft_assembly.yaml --profile ./profile/config_local/
```

### Output
The output file should contain:
- `c6_draft_assembly/sample_assembly/{sample}/assembly/{sample}.{hap1,hap2}.fasta`: Draft diploid assembly of each sample
