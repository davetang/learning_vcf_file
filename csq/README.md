## README

Problem:

> However, current predictors analyse variants as isolated events, which can
lead to incorrect predictions when adjacent variants alter the same codon, or
when a frame-shifting indel is followed by a frame-restoring indel.

BCFtools/csq is a fast program for haplotype-aware consequence calling which
can take into account known phase.

There are several popular existing programs for variant annotation including:

1. Ensembl Variant Effect Predictor (VEP)
2. SnpEff
3. ANNOVAR

but they do not take phasing into account.

## BCFtools/csq

`bcftools csq` requires a phased VCF, a GFF3 file with gene predictions, and a
reference FASTA file.

```
About: Haplotype-aware consequence caller.
Usage: bcftools csq [OPTIONS] in.vcf

Required options:
   -f, --fasta-ref FILE              Reference file in fasta format
   -g, --gff-annot FILE              GFF3 annotation file

CSQ options:
   -B, --trim-protein-seq INT        Abbreviate protein-changing predictions to max INT aminoacids
   -c, --custom-tag STRING           Use this tag instead of the default BCSQ
   -l, --local-csq                   Localized predictions, consider only one VCF record at a time
   -n, --ncsq INT                    Maximum number of per-haplotype consequences to consider for each site [15]
   -p, --phase a|m|r|R|s             How to handle unphased heterozygous genotypes: [r]
                                       a: take GTs as is, create haplotypes regardless of phase (0/1 -> 0|1)
                                       m: merge *all* GTs into a single haplotype (0/1 -> 1, 1/2 -> 1)
                                       r: require phased GTs, throw an error on unphased het GTs
                                       R: create non-reference haplotypes if possible (0/1 -> 1|1, 1/2 -> 1|2)
                                       s: skip unphased hets
Options:
   -e, --exclude EXPR                Exclude sites for which the expression is true
       --force                       Run even if some sanity checks fail
   -i, --include EXPR                Select sites for which the expression is true
       --no-version                  Do not append version and command line to the header
   -o, --output FILE                 Write output to a file [standard output]
   -O, --output-type b|u|z|v|t[0-9]  b: compressed BCF, u: uncompressed BCF, z: compressed VCF
                                     v: uncompressed VCF, t: plain tab-delimited text output, 0-9: compression level [v]
   -r, --regions REGION              Restrict to comma-separated list of regions
   -R, --regions-file FILE           Restrict to regions listed in a file
       --regions-overlap 0|1|2       Include if POS in the region (0), record overlaps (1), variant overlaps (2) [1]
   -s, --samples -|LIST              Samples to include or "-" to apply all variants and ignore samples
   -S, --samples-file FILE           Samples to include
   -t, --targets REGION              Similar to -r but streams rather than index-jumps
   -T, --targets-file FILE           Similar to -R but streams rather than index-jumps
       --targets-overlap 0|1|2       Include if POS in the region (0), record overlaps (1), variant overlaps (2) [0]
       --threads INT                 Use multithreading with <int> worker threads [0]
   -v, --verbose INT                 Verbosity level 0-2 [1]

Example:
   bcftools csq -f hs37d5.fa -g Homo_sapiens.GRCh37.82.gff3.gz in.vcf

   # GFF3 annotation files can be downloaded from Ensembl. e.g. for human:
   ftp://ftp.ensembl.org/pub/current_gff3/homo_sapiens/
   ftp://ftp.ensembl.org/pub/grch37/release-84/gff3/homo_sapiens/
```

The program begins by parsing gene predictions in the GFF3 file, then streams
through the VCF file using a fast region lookup at each site to find overlaps
with regions of supported genomic types (exons, CDS, UTRs or general
transcripts). For more details read the paper (see [Further
reading](#further-reading).

## Further reading

* [BCFtools/csq: haplotype-aware variant consequences
](https://academic.oup.com/bioinformatics/article/33/13/2037/3000373)
