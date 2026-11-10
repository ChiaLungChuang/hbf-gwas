# resources/

## long_range_ld_grch37.tsv

Regions of long-range, high linkage disequilibrium excluded before PCA
(`scripts/03_pca.sh`). GRCh37, 1-based inclusive coordinates, 24 autosomal
regions. Columns: `chrom start end id verified note`.

**Source.** Transcribed on 2026-11-10 from the GRCh37 table of the Genome
Analysis Wiki page "Regions of high linkage disequilibrium (LD)", revision
15179: <https://genome.sph.umich.edu/wiki/Regions_of_high_linkage_disequilibrium_(LD)>.
That page attributes the regions to Price et al. (2008), "Long-Range LD Can
Confound Genome Scans in Admixed Populations", *Am. J. Hum. Genet.*, and to
Weale (2010), "Quality Control for Genome-Wide Association Studies", in
*Genetic Variation: Methods and Protocols*, Methods Mol. Biol. 628,
doi:10.1007/978-1-60327-367-1_19. The page prints the Price et al. volume and
pages as "86, 127-147"; that citation was not checked here.

**Every entry is `verified = false`.** The wiki's GRCh37 coordinates are
rounded and some differ from its own build 36 table by more than a lift-over
would explain (e.g. chr6 25-35 Mb vs 25.5-33.5 Mb in build 36). Nothing here
was compared with the primary paper. The chrX regions on the wiki's build 36
table are not included: the pipeline is autosomal only.

The ids `lrld01`-`lrld24` were assigned here; they are not from the source.
