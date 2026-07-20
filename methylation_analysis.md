Methylation bioinformatic pipeline
================
Jacob Cassens
2026-01-05

- [Methylation](#methylation)

### Methylation

#### Concatenate fastq files

As a pre-processing step, all fastq files generated from the sequencing
experiment were concatenated into a single fastq file for each
individual. Downstream analyses use the concatenated fastq file.

``` bash
cat *.fastq.gz > iscap_MN_meth.fastq.gz
```

#### Mapping reads to reference genome

Raw reads were mapped to the [PalLabHifi reference
genome](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_016920785.2/)
using [Minimap2](https://github.com/lh3/minimap2). Mapped sam files were
converted to bam files, sorted, and indexed using
[samtools](http://www.htslib.org/).

``` bash
# Align to reference
minimap2 -ax map-ont pal_reference.fasta iscap_MN_meth.fastq.gz -t 32 > iscap_MN_align_meth_pal.sam
# Convert to sam to bam file 
samtools view -S -b iscap_MN_align_meth_pal.sam > iscap_MN_meth_mapped_pal.bam
# Sort mapped bam file to removed unmapped reads
samtools view -b -F 4 iscap_MN_meth_mapped_pal.bam > iscap_MN_meth_aligned_pal.bam
```

#### Converting to bedMethyl format

The mapped BAM file containing modified base information was converted
to bedMethyl format using
[Modkit](https://github.com/nanoporetech/modkit) to extract methylation
counts.

``` bash
modkit pileup --log-filepath modkit.log -t 28 -r pal_reference.fasta --cpg --combine-strands --only-tabs --prefix iscap_MN-HP --partition-tag HP iscap_MN_meth_aligned_pal.bam modkit_longphase
```

modkit output contains both 5mC and 5hmC. We extracted 5mC counts with
awk -v OFS=“ ‘\$4 ==“m”’ HP.bed.

#### Methylkit

##### Filtering and normalization

The following steps provide an example for one comparison group (MN and
PA versus TX), although this analysis was repeated for three other
comparison groups.

``` r
# Load packages
library(tidyverse)
library(methylKit)
## 5mC analysis comparing Minnesota and Pennsylvania individuals to the Texas individual
# List the modkit files
file.list.IS <- list("/Volumes/Xtreme/aim2/methylation/input_files/IS_MN_modkit_5mc.bed.gz",
                      "/Volumes/Xtreme/aim2/methylation/input_files/IS_PA_modkit_5mc.bed.gz",
                      "/Volumes/Xtreme/aim2/methylation/input_files/IS_TX_modkit_5mc.bed.gz")
sample.id.IS <- list("IS_MN", "IS_PA", "IS_TX")
# Specify the modkit column IDs
modkit_cols <- list(
  fraction = FALSE,
  chr.col = 1,
  start.col = 2,
  end.col = 3,
  coverage.col = 5,
  strand.col = 6,
  freqC.col = 11
)
# Read in methylation data
methylationData.IS <- methRead(file.list.IS,
                                sample.id = sample.id.IS,
                                assembly = "ASM1692078v2",  # Specify your reference genome assembly
                                treatment = c(1, 1, 0),
                                pipeline = modkit_cols,
                                header = F,
                                dbtype = "tabix",
                                dbdir = "/Volumes/Xtreme/aim2/methylation/methylkit_out_IS/tabix",
                                context = "CpG",  # Context of methylation
                                mincov = 10)  # Minimum coverage to filter sites

# Merging and filtering
mk_meth.IS <- methylationData.IS |>
  filterByCoverage(lo.count = 10,
                   hi.perc = 99.9,
                   save.db = F) |>
  normalizeCoverage(save.db = F) |>
  unite(save.db = F)
```

##### Call differentially methylated regions (DMRs)

Count filtering, tiling into 100 bp windows, and differential
methylation analysis were performed with
[MethylKit](https://github.com/al2na/methylKit) according to [Flack et
al.,
2024](https://academic.oup.com/g3journal/article/14/8/jkae113/7683801).

``` r
# DMRs
mk_meth.all <- methylationData.all |>
  unite(min.per.group = 1L, save.db = FALSE)
tile_meth.all <- tileMethylCounts(mk_meth.all,
                              win.size = 100,
                              step.size = 100,
                              mc.cores = 8,
                              cov.bases = 10,
                              save.db = F)
tile_diff.all <- calculateDiffMeth(tile_meth.all,
                               slim = F,
                               save.db = F)
tile_delta.all <- getMethylDiff(tile_diff.all,
                            difference = 50,
                            qvalue = 0.05)
# Summary statistics
dmr_summary.all <- getData(tile_delta.all) |>
  summarize(n_DMR = n(),
            mean_abs_delta = mean(abs(meth.diff)),
            min_delta = min(meth.diff),
            max_delta = max(meth.diff)) |>
  mutate(across(everything(), ~ round(.x, digits = 1))) |>
  bind_cols(data.frame(windows_tested = nrow(getData(tile_diff.all)),
                       window_size = 100)) |>
  t() |>
  data.frame() %>%
  rename_with(function(x){x = "value"}) |>
  rownames_to_column(var = "stat")
# Write table with DMR information
write.table(tile_delta.IS_all[, c("chr", "start", "end")],
            file = "dmr_data_IS.bed",
            sep = "\t",
            quote = FALSE,
            col.names = FALSE,
            row.names = FALSE)
```

#### Visualization

Differential methylation was visualized by the fourteen longest
scaffolds in the
[PalLabHifi](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_016920785.2/)
Ixodes scapularis genome using
[karyotypeR](https://bernatgel.github.io/karyoploter_tutorial/). Once
again, this visualization is for one of the comparison groups, and the
analysis was repeated for three other comparison groups.

``` r
# Load packages
library(karyoploteR)
# Read in raw dmr data
dmr_all <- read_excel("/Volumes/Xtreme/aim2/nuc_genome/excel/5mc_dmr_analysis.xlsx", sheet="dmr_data_IS_all")
# Convert raw dmr data to GRanges
dmr_granges <- makeGRangesFromDataFrame(
  dmr_all,
  seqnames.field = "chr",  # Column with chromosome/scaffold names
  start.field = "start",   # Column with start positions
  end.field = "end",       # Column with end positions
  keep.extra.columns = TRUE  # Keep other columns if needed
)
# Read in significant DMRs
dmr_sig <- read_excel("/Volumes/Xtreme/aim2/nuc_genome/excel/5mc_dmr_analysis.xlsx", sheet="dmr_data_IS_sig")
# Convert significant DMRs to GRanges
dmr_sig_granges <- makeGRangesFromDataFrame(
  dmr_sig,
  seqnames.field = "chr",  # Column with chromosome/scaffold names
  start.field = "start",   # Column with start positions
  end.field = "end",       # Column with end positions
  keep.extra.columns = TRUE  # Keep other columns if needed
)
# Read in gene information
iscap_genes <- read_excel("/Volumes/Xtreme/aim2/nuc_genome/excel/high_impact_variants.xlsx", sheet="PalLabHifi_gene_annotation_IDs")
# Convert gene information to GRanges
iscap_genes_granges <- makeGRangesFromDataFrame(
  iscap_genes,
  seqnames.field = "CHROM",  # Column with chromosome/scaffold names
  start.field = "BIN_START",   # Column with start positions
  end.field = "BIN_END",       # Column with end positions
  keep.extra.columns = TRUE  # Keep other columns if needed
)
# Create the plot
kp <- plotKaryotype(genome = iscap.genome, chromosomes = c("NW_024609835", "NW_024609846", "NW_024609857", "NW_024609868", "NW_024609879",
                                                           "NW_024609880", "NW_024609881", "NW_024609883", "NW_024609882", "NW_024609884",
                                                           "NW_024609836", "NW_024609839", "NW_024609837", "NW_024609838"), plot.type = 2)
kpAddBaseNumbers(kp,
                 tick.dist = 50000000,  # Distance between ticks (adjust based on your genome size)
                 cex = 0.6,             # Text size for labels
                 add.units = TRUE,      # Adds 'Mb' units to the labels
                 digits = 2,            # Number of digits after the decimal in labels
                 minor.ticks = 0,
                 r0=0, r0=1)  
kpPlotHorizon(kp, data = dmr_granges, 
              y = dmr_granges$`meth.diff.MN-PAvTX`,
              y0 = 0,
              col = c("#27AEF9", "#E41A1C"), # You can adjust the color scheme
              border = NA,
              data.panel = 1)
kpPlotDensity(kp, data=dmr_sig_granges, data.panel = 2, col="#FFB000")
kpPlotDensity(kp, data=iscap_genes_granges, data.panel = 2, col="#AA88FF")
```
