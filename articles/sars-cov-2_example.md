# SARS-CoV-2 example

In this vignette, we show an example of how to modify the arguments to
simulate SARS-CoV-2 S gene evolution. Note that this is purely an
example and we are not experts in this area, so if you would like to
simulate SARS-CoV-2 evolution for a research project, we would recommend
looking into the specific parameter values yourself.

If you would like more detail about how to generate inputs or more
details about the
[`run_wavess()`](https://molevolepid.github.io/wavess/reference/run_wavess.md)
function, please see
[`vignette("prepare_input_data")`](https://molevolepid.github.io/wavess/articles/prepare_input_data.md)
and
[`vignette("run_wavess")`](https://molevolepid.github.io/wavess/articles/run_wavess.md),
respectively.

## Load libraries

First, we have to load the wavess library and a few additional libraries
that we’ll use for this vignette:

``` r

library(wavess)
library(ape)
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following object is masked from 'package:ape':
#> 
#>     where
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
library(tibble)
library(readr)
library(tidyr)
```

## Create Python virtual environment

Next, we must create a Python virtual environment. The virtual
environment is required because the underlying code of
[`run_wavess()`](https://molevolepid.github.io/wavess/reference/run_wavess.md)
is written in Python. You only have to create the virtual environment
once on each machine.

``` r

create_python_venv()
#> Warning in create_python_venv(): Skipped installation of the following
#> packages: numpyscipypandas Use `force` to force installation or update.
```

## Prepare input data

Next, we must generate the input data. If you’re interested in learning
more about how to prepare the input data, please see the corresponding
vignette (see
[`vignette("prepare_input_data")`](https://molevolepid.github.io/wavess/articles/prepare_input_data.md)).

We will only simulate 60 generations with sampling every day so the
simulation doesn’t take too long to run.

``` r

pop <- define_growth_curve(n_gens = 60 / 0.125)
samp <- define_sampling_scheme(sampling_frequency_active = 1, sampling_frequency_latent = 0) |> filter(day <= 60)
```

Let’s load in the founder sequence first:

``` r

fa_file <- system.file("extdata", "EPI_ISL_402124-S.fasta", package = "wavess")
founder <- extract_seqs(read.dna(fa_file, format = "fasta"),
  founder = "EPI_ISL_402124-S"
)
```

We will need to change the generation time, mutation rate, and
recombination rate.

Unlike HIV, SARS-CoV-2 has no latent reservoir, so we will turn this
feature off in the model.

generation time - 3 hours –\> 0.125 days  
<https://www.pnas.org/doi/10.1073/pnas.2406303121>

mutation rate - 1.5e-6 per site per generation  
<https://www.nature.com/articles/s41467-025-61555-x>

recombination rate - 2e-6 per site per year –\> 2e-6/(365/0.125) =
6.8e-10 per site per generation  
<https://www.nature.com/articles/s41467-022-31749-8>

``` r

gen_time <- 0.125
mut_rate <- 1.5e-6
rec_rate <- 6.8e-10
```

For simplicity, we will not model conserved sites or replicative
fitness, and we will use a Jukes Cantor model of nucleotide subtitution
rather than a SARS-CoV-2 specific Q matrix, although these would make
the simulation output more realistic.

``` r

q_mat <- matrix(
  c(
    -0.75, 0.25, 0.25, 0.25,
    0.25, -0.75, 0.25, 0.25,
    0.25, 0.25, -0.75, 0.25,
    0.25, 0.25, 0.25, -0.75
  ),
  ncol = 4
)
rownames(q_mat) <- c("A", "C", "G", "T")
colnames(q_mat) <- c("A", "C", "G", "T")
```

To get CTL epitopes for the founder amino acid sequence, we used the
IEDB T-cell prediction tool for MHC class I molecules with HLA alleles
A*02:01, A*01:01, B*07:02, B*08:01 and pMHC immunogenicity prediction,
both with the default settings. We filtered them in the same way as we
do for the HIV CTL epitopes. Website:
<https://nextgen-tools.iedb.org/pipeline?tool=tc1>

``` r

# get the amino acid sequence for use with IEDB
founder_aa <- ape::trans(ape::as.DNAbin(strsplit(founder$founder, split = "")[[1]]))
```

``` r

peptide_table_file <- system.file("extdata", "sars-cov-2-spike_peptide_table.tsv", package = "wavess")
ctl_epitopes <- read_tsv(peptide_table_file) %>%
  # arrange(`netmhcpan_el percentile`)
  filter(`netmhcpan_el percentile` < 0.1 & `immunogenicity score` > 1 / 90) %>%
  mutate(
    start = (start - 1) * 3,
    start = ifelse(`seq #` == 1, start, start + nchar(founder_ref_pol$founder)),
    days_to_full_potency = 1 / `immunogenicity score`
  ) %>%
  group_by(`seq #`, peptide, start) %>%
  summarize(days_to_full_potency = round(min(days_to_full_potency))) %>% # if multiple HLAs recognize the same epitope
  arrange(`seq #`, start) %>%
  rowwise() %>%
  mutate(
    anchor2 = strsplit(peptide, "")[[1]][2],
    anchor9 = strsplit(peptide, "")[[1]][9]
  ) %>%
  ungroup() %>%
  pivot_longer(starts_with("anchor"), names_to = "escape_position", values_to = "recognized_aa") %>%
  mutate(escape_position = as.numeric(gsub("anchor", "", escape_position))) %>%
  select(-c(peptide, `seq #`))
#> Rows: 5060 Columns: 13
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: "\t"
#> chr (4): peptide, allele, netmhcpan_el core, netmhcpan_el icore
#> dbl (9): seq #, start, end, peptide length, peptide index, median binding pe...
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
#> `summarise()` has regrouped the output.
```

We downloaded SARS-CoV-2 spike protein antibody epitopes from the LANL
database and used them in the same way as we use the HIV antibody
epitopes.  
Website (Spike feature download):
<https://cov.lanl.gov/components/sequence/COV/neutralization/download_db.comp>

``` r

antibody_file <- system.file("extdata", "sars-cov-2-spike-antibody-features.txt", package = "wavess")
antibody_epitopes <- read_tsv(antibody_file) %>%
  pull(Position) %>%
  get_epitope_frequencies() %>%
  sample_epitopes()
#> Rows: 467 Columns: 6
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: "\t"
#> chr (5): Antibody, Feature type, Study, Method, Annotation
#> dbl (1): Position
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
#> 26 resamples required
```

## SARS-CoV-2 simulation with `run_wavess()`

``` r

spike_sim <- run_wavess(
  inf_pop_size = pop,
  samp_scheme = samp,
  generation_time = gen_time,
  mut_rate = mut_rate,
  q = q_mat,
  recomb_rate = rec_rate,
  founder_seqs = rep(founder$founder, 10),
  b_epitope_locations = antibody_epitopes,
  t_epitope_locations = ctl_epitopes,
  act_to_lat = 0 # no latency
)
```

Summary:

``` r

spike_sim$counts %>%
  select(-c(contains("latent"), contains("replicative"), contains("conserved"))) # not modeling these so they are uninteresting
#> # A tibble: 61 × 8
#>    generation active_cell_count number_mutations number_recombinations
#>         <int>             <int>            <int>                 <int>
#>  1          0                10                0                     0
#>  2          8                80                0                     0
#>  3         16               529                0                     0
#>  4         24              1605               15                     0
#>  5         32              1970               12                     0
#>  6         40              1999                9                     0
#>  7         48              2000                9                     0
#>  8         56              2000               10                     0
#>  9         64              2000                6                     0
#> 10         72              2000                9                     0
#> # ℹ 51 more rows
#> # ℹ 4 more variables: mean_fitness_active <dbl>, mean_b_immune_active <dbl>,
#> #   mean_t_immune_active <dbl>, mean_immune_active <dbl>
```

Fitness:

``` r

spike_sim$fitness
#> # A tibble: 1,220 × 8
#>    generation seq_id   b_immune t_immune conserved replicative overall immune
#>    <chr>      <chr>       <dbl>    <dbl>     <dbl>       <dbl>   <dbl>  <dbl>
#>  1 founder    founder0        1        1         1           1       1      1
#>  2 founder    founder1        1        1         1           1       1      1
#>  3 founder    founder2        1        1         1           1       1      1
#>  4 founder    founder3        1        1         1           1       1      1
#>  5 founder    founder4        1        1         1           1       1      1
#>  6 founder    founder5        1        1         1           1       1      1
#>  7 founder    founder6        1        1         1           1       1      1
#>  8 founder    founder7        1        1         1           1       1      1
#>  9 founder    founder8        1        1         1           1       1      1
#> 10 founder    founder9        1        1         1           1       1      1
#> # ℹ 1,210 more rows
```

Sequences from active cells:

``` r

spike_sim$seqs_active
#> 1220 DNA sequences in binary format stored in a matrix.
#> 
#> All sequences of same length: 3822 
#> 
#> Labels:
#> founder0
#> founder1
#> founder2
#> founder3
#> founder4
#> founder5
#> ...
#> 
#> Base composition:
#>     a     c     g     t 
#> 0.295 0.189 0.184 0.333 
#> (Total: 4.66 Mb)
```

If you’d like to learn more about how to analyze the output, please
check out the post-processing vignette
[`vignette("analyze_output")`](https://molevolepid.github.io/wavess/articles/analyze_output.md).
