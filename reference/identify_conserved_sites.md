# Identify conserved sites

Identify conserved sites relative to the founder sequence. Note that the
positions returned assume that the start of the alignment is also the
start of the founder sequence in the simulation.

## Usage

``` r
identify_conserved_sites(
  aln,
  founder,
  thresh = 0.99,
  ref = NULL,
  founder_aln = NULL
)
```

## Arguments

- aln:

  Alignment in [ape::DNAbin](https://rdrr.io/pkg/ape/man/DNAbin.html)
  format (e.g. read in with
  [`ape::read.dna()`](https://rdrr.io/pkg/ape/man/read.dna.html) or
  [`ape::read.FASTA()`](https://rdrr.io/pkg/ape/man/read.dna.html)) that
  includes the founder sequence or a reference sequence, plus a
  representative set of sequences for your genome region of interest

- founder:

  Name of the founder sequence in the input alignment or in the optional
  founder-specific alignment (`founder_aln`)

- thresh:

  Conserved site threshold. A position is considered to be conserved if
  \>thresh proportion of sequences in the alignment are the same base
  (default: 0.99)

- ref:

  Name of the reference sequence in the input alignment, required if the
  founder sequence is not in `aln` (default: `NULL`)

- founder_aln:

  Alignment including the reference and founder sequences, required if
  the founder is not present in `aln`. **NOTE:** This alignment and
  `aln` are assumed to begin at the same position in the reference
  sequence (default: `NULL`)

## Value

Tibble including the following columns:

- `founder_pos`: founder position

- `founder_base`: founder base

- `consensus_base`: consensus base

- `consensus_prop`: proportion of sequences that had that base at that
  position

- `conserved`: whether or not the position is conserved (Yes means
  conserved, No means not conserved, NA means the conserved position is
  a gap ('-')) When using a reference, `NA` in the consensus columns
  indicates that that position was an insertion relative to the
  reference. All positions are indexed at 0.

## Examples

``` r
gp120_flt_2022 <- slice_aln(hiv_env_flt_2022, start = 1, end = 2517)
gp120_hxb2_cons_founder <- slice_aln(hxb2_cons_founder, start = 6225, end = 7757)
identify_conserved_sites(
  gp120_flt_2022,
  "B.FR.83.HXB2_LAI_IIIB_BRU.K03455"
)
#> # A tibble: 1,023 × 5
#>    founder_pos founder_base consensus_base consensus_prop conserved
#>          <dbl> <chr>        <chr>                   <dbl> <chr>    
#>  1           0 a            a                         1   Yes      
#>  2           1 t            t                         1   Yes      
#>  3           2 g            g                         1   Yes      
#>  4           3 a            a                         1   Yes      
#>  5           4 g            g                         1   Yes      
#>  6           5 a            a                         1   Yes      
#>  7           6 g            g                         1   Yes      
#>  8           7 t            t                         0.8 No       
#>  9           8 g            g                         1   Yes      
#> 10           9 a            a                         1   Yes      
#> # ℹ 1,013 more rows
identify_conserved_sites(gp120_flt_2022,
  "B.US.2011.DEMB11US006.KC473833",
  ref = "B.FR.83.HXB2_LAI_IIIB_BRU.K03455",
  founder_aln = gp120_hxb2_cons_founder
)
#> # A tibble: 1,473 × 5
#>    founder_pos founder_base consensus_base consensus_prop conserved
#>          <dbl> <chr>        <chr>                   <dbl> <chr>    
#>  1           0 a            a                         1   Yes      
#>  2           1 t            t                         1   Yes      
#>  3           2 g            g                         1   Yes      
#>  4           3 a            a                         1   Yes      
#>  5           4 g            g                         1   Yes      
#>  6           5 a            a                         1   Yes      
#>  7           6 g            g                         1   Yes      
#>  8           7 c            t                         0.8 No       
#>  9           8 g            g                         1   Yes      
#> 10           9 a            a                         1   Yes      
#> # ℹ 1,463 more rows
```
