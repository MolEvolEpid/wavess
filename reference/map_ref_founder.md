# Map reference and founder sequence positions

Map reference and founder sequence positions

## Usage

``` r
map_ref_founder(aln, ref, founder)
```

## Arguments

- aln:

  Alignment in [ape::DNAbin](https://rdrr.io/pkg/ape/man/DNAbin.html)
  format (e.g. read in with
  [`ape::read.dna()`](https://rdrr.io/pkg/ape/man/read.dna.html) or
  [`ape::read.FASTA()`](https://rdrr.io/pkg/ape/man/read.dna.html)) that
  includes the reference sequence and the founder sequence

- ref:

  Name of the reference sequence in the input alignment

- founder:

  Name of the founder sequence in the input alignment

## Value

Tibble including columns for:

- Alignment position (`alignment_pos`)

- Reference base position (`ref_pos`)

- Founder base position (`founder_pos`)

- Reference base (`ref_base`)

- Founder base (`founder_base`) All positions are indexed at 0.

## Examples

``` r
map_ref_founder(
  hxb2_cons_founder,
  "B.FR.83.HXB2_LAI_IIIB_BRU.K03455",
  "B.US.2011.DEMB11US006.KC473833"
)
#> # A tibble: 9,022 × 5
#>    alignment_pos ref_pos founder_pos ref_base founder_base
#>            <dbl>   <dbl>       <dbl> <chr>    <chr>       
#>  1           646     646           0 c        c           
#>  2           647     647           1 a        a           
#>  3           648     648           2 g        g           
#>  4           649     649           3 g        g           
#>  5           650     650           4 g        g           
#>  6           651     651           5 a        a           
#>  7           652     652           6 c        c           
#>  8           653     653           7 c        t           
#>  9           654     654           8 t        t           
#> 10           655     655           9 g        g           
#> # ℹ 9,012 more rows
```
