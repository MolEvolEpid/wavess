# Slice alignment

Slice alignment

## Usage

``` r
slice_aln(aln, start, end, seqs = labels(aln))
```

## Arguments

- aln:

  alignment

- start:

  start position in alignment

- end:

  end position in alignment

- seqs:

  sequences to keep (default: labels(aln), i.e. all sequences)

## Value

sliced alignment in
[ape::DNAbin](https://rdrr.io/pkg/ape/man/DNAbin.html) format

## Examples

``` r
slice_aln(hxb2_cons_founder, 1, 100)
#> 3 DNA sequences in binary format stored in a matrix.
#> 
#> All sequences of same length: 100 
#> 
#> Labels:
#> B.FR.83.HXB2_LAI_IIIB_BRU.K03455
#> CON_B(1295)
#> B.US.2011.DEMB11US006.KC473833
#> 
#> Base composition:
#>    a    c    g    t 
#> 0.32 0.28 0.19 0.21 
#> (Total: 300 bases)
```
