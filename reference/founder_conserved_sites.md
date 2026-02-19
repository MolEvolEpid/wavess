# Example founder conserved sites

Conserved sites (indexed at 1) for the gp120 DEMB11US006 sequence, which
is used as an example founder sequence throughout the package. These
sites were identified using the
[`identify_conserved_sites()`](https://molevolepid.github.io/wavess/reference/identify_conserved_sites.md)
function with the full HIV ENV filtered alignment from the LANL HIV
database. See `hiv_env_flt_2022` for more details and the link to the
entire alignment.

## Usage

``` r
founder_conserved_sites
```

## Format

### `founder_conserved_sites`

Vector of conserved nucleotides for the gp120 gene of DEMB11US006
(founder), named by the sequence position (indexed at 0)

## Source

Full alignment used:
<https://www.hiv.lanl.gov/content/sequence/NEWALIGN/align.html> Founder
sequence:
<https://www.sciencedirect.com/science/article/pii/S0022175914000143?via%3Dihub>
<https://www.hiv.lanl.gov/components/sequence/HIV/asearch/query_one.comp?se_id=KC473833>
