# wavess 1.2.0

* Add support for modeling both B-cell and T-cell immune responses in the
  within-host simulator, including separate inputs for B-cell and T-cell
  epitopes and corresponding fitness summaries in the output.
* Update the R `run_wavess()` wrapper, input checks, vignettes, and Python
  config files to reflect the new immune modeling options.

# wavess 1.1.0

* Allow for variable recombination rate across the sequence in the `run_wavess()` 
R function and the `run_wavess.py` Python script. This means you can now model
recombination hotspots, two genes without the intervening sequence, and reassortment.

# wavess 1.0.0

* Initial version.

