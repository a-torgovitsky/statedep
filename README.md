## Nonparametric Inference on State Dependence in Unemployment

This repository contains code for reproducing the empirical estimates and Monte Carlo simulations in [Torgovitsky, A.](https://a-torgovitsky.github.io/), ["Nonparametric Inference on State Dependence in Unemployment," _Econometrica_, Vol. 87, No. 5 (September, 2019)](https://a-torgovitsky.github.io/statedep.pdf)

### Important

The code included in the supplemental material for _Econometrica_ is from May 3,
2019.<br/>
Please download the most recent version of the code from [the GitHub repository][GitHub].

### Software Requirements

* [MATLAB](https://www.mathworks.com/products/matlab.html). No special toolboxes are required.

* [A Mathematical Programming Language (AMPL)](http://ampl.com/). The trial/student version of AMPL is size restricted.
 Some limited versions of the code might run with the trial version, but a full license is required to reproduce the results in the paper.

* A linear programming solver for AMPL.  The default is [CPLEX](https://www.ibm.com/analytics/cplex-optimizer).
 The default can be changed by passing e.g.  `Settings.Solver = 'gurobi'`
 when calling `./src/DPO.m`.

* The [AMPL-MATLAB API](http://ampl.com/api/latest/matlab/getting-started.html)

* Linux (or perhaps OSX). I coded this on a Linux system and made no attempt to be platform-independent.
 However, the code is primarily in MATLAB, so should be mostly platform-independent.
 Some file operations are used for recording the results.
 These would be likely sources of issues for other operating systems, but should
 be easy enough to fix.

### Reproducing the Results

* **Important first step:**
 Open `./cfg/Config.m`.
    - Change `ResultsPath` to a directory where you want results to be saved.
    - Change `AMPLAPISetupPath` to the location where you installed the AMPL-MATLAB API.

* The primary code for the DPO model is contained is `./src/DPO.m` and the routines called from within.
 It contains many options, which are given default values in the structure called `Settings` that is defined at the top of that file.
 The code for the comparison parametric dynamic binary response (PDBR) model is
 contained in `./src/PDBR.m`.

* The directory `./bin/` contains a file called `RunSIPP.m` that generates the
  empirical results in the paper.

  - Selecting `SimSet = main` and `SimNum = n` will produce column `n` (for `n` between 1 and 12) in Table 2 of the paper.

  - Selecting `SimSet = sigma` or `SimSet = sigma-young` and `SimNum = n` will produce results for the `n`th gridpoint in Figures 2 and 3 (`n` between 1 and 9).

  - Selecting `SimSet = extra` will produce results for Appendix Table S4 in the
supplemental appendix.

* Running all of the empirical results in the paper will take a long
  time, primarily due to the procedure for constructing confidence regions.

  - Confidence regions can be turned off by changing `Settings.BuildConfidenceRegions = 1` to `Settings.BuildConfidenceRegions = 0` in the function `LoadSpec` in `./bin/RunSIPP.m`.

  - The number of bootstrap replications can be reduced by lowering `Settings.B
     = 250` to some smaller number, also in the function `LoadSpec` in
     `./bin/RunSIPP.m`.

  - Multiple results for each `SimSet` can be produced simultaneously by using
    the file `./bin/BatchRunSIPP.m`
    This is basically a poor-man's parallel that opens up multiple MATLAB
    threads. (Unfortunately, the AMPL-MATLAB API is not easy to parallelize,
    which is why I am using this crude workaround.)
    Using this function, all of the results of (e.g.) Table 1 can be
    produced with the command `BatchRunSIPP('your-save-dir', 'main')`.
    Note that this will open 12 MATLAB and AMPL instances at one time, which will strain a typical system.

* The directory `./bin/` also contains a file called `RunMonteCarlo.m` that
  generates the simulation results for the Monte Carlos reported in the
  supplemental Appendix.

  - `SimNumber = 1` produces the results for Table S1 and Figure S1.
  - `SimNumber = 2` produces the results for Table S2.
  - `SimNumber = 3` produces the results for Table S3.

* The directory `./bin/` contains a batching file for the Monte Carlos called
  `BatchRunMonteCarlo.m`. This opens three MATLAB threads that produce results
  for three different sample sizes for `SimNumber = 1` or `3`.

### Reproducing the Data

* The cleaned data used for both the empirical results and simulations is contained in
  `./data/sipp08.tsv` and `./data/sipp08-young.tsv`. These files are included
  with the repository. There is also a wide-form `./data/sipp08-wide.tsv` that I
  use for producing the table of summary statistics (Table 1).

* I have included a Bash script `./data/DownloadAndCleanSIPP.sh` that downloads the raw 2008 SIPP data from [the NBER page](http://www.nber.org/data/survey-of-income-and-program-participation-sipp-data.html), converts it to Stata format and then creates my extract:

  - The script uses the Stata dictionary and do files provided by the NBER, except
    that I have edited the beginning of the do files to remove their annoying
    directory hardcoding. These files are included with the repository in
    `./data`.

  - By default, the script deletes the ~8--10GB of SIPP data once my extract is created. This can be changed by running `./data/DownloadAndCleanSIPP.sh nocleanup` instead.

  - Part of the script involves running `./data/CleanSIPP.do` in Stata, which implements the sample selection rules discussed in the paper.

### Reproducing the Tables and Plots

The directory `./post` contains some Python scripts and LaTeX templates used to
turn the empirical and simulation results into tables and figures.

For the empirical results:
  - Table 1 (summary statistics) is generated by running
    `./post/BuildSumStatsTable.py ./data/sipp08-wide.tsv destdir` where `destdir` is
    the output location.

  - Table 2 (main empirical results) is generated by `./post/BuildResultsTable.py
    simdir/results/main` where `simdir` is the location of a simulation
    directory and `main` is the directory name that is created when the `SimSet`
    variable (see above) is set to `main` to produce the main results.
    In order to fully reproduce Table 2, all `SimNum`'s (1 through 12) should
    have been called so that there are directories `simdir/results/main/001`
    through `simdir/results/main/012`.
    The completed table will be located in `/simdirs/results/main/`.

  - Figure 2 (sensitivity analysis) is generated by
    `./post/BuildResultsTable.py simdir/results/sigma` where `simdir` is the
    location of a simulation directory and `sigma` is the directory name that
    is created when the `SimSet` variable (see above) is set to `sigma` to
    produce the sensitivity analysis.

  - Figure 3 (sensitivity analysis with young sample) is generated like Figure 2
    but with `sigma-young` in place of `sigma`.

  - Figure S4 (extra results) is generated like Table 2 but with `extra` in
    place of `main`.

For the Monte Carlo simulations:
  - Figure S1 and Table S1 are generated by running `./post/BuildMCEstTable.py
      simdir/results/` where `simdir` is the location of a simulation directory
      containing the results from `BatchRunMonteCarlo.m` with `SimNumber = 1`.

  - Similarly, figure S3 is generated by running `./post/BuildMCEstTable.py
      simdir/results/` after `BatchRunMonteCarlo.m` with `SimNumber = 3`.

  - Figure S2 is generated by running `./post/BuildMCTestTable.py
      simdir/results/` when after running `RunMonteCarlo.m` with `SimNumber =
      2`.

### My Software Versions

The results in the published paper were run with:

* MATLAB version 9.3.0.713579 (R2017b)
* AMPL version 20180927
* Gurobi version 8.1.0
* AMPL-MATLAB API version 1.4.0
* Stata 13.1 (for data cleaning only)

### Software Acknowledgments

The code uses two user-written MATLAB functions both located in `./ext`:

* [csvimport.m](https://www.mathworks.com/matlabcentral/fileexchange/23573-csvimport?focused=5196836&tab=function) written by Ashish Sadanandan

* [GaussQuad](https://www.mathworks.com/matlabcentral/fileexchange/26737-legendre-laguerre-and-hermite-gauss-quadrature?s_tid=prof_contriblnk) written by Geert Van Damme

### Problems or Bugs?

Please use [GitHub] to open an issue and I will be happy to look into it.

[GitHub]: http://www.github.com/a-torgovitsky/statedep
