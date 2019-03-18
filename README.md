
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MethylCal

This package allows the calibration of methylation levels using
**MethylCal**, a Bayesian calibration tool based on INLA ([Rue et
al., 2009](https://doi.org/10.1111/j.1467-9868.2008.00700.x)). The
package permits the visualisation of the calibration curve derived from
standard controls with known methylation percentages for a selected
assay. It also allows the correction of the methylation levels in case
and control samples. If cases are provided, MethylCal performs a
differential methylation test to detect hyper- and hypo-methylated
samples. Besides MethylCal method, the package also includes the
calibration and the correction provided by two alternative methods
([Warnecke et
al., 1997](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC147052/);
[Moskalev et al., 2011](https://www.ncbi.nlm.nih.gov/pubmed/21486748)).

## Installation

**MethylCal** has two dependencies,
[**lattice**](https://cran.r-project.org/web/packages/lattice/index.html)
and
[**latticeExtra**](https://cran.r-project.org/web/packages/latticeExtra/index.html),
that are automatically downloaded and installed. Another dependency
[**R-INLA**](http://www.r-inla.org/) needs to be downloaded and
installed manually. This is achieved by typing in the R command
    line

    install.packages("INLA", repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"), dep = TRUE)

Once **R-INLA** has been installed (it may take a couple of minutes),
installation of **MethylCal** requires the following steps:

1.  Install the [**devtools**](https://github.com/r-lib/devtools)
    package. This can be done from
    [**CRAN**](https://cran.r-project.org/). Invoke R and then type
    
    ``` 
     install.packages("devtools")
    ```

2.  Load the **devtools** package
    
    ``` 
     library(devtools)
    ```

3.  Install **MethylCal** package by
    typing
    
    ``` 
     devtools::install_github("lb664/MethylCal", build_opts = c("--no-resave-data", "--no-manual"))
    ```

4.  Finally, load the **MethylCal** package
    
    ``` 
     library(MethylCal)
    ```

## Example

The example below allows to replicate **Figure 3** C\&D in Ochoa et.
(2019).

    data(BWS_data)
    AMP = c(0, 25, 50, 75, 100)
    data = Formatting(BWS_data, AMP = AMP, n_Control = 15, n_Case = 18)
    corr_data = MethylCalCorrection(data, Target = "KCNQ1OT1", n_Control = 15, n_Case = 18, printing = FALSE)

The output obtained by running `MethylCalCorrection` is transformed into
a table with the same data format used in the uploaded data set with the
function

    output_data = TableFormatting(corr_data, n_Control = 15, n_Case = 18)

## Further information

Further information can be found in the accompanying vignette by typing

    vignette("MethylCal")
