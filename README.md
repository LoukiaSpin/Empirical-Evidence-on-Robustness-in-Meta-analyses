***
<img src="man/Robustness_Index.png" align="right" width="500" height="200"  />

# How robust are findings of pairwise and network meta-analysis when missing participant outcome data occur?: an empirical study 

## Description of the repository

The repository offers the typical structure of separate folders for data, and R (code/scripts), respectively.
* The _data_ folder includes four RData and text files. The RData files __binary_NMA__, __binary_PMA__, __continuous_NMA__, and __continuous_PMA__ refer to the datasets with the analysed pairwise meta-analyses (PMA) and network meta-analyses (NMA) on selected binary and continuous outcomes. The remaining RData files refer to the mean and standard deviation of selected empirically-based prior distributions for the between-trial variance that align with the outcome and intervention-comparison type of each analysis. Finally, the text files contain the posterior summaries from the Bayesian analysis (via the [R2jags](https://github.com/suyusung/R2jags/issues/) package) under all pre-defined scenarios for the missingness mechanism;
* The _R_ folder includes two analysis scripts (__A.Run empirical analysis_PMA & NMA.R__ and __Β.Determine robustness_PMA & NMA.R__). The first script sources the RData files with the datasets and performs the Bayesian analyses (one-stage random-effects PMA and NMA for continuous and binary outcomes with incorporation of the pattern-mixture model). The second script calculates the robustness index and apply our proposed decision framework for robustness of the primary analysis results for each analysis. Furthermore, it produces 1) the heatmap of robustness index for all possible comparisons of the network and 2) the panel of density plots of the summary log OR under all missingness scenarios and the Kullback-Leibler divergence. The network of Liu et al. (the motivating example in our article) has been used for that purpose. Finally, the remaining R scripts refer to the necessary functions to apply the aforementioned analysis scripts.<br>

[JAGS](http://mcmc-jags.sourceforge.net/) must be installed to employ the [R2jags](https://github.com/suyusung/R2jags/issues/) package. After downloading/cloning the repo, the user can use the .Rproj file to source all code.

The next sections briefly illustrate the functions of this repository.

## Output 

Prerequisite R packages: [R2jags](https://CRAN.R-project.org/package=R2jags), [dplyr](https://CRAN.R-project.org/package=dplyr), [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html), and [reshape2](https://cran.r-project.org/web/packages/reshape2/index.html)

### Data preparation 

To prepare the data in the proper format to use [R2jags](https://CRAN.R-project.org/package=R2jags), we have developed the function `data.preparation()` which has the following syntax:

```r
data.preparation(data, measure)
```

#### Explaining the arguments

* data: The input is a data-frame of a one-trial-per-row format containing arm-level data of each trial. This format is widely used for BUGS models. The columns of `data` refer to the following elements for a continuous outcome:
__t__ which is an intervention identifier;
__y__, the observed mean value of the outcome;
__sd__, the observed standard deviation of the outcome;
__m__, the number of missing participant outcome data. If a trial does **not** report this information for any investigated arm, insert `NA` in the corresponding arm(s); and
__n__, the number of participants randomised on the assigned intervention.
For a binary outcome, the columns of the data-frame `data` refer to the following elements:
__t__, the intervention identifier;
__r__, the number of observed events;
__m__, the number of missing participant outcome data. If a trial does **not** report this information for any investigated arm, insert `NA` in the corresponding arm(s); and
__n__, the number of participants randomised on the assigned intervention.
All elements appear in `data` as many times as the maximum number of interventions compared in a trial.
* measure: A character string indicating the effect measure with values `OR`, `MD`, `SMD`, or `ROM` for the odds ratio, mean difference, standardised mean difference and ratio of means, respectively.

#### Output of the function

This function returns a list with the necessary data to run [R2jags](https://CRAN.R-project.org/package=R2jags) through the function `sensitivity.analysis.mod()`.

### Sensitivity analysis 

We have developed the function `sensitivity.analysis.mod()` to run automatically Bayesian random-effects pairwise and network meta-analysis (with incorporation of the pattern-mixture model) for each pre-defined missingness scenario. The function has the following syntax:

```r
sensitivity.analysis.mod(data, measure, model, assumption, heter.prior, mean.misspar, var.misspar, D, n.chains, n.iter, n.burnin, n.thin)
```

#### Explaining the arguments

* data: See function `data.preparation()`.
* measure: See function `data.preparation()`.
* assumption: A character string indicating the analysis model with values `RE`, or `FE` for the random-effects and fixed-effect model, respectively.
* assumption: A character string indicating the structure of the informative missingness parameter. Set `assumption` equal to one of the following:  `IDE-ARM`, `IDE-TRIAL`, `IDE-COMMON`, `HIE-ARM`, `HIE-TRIAL`, `HIE-COMMON`, `IND-CORR`, or `IND-UNCORR`.
* heter.prior: A list of three elements with the following order: 1) a character string indicating the distribution with (currently available) values `halfnormal`, `uniform`, `lognormal`, or `logt`; and 2) two numbers for each parameter of the selected distribution. For `halfnormal`, `lognormal`, and `logt` these numbers refer to the mean and precision, respectively. For `uniform`, these numbers refer to the minimum and maximum value of the distribution. The `heter.prior` argument informs the `heterogeneity.param.prior()` function that determines the distribution to be used in the `prepare.model()` function (see section _Important details_).
* mean.misspar: A real number for the mean of the normal distribution of the selected informative missingness parameter (see argument `assumption`). This argument informs the `missingness.param.prior()` function to determine the 'mode' of `mean.misspar` as vector or scalar according to the argument `assumption`.
* var.misspar: A positive non-zero number for the variance of the normal distribution of the selected informative missingness parameter (see argument `assumption`).
* D: A binary number for the direction of the outcome. Set `D = 1` for a positive outcome and `D = 0` for a negative outcome. 
* n.chains: An integer specifying the number of chains for the MCMC sampling; an argument of the `jags()` function in [R2jags](https://CRAN.R-project.org/package=R2jags).
* n.iter: An integer specifying the number of Markov chains for the MCMC sampling; an argument of the `jags()` function in [R2jags](https://CRAN.R-project.org/package=R2jags).
* n.burnin: An integer specifying the number of iterations to discard at the beginning of the MCMC sampling; an argument of the `jags()` function in [R2jags](https://CRAN.R-project.org/package=R2jags).
* n.thin: An integer specifying the thinning rate for the MCMC sampling; an argument of the `jags()` function in [R2jags](https://CRAN.R-project.org/package=R2jags).

#### Output of the function

This function returns a list with two elements that contain an R2jags output: __EM__ which is the estimated effect measure of all possible comparisons of interventions under each scenario, and __tau__ which is the between-trial standard deviation (assumed to be common for all observed comparisons) under each scenario.

#### Important details
This function calls the `prepare.model()` function that contains all models in the BUGS language. The function `sensitivity.analysis.mod()` is found in the script __A.Run empirical analysis_PMA & NMA.R__.

### Robustness Index 

We have developed the function `robustness.index()` to calculate the __robustness index__ for the summary treatment effect of each comparison. The function has the following syntax:

```r
robustness.index(ES.mat, threshold, primary.scenar, nt)
```

#### Explaining the arguments

* ES.mat: The input is the element __EM__ from the function `sensitivity.analysis.mod()`. 
* threshold: The threshold of robustness. We suggest using 0.17 and 0.28 for a continuous and a binary outcome, respectively.
* primary.scenar: A number to indicate the primary analysis (here, the missing at random assumption). We have considered `primary.scenar = 13`, namely, the 13-th scenario.
* nt: The number of investigated interventions. In the case of pairwise meta-analysis, we have `nt = 2`. In the case of network meta-analysis, `nt` equals the number of interventions in the investigated network.

#### Output of the function

The function `robustness.index()` returns a list of three items:

1. a vector with the robustness index calculated for all pairwise comparisons, 
2. a list of the Kullback-Leibler divergence measure for all alternative scenarios for every pairwise comparison, and
3. a vector with the character strings _robust_ or _frail_ that indicate whether the primary analysis results for the corresponding comparison is robust or frail to the alternative missingness scenarios.

#### Important notes

The robustness index is specific to the pairwise comparison. Therefore, we can calculate only _one robustness index_ for a pairwise meta-analysis, but _as many as the number of possible comparisons_ in the network meta-analysis. For instance, in a network of four interventions, we have six possible comparisons, and hence, we can calculate a total of six robustness indeces. 

The function `robustness.index()` is found in the script __Β.Determine robustness_PMA & NMA.R__.

### Heatmap of Robustness Index 

To create the heatmap of robustness index, we have developed the function `HeatMap.AllComparisons.RI()` which has the following syntax:

```r
HeatMap.AllComparisons.RI(RI, drug.names, threshold)
```

#### Explaining the arguments 

* RI: A vector with the robustness index calculated for all pairwise comparisons. You need to use, first, the `robustness.index()` function. 
* drug.names: A vector with the names of the interventions compared. The interventions should be in the same order as in the argument `data` of the function `sensitivity.analysis.mod()`. This is important particularly in the case of network meta-analysis so that you do not match the robustness index with the wrong comparisons.
* threshold: The threshold of robustness. We suggest using 0.17 and 0.28 for a continuous and a binary outcome, respectively.

#### Output of the function

The function `HeatMap.AllComparisons.RI()` returns a lower triangular heatmap matrix which should be read from left to right. Each cell illustrates the robustness index for the corresponding pairwise comparison. A robustness index below the `threshold` implies present robustness for the corresponding comparisons (__green__ cells), whereas a robustness index equal or above the `threshold` implies lack of robustness (__red__ cells). 

The function `HeatMap.AllComparisons.RI()` is found in the script __Β.Determine robustness_PMA & NMA.R__.

### Density plots on the effect measure and Kullback-Leibler divergence under all missingness scenarios 

We have developed the function `KLD.plots()` to obtain a panel with density plots on the effect measure and Kullback-Leibler divergence under all missingness scenarios. The function has the following syntax:

```r
KLD.plots(ES.mat, primary.scenar, compar, outcome, drug.names)
```

#### Explaining the arguments 

* ES.mat: The input is the element __EM__ from the function `sensitivity.analysis.mod()`. 
* primary.scenar: A number to indicate the primary analysis (here, the missing at random assumption). We have considered `primary.scenar = 13`, namely, the 13-th scenario.
* compar: A number that indicates the comparison of interest. The function `possible.comparisons.id()` can be used to identify the number that corresponds to the comparison of interest. This function is relevant and useful for a network of interventions.
* outcome: A character string with values `binary` and `continuous` that indicates the outcome type.
* drug.names: A vector with the names of the interventions compared. The interventions should be in the same order as in the argument `data` of the function `sensitivity.analysis.mod()`. This is important particularly in the case of network meta-analysis so that you do not match the robustness index with the wrong comparisons.

#### Output of the function

The function `KLD.plots()` returns a five-by-five panel that illustrates three density plots simultaneously: one on the effect measure under the missing at random assumption (the primary analysis), one on the effect measure under an alternative informative missingness scenario, and one on the Kullback-Leibler divergence of these two missngness scenarios. 
  
The function `KLD.plots()` is found in the script __Β.Determine robustness_PMA & NMA.R__.
