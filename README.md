# ineqReSample
This R package provides tools to correct citation inequality measures for marginals bias using resampling (Kim, Adolph, West, and Stovel, 2020).

**Background**:  Kim et al. (2020) show that measures of inequality are subject to "marginals bias," which impairs the ability to compare measures across time or contexts. Consider the example of citation distributions, which describe the relative shares of incoming citations sent to each member of a population of published papers, e.g., for a specific field and year. Inequality measures like the Gini coefficient (and even typically robust measures such as quantiles of the distribution) are only comparable across disciplinary fields or time periods if both the total number of papers potentially receiving citations and the total number of citations sent to this population remain constant. If these "marginals" are rising, inequality may be understated, and if they are falling, inequality may be overstated.

**Method**: This function corrects inequality measures for marginals bias by resampling the observed citations to have common marginals fixed to some reference level (the method suggested by Kim et al, 2020).  Users may request either that inequality measures be adjusted to be comparable to those of a base period  or to a fixed, user-provided count of papers and citations.  While these functions are designed to adjust citation inequality metrics specifically, they should have more general utility for any inequality measure that pertains to the concentration of “signals” sent by one population to another.

**Output**: Available corrected inequality measures include the percentage of ever cited papers, the Gini coefficient, the Herfindahl-Hirschman index, and quantile measures.

This package was created by Christopher Adolph and Lanu Kim (maintainer).

**Package citation**:  

Christopher Adolph and Lanu Kim.  2020.  "ineqReSample." R Package.  Version 0.1. https://github.com/lanukim/ineqReSample

**References**:

Lanu Kim, Christopher Adolph, Jevin West, and Katherine Stovel.  Forthcoming.  "The Influence of Changing Marginals on Measures of Inequality in Scholarly Citations: Evidence of Bias and a Resampling Correction." Sociological Science.

Please report any package related issues here: https://github.com/lanukim/ineqReSample/issues

