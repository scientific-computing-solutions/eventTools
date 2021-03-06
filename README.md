# eventTools

[![Build Status](https://travis-ci.org/scientific-computing-solutions/eventTools.svg?branch=master)](https://travis-ci.org/scientific-computing-solutions/eventTools)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/eventTools)](https://cran.r-project.org/package=eventTools)
[![Coverage Status](https://coveralls.io/repos/scientific-computing-solutions/eventTools/badge.svg?branch=master&service=github)](https://coveralls.io/github/scientific-computing-solutions/eventTools?branch=master)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/scientific-computing-solutions/eventTools?branch=master&svg=true)](https://ci.appveyor.com/project/scientific-computing-solutions/eventTools)
[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/eventTools)](https://cran.r-project.org/package=eventTools)

Extension to the eventPrediction package for N-piecewise Weibull and Weibull mixture models.

The eventTools package contains some extensions of functionality for the eventPrediction
package. It will use the same interface when possible and the idea is to mimic the
command-sequences used and propagate the functionality (accrual processes, dropouts etc).
The first add-in is a method to allow for changes in the event rates at specific
time-points resulting from, for example, a lag-period.  The second add-in implements
a Weibull mixture model for dealing with subgroups. The data is aggregated during
an ongoing trial and will be assumed blinded. This is related to the Predict from
data functionality and is available in version $\ge 1.0$ of the eventTools package.
The plan is to include more functionality to this package with iterative releases.

## Contributors
Dalevi, Daniel (maintainer); Ruau, David;

## Installation

To install the development version from GitHub:
```R
install.packages("devtools")
# eventTools use eventPrediction
devtools::install_github("scientific-computing-solutions/eventPrediction", 
                         build_vignettes = TRUE)
# We spent a lot of time developing the vignettes. We recommend the read but 
# building them from source takes some time
devtools::install_github("scientific-computing-solutions/eventTools", 
                         build_vignettes = TRUE)
```
