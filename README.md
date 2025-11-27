Wrapped Gaussian Process Modeling
=================================

Description
-----------

This repository contains the source code for an R package designed to fit wrapped Gaussian Processes (WGPs). WGPs are appropriate for when the response of interest is an angle, in which our model must fit to data on a circle. This is common when working in applications that deal with directional data, such as meteorology or physics.

The model implemented in this package is discussed in the following paper in review:

> Andrew Cooper, et al. "Robust Wrapped Gaussian Process Inference for Noisy Angular Data." 2025.

Installation
------------

wrapgp is written in R. You can install the package locally by cloning this repository, going to the folder where it was saved, and running in R: 

```r
install.packages("wrapgp")
```

License
-------

This source code is protected under copyright license O# (O4843). Complete license information can be found in the file LICENSE.md.

Maintaners
----------
- Andrew Cooper, <ahcooper@lanl.gov>
- Justin Strait, <jstrait@lanl.gov>
- Mary Frances Dorn, <mfdorn@lanl.gov>
- Brendon Parsons, <bparsons@lanl.gov>
- Alessandro Cattaneo, <cattaneo@lanl.gov>
