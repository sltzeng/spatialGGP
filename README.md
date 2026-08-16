# Monte Carlo Comparison of Scalable Approximations for Generalized Gaussian Processes in Geostatistical Count Data

This repository contains the R code for the simulation study presented in the associated paper.

The code provides a benchmark comparison of several scalable methods for fitting generalized Gaussian process models to geostatistical count data. The simulation code is provided to facilitate the reproducibility of the Monte Carlo experiments presented in the paper and to provide a benchmark framework for comparing the statistical and computational performance of these approaches.

## Requirements
**Please make sure that your R version is:**

> R >= 4.4.3


**The simulation code relies on the following R packages:**

> spaMM,
> INLA,
> sp,
> sf,
> geoR,
> MASS,
> fields,
> glmmTMB,
> fmesher,
> inlabru,
> spStack,
> callr

All packages other than INLA can be installed through the usual R package installation process:

```
install.packages(c(
  "spaMM",
  "inlabru",
  "glmmTMB",
  "spStack",
  "sp",
  "sf",
  "fmesher",
  "geoR",
  "callr",
  "fields"
))
```

For installation instructions for INLA, please refer to the official INLA website: <br>
[https://www.r-inla.org/download/index.html](https://www.r-inla.org/download/index.html)

## Methods Compared

The simulation code compares the following approaches:

| Method | R package |
|:---:|:---:|
| spaMM | spaMM |
| INLA | INLA |
| glmmTMB | glmmTMB |
| spStack | spStack |
| ensGP | implemented in-house |

The methods are evaluated under the simulation settings considered in the paper, allowing their statistical and computational performance to be compared within a common Monte Carlo framework.

## Important Note on glmmTMB
The glmmTMB method has been tested and validated on Linux platforms, where it is known to run reliably. Therefore, **we strongly recommend using Linux when running the glmmTMB part of the simulation.**

On other operating systems, including Windows and macOS, the underlying TMB/C++ implementation may terminate unexpectedly. To prevent such failures from interrupting the entire simulation procedure, the code may return NA values when the corresponding subprocess fails.

Consequently, if messages such as "callr subprocess failed" appear repeatedly when running the glmmTMB method, the resulting values should be regarded as unreliable and interpreted with caution.

## Variable Names and Notation
Note that the variable names used in the simulation code are not identical to the notation used in the paper.

The correspondence is as follows:

| Paper notation | Code notation | Description |
|---|---|---|
| $\mathbf{s}$ | `s = (x, y)` | Spatial location |
| $s_1, s_2$ | `x, y` | Two spatial coordinates |
| $y$ | `z` | Observed count |
| $z$ | `latent field` | Latent spatial process |

In the paper:

- $\mathbf{s}$ denotes the spatial location;
- $y$ denotes the observed count;
- $z$ denotes the latent spatial field.

In the code:

- `s = (x, y)` denotes the spatial location, where x and y are the two spatial coordinates corresponding to $s_1$ and $s_2$ in the paper;
- `z` denotes the observed count, corresponding to $y$ in the paper.
- This difference in notation should be taken into account when relating the simulation code to the mathematical notation in the paper.

## Citation
> If you use this code in your research, please cite the associated paper.
> 
> ShengLi Tzeng and Meng-Ting Huang. (2026+) Monte Carlo comparison of scalable approximations for generalized Gaussian processes in geostatistical count data. Journal of the Chinese Statistical Association. Accepted.

## Disclaimer
This repository contains research code intended primarily for simulation, benchmarking, and reproducibility purposes. The code has not necessarily been optimized for general-purpose use outside the simulation settings considered in the paper.
