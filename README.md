# SimNet

<img src="man/figures/logo.png" width="150"/>

The **SimNet** R package provides a flexible framework to **SIM**ulate and infer disease transmission in populations connected through temporal contact **NET**works. It is designed to evaluate how different contact structures shape stochastic transmission processes and to identify the most likely contact networks and individual pathways underlying observed infection events. The package supports workflows such as:

-   Creating and managing simulated populations with individual infection-status trajectories
-   Simulating temporally correlated weighted contact networks
-   Simulating homogeneous-mixing contact networks as naïve transmission scenarios
-   Simulating directed object-mediated contact networks based on sequential use of shared resources
-   Simulating stochastic disease transmission through contact networks using SIS, SIR, and SIRS models
-   Estimating infection (beta) and recovery (gamma) parameters from temporal infection-status data
-   Comparing alternative candidate contact networks according to their probability of explaining observed infection events
-   Identifying transmission events that are incompatible with particular contact-network structures

🔗 **Project supporting the development of this package**:

-   [DigiGuard project](https://www.slu.se/en/research/research-catalogue/projekt/d/digiguard-project/)

<img src="man/figures/DGLogo.png" width="200"/>

📍 **Developed at:** Swedish University of Agricultural Sciences (SLU), Uppsala, Sweden

📅 **Version 1.0.0 release date:** 31 August 2026

------------------------------------------------------------------------

## Installation

You can install the development version of `SimNet` from GitHub using `remotes`:

``` r
# install.packages("remotes")
remotes::install_github("Hector-Marina/SimNet", build_vignettes=TRUE)
```

## Tutorial

A complete tutorial covering the main `SimNet` workflow is included as a package vignette. After installing `SimNet` with the vignettes enabled, it can be opened directly from R using:

``` r
browseVignettes("SimNet")
```

------------------------------------------------------------------------

## Authors

-   Hector Marina [![ORCID iD](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png)](https://orcid.org/0000-0001-9226-2902) **(Maintainer)**

-   Javier Sanchez [![ORCID iD](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png)](https://orcid.org/0000-0003-2605-8094)

-   Lars Rönnegård [![ORCID iD](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png)](https://orcid.org/0000-0002-1057-5401)

Any suggestions, bug reports, forks and pull requests are appreciated. Get in touch.

------------------------------------------------------------------------

## Citation

If you use `SimNet` in your research, please cite:

> *Marina, H., Sanchez, J., Woudstra, J., Doeschl-Wilson, A., Nielsen, PP., Nyman, A., & Rönnegård, L. (2026).* Utilizing data from real-time location systems to find the most likely network of transmission given observed infections. *Under revision*

------------------------------------------------------------------------

## 📖 Versioning

The `SimNet` package uses [semantic versioning](https://semver.org/).

------------------------------------------------------------------------

## 📜 License

The `SimNet` package is licensed under the [GPLv3](https://github.com/stewid/SimNet/blob/main/LICENSE).
