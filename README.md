# MCMC for Spike and Slab Survival Models

Giibs sampler for spike and slab survival models.

## Installation

```
git clone https://github.com/mkomod/mcmc_ss_surv
```

## Details

See [report](./sampler.pdf). 

Summary

 - Spike and slab consists of a Laplace slab and Dirac spike
 - Sampling requires a Metropolis-Hastings sampler for beta


## Example

Included is an example on simulated data, with 40% cenosred data. We ran our smaple for 10,000 iterations

![densities](./figures/grid.pdf)
