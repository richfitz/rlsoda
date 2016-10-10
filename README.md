# rlsoda

[![Project Status: WIP - Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)
[![Travis-CI Build Status](https://travis-ci.org/richfitz/rlsoda.svg?branch=master)](https://travis-ci.org/richfitz/rlsoda)
[![Coverage Status](https://img.shields.io/codecov/c/github/richfitz/rlsoda/master.svg)](https://codecov.io/github/richfitz/rlsoda?branch=master)

Highly experimental.

## Installation

For now, run `./setup.sh` in the root to clone/update the upstream liblsoda sources.  Yes, this could be done with submodules, but they're a total faff.  Soon I'll include the sources directly as we'll need that for travis really, and to work with things like `install_github`.

After that, running `make install` will install the package into your default R library.

## The interface

The interface here is related to `dde` and differs slightly from `deSolve`.  Probably the `dde` / `rlsoda` differences can be smoothed over.  The interface is tuned so that all modifications to the integrated data are optional so that where speed is really important (e.g., running an MCMC chain with the output of the ODE) unnecessary modifications are not done.  This will require some modifications for functions expecting `deSolve`-like output.

## The library

This package uses [Simon Frost](https://github.com/sdwfrost)'s [`liblsoda`](https://github.com/sdwfrost/liblsoda).  For ease of installation, this is bundled into the package here.

`liblsoda` is a thread-safe C rewrite of the venerable Fortran `LSODA` algorith (from the Linda Petzold and Alan Hindmarsh).
