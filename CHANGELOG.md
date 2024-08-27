# SingleCell10x.jl changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed

* Fix bug when loading strange .h5 and .mtx files that are not sorted.

## [0.2.1] - 2023-09-13

### Fixed

* Add compat with HDF5.jl v0.17

## [0.2.0] - 2023-08-16

### Breaking

* `RawCSC` is replaced by `RawIJV` and no longer exported.
* It is no longer possible to pass a user-defined function to matrix loading functions.
* If type parameters are not specified, they now default to what is present in the source file.

### Fixed

* Significantly decreased time and memory usage when loading from .h5 files.
