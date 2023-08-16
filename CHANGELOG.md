# SingleCell10x.jl changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Breaking

* `RawCSC` is replaced by `RawIJV` and no longer exported.
* It is no longer possible to pass a user-defined function to matrix loading functions.
* If type parameters are not specified, they now default to what is present in the source file.

### Fixed

* Significantly decreased time and memory usage when loading from .h5 files.
