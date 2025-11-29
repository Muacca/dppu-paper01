# DPPU Paper 1: Microscopic Handles in Teleparallel Gravity

This repository contains the LaTeX source and supplementary SymPy code for the first DPPU paper ("microscopic handles in teleparallel gravity").

## Contents

- `paper/`: LaTeX source of the paper (main text and appendices).
- `code/`: SymPy scripts and logs used to check the Phase 3 stiffness and related cancellations.

## Build

A typical build with `latexmk`:

```bash
cd paper
latexmk -pdf main.tex
````

## Code and data availability

All SymPy scripts and logs used for the calculations in Appendix E are provided in the `code/` directory.

The canonical repository URL is:

`https://github.com/Muacca/dppu-paper01`

## License

* All **code** under `code/` is licensed under the MIT License (see `LICENSE-code`).
* All **text, figures, and LaTeX sources** under `paper/` are licensed under the Creative Commons Attribution 4.0 International License (CC BY 4.0; see `LICENSE-text`).


