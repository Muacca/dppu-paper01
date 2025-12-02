# DPPU Paper 1: Microscopic Handles in Teleparallel Gravity

This is the English README.  
For the Japanese version, see:  
- ➡️ [日本語版 README](./README_ja.md)
- ➡📄 paper PDF: [main.pdf](./main.pdf)
- DOI: [10.5281/zenodo.17762204](https://doi.org/10.5281/zenodo.17762204)

---

## Overview

This repository hosts the LaTeX source and supplementary SymPy code for the first DPPU paper,  
“Microscopic Handles in Teleparallel Gravity”.

In this paper we work within teleparallel gravity (TEGR) and consider microscopic spacetime
handles of topology $S^2 \times S^1$ carrying an integer torsional monopole charge $q$.
We investigate how far geometry alone can account for:

- a stable handle radius and a mass scale with $E \propto |q|$,
- chirality selection through the Nieh–Yan term, and
- a precession-mode stiffness $k(q) \propto \omega^2 |q|$ with classical critical exponent $\gamma = 1$.

---

## Repository layout

- `paper/`  
  LaTeX sources of the paper (main text, appendices, references, figures).

- `code/`  
  SymPy scripts and logs used in Appendix E.  
  These cover, for example, the Phase 3 stiffness calculation and the explicit
  verification of the $q^2 \varepsilon^2$ cancellations.

---

## Build instructions

We assume `latexmk` is available in your LaTeX environment.  
A typical build looks like:

```bash
cd paper
latexmk -pdf main.tex
````

Standard LaTeX packages such as `amsmath`, `amssymb`, `graphicx`, and `hyperref`
are expected to be installed.

---

## Code and data availability

All SymPy scripts and raw text logs used for the computations in Appendix E are
included in the `code/` directory. They implement, for example:

* perturbative expansion of the torsion tensor,
* automatic extraction and cancellation checks of the $q^2 \varepsilon^2$ terms, and
* verification of the scaling $k(q) \propto \omega^2 |q|$.

The canonical repository URL is:

`https://github.com/Muacca/dppu-paper01`

---

## License

This repository uses different licenses for code and text:

* All **code** under `code/` is released under the MIT License
  (see `LICENSE-code`).
* All **text, figures, and LaTeX sources** under `paper/` are released under the
  Creative Commons Attribution 4.0 International License (CC BY 4.0;
  see `LICENSE-text`).




