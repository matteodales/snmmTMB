# Code supplement: *A Semiparametric Nonlinear Mixed Effects Model with Penalized Splines Using Automatic Differentiation*

This folder contains all files needed to reproduce the figures and tables of the manuscript and its online appendices.
The method is referred to as **snmmTMB** ; the two competitors are
**assist** (`assist::snm()`, Ke and Wang 2001) and **snmmAGQ** (our
implementation of the adaptive Gauss–Hermite quadrature estimator of Elmi et
al. 2011, in `src/snmmAGQ_functions_*.R`).

---

## Layout

```
code/
├── README.md                     this file
├── run_all.R                     master script to run all simulations
├── install_packages.R            installs the CRAN packages used (any version)
├── renv.lock                     exact package versions used for the paper (renv::restore())
│
├── src/                          compiled code and shared R functions
│   ├── starting_points.cpp                          penalized-spline GAM for the starting values (Section 2.3.2)
│   ├── snmmTMB_likelihood_sinecurve_simulation.cpp  snmmTMB likelihood, sine-curve model (Section 3)
│   ├── snmmTMB_likelihood_bellcurve_simulation.cpp  snmmTMB likelihood, bell-curve model (Section 3)
│   ├── snmmTMB_likelihood_smocc_application.cpp     snmmTMB likelihood, SMOCC height model (Section 4)
│   ├── snmmTMB_derivpenalty.cpp                     penalized spline with monotonicity penalty (Appendix D)
│   ├── snmmAGQ_functions_sinecurve.R                snmmAGQ estimator, sine-curve model
│   └── snmmAGQ_functions_bellcurve.R                snmmAGQ estimator, bell-curve model
│
├── simulation_sinecurve/         Section 3, sine curve: Figures 1 and 2
│   ├── sinecurve_example.R                     one replication, all three methods, plotted
│   ├── sinecurve_server_compile.R              compiles the two TMB models
│   ├── sinecurve_server_simulation_<method>.R  the Monte Carlo study, one script per method
│   ├── sinecurve_server_run_<method>.sh        cluster launchers (Linux, `screen` + `taskset`)
│   ├── sinecurve_server_merge_and_plots.R      merges the job files, draws Figures 1 and 2
│   ├── results/sinecurve_server_<method>_results.RDS   merged results
│   ├── Figure1.pdf, Figure2.pdf
│   └── scalability/                Section 3, timing study: Figure 3
│       ├── sinecurve_scalability.R, sinecurve_server_run_scalability.sh
│       ├── sinecurve_scalability_merge_and_plots.R
│       ├── results/sinecurve_scalability_results.RDS
│       └── Figure3.pdf
│
├── simulation_bellcurve/         Section 3, bell curve: Figure 4  (same file pattern as above)
│   ├── bellcurve_example.R, bellcurve_server_compile.R
│   ├── bellcurve_server_simulation_<method>.R, bellcurve_server_run_<method>.sh
│   ├── bellcurve_server_merge_and_plots.R
│   ├── results/bellcurve_server_<method>_results.RDS
│   └── Figure4.pdf
│
├── application_smocc/            Section 4, SMOCC height data: Table 1, Figures 5 and 6
│   ├── data/smocc_200.csv                      the data (brokenstick::smocc_200)
│   ├── smocc_model_snmmTMB.R                   snmmTMB fit          -> results/smocc_snmmTMB_results.RDS
│   ├── smocc_model_sitar.R                     sitar fit + bootstrap -> results/smocc_sitar_results.RDS
│   ├── smocc_comparison_snmmTMB_vs_sitar.R     Table 1 (results/smocc_comparison_summary.txt) and Figure5.pdf
│   ├── smocc_parambootstrap.R                  parametric bootstrap -> results/smocc_bootstrap_results.RDS, Figure6.pdf
│   └── Figure5.pdf, Figure6.pdf
│
└── appendices/
    ├── appendixB/                non-Gaussian sine-curve variants: Supplementary Figures 1 and 2
    │   ├── t_residuals/          t-distributed residuals
    │   └── t_randomeffects/      t-distributed random effects
    ├── appendixC/                effective degrees of freedom diagnostic: Supplementary Figure 3
    │   ├── bellcurve_edf_diagnostic.R
    │   ├── results/bellcurve_edf_diagnostic_results.RDS
    │   └── Supplementary_Figure3.pdf
    └── appendixD/                monotonicity penalty illustration: Supplementary Figure 4
        ├── snmm_monotonesimulation.R
        └── Supplementary_Figure4.pdf
```



## Running the supplement

### 1. Install

R ≥ 4.1 and a C++ compiler (Rtools on Windows, Xcode command-line tools on
macOS, `g++` on Linux) are needed for `TMB::compile()`. Then either

```r
# quick: current CRAN versions
Rscript install_packages.R
```

or, for the exact versions used for the paper,

```r
install.packages("renv"); renv::restore(lockfile = "renv.lock")
```

### 2. `run_all.R`

Open `run_all.R` in RStudio, set the working directory to this folder
(Session > Set Working Directory > To Source File Location), choose what to
run:

|  | What it runs | Approximate time |
|---|---|---|
| `figures <- TRUE` | every `*_merge_and_plots.R` script plus the application and appendix plot scripts: all figures and Table 1 redrawn from the available results | ~5 min the first time (compiles two TMB models), ~2 min after |
| `examples <- TRUE` | `sinecurve_example.R` and `bellcurve_example.R`: one simulated dataset each, fitted with all three methods and plotted against the truth | ~1 min |
| `application <- TRUE` | `smocc_model_snmmTMB.R` and `smocc_model_sitar.R`: the SMOCC fits from `data/smocc_200.csv`; `figures` then rebuilds Table 1 and Figures 5–6 from them | ~7 min |
| `simulations <- TRUE` | the job scripts of the five simulation studies | ~2 400 CPU-hours: use a cluster, see below |

Each script runs in its own `Rscript` session from its own folder. The run
stops at the first script that fails.

Two figures have a recompute switch inside their own script rather than in
`run_all.R`: `run_simulation <- TRUE` in
`appendices/appendixC/bellcurve_edf_diagnostic.R` recomputes Supplementary
Figure 3 (~25 min), and `run_bootstrap <- TRUE` in
`application_smocc/smocc_parambootstrap.R` reruns the 500-replicate bootstrap
of Figure 6 (~7 h, resumes if interrupted).

### 3. On a cluster

The `*_server_run_*.sh` launchers start `NJOBS` independent `Rscript` processes
(`Rscript <script> <job> <njobs>`) in detached `screen` sessions, pinned to
one core each (Linux; `screen`, `taskset`, single-threaded BLAS). Job `j`
handles replications `(setting, iter)` with `((task - 1) mod njobs) + 1 == j`
and writes only `results/<prefix>_job<j>_of<njobs>.RDS`, flushed after every
replication, so an interrupted job resumes where it stopped. Run
`*_server_compile.R` once first, and the `*_merge_and_plots.R` script with
`do_merge <- TRUE` once all jobs are done. The results do not depend on
`njobs`.

---

## Simulation studies (Figures 1–4, Supplementary Figures 1–2)

| Study | Folder | Settings | Replications per setting | Merged result files | CPU time |
|---|---|---|---|---|---|
| Sine curve (Fig. 1–2) | `simulation_sinecurve/` | 16 = n ∈ {10, 50} × m ∈ {10, 20} × σ ∈ {0.4, 1} × D ∈ {diag(0.25, 0.16, 0.04), diag(1, 0.25, 0.16)} | 500 | `sinecurve_server_<method>_results.RDS` | ~610 h |
| Scalability (Fig. 3) | `simulation_sinecurve/scalability/` | 17 cells: n ∈ {10, 20, 40, 60, 80, 100} at m = 10; m ∈ {5, 10, 20, 40, 80} at n = 20; K ∈ {8, 10, 15, 20, 30, 45} at (20, 10), snmmTMB and snmmAGQ only | 50 | `sinecurve_scalability_results.RDS` | ~65 h (single-threaded) |
| Bell curve (Fig. 4) | `simulation_bellcurve/` | 8 = n ∈ {10, 50} × m ∈ {10, 20} × σ ∈ {0.2, 0.4}, D = [2 1; 1 2] | 500 | `bellcurve_server_<method>_results.RDS` | ~180 h |
| t(3) residuals (Suppl. Fig. 1) | `appendices/appendixB/t_residuals/` | the 16 sine-curve settings | 500 | `sinecurve_server_t_residuals_<method>_results.RDS` | ~820 h |
| t(3) random effects (Suppl. Fig. 2) | `appendices/appendixB/t_randomeffects/` | the 16 sine-curve settings | 500 | `sinecurve_server_t_randomeffects_<method>_results.RDS` | ~730 h |
---

CPU time is the sum of the per-replication fit times stored in the merged
results (`metrics$time`).

---

## Application (Table 1, Figures 5–6)

The data are the 200 children of `brokenstick::smocc_200` (Van Buuren,
`brokenstick` R package), exported to
`application_smocc/data/smocc_200.csv`:
1 942 height measurements on 200 subjects. The scripts drop rows with missing
values, leaving 1 904 observations (`smocc_comparison_summary.txt`).

* `smocc_model_snmmTMB.R` fits the model of Section 4 (about 2 min including
  the bands); `smocc_model_sitar.R` fits `sitar::sitar()` with the spline df
  chosen by BIC over 3–20 and bootstraps its curve (200 replicates, about
  5 min); `smocc_comparison_snmmTMB_vs_sitar.R` writes Table 1 to
  `results/smocc_comparison_summary.txt` (and prints it) and draws Figure 5
  from the two saved fits (about 15 s, plus 2 min compiling the TMB model the
  first time).
* `smocc_parambootstrap.R` simulates 500 datasets from the fitted snmmTMB
  model, refits each (about 50 s per replicate, 7 h in total) and draws
  Figure 6. With `run_bootstrap <- FALSE` (the shipped default) it only
  redraws the figure from `results/smocc_bootstrap_results.RDS`, in seconds.

---

## Appendices

* **Appendix B** (`appendices/appendixB/`): the sine-curve study with t-distributed residuals and random effects.
  The scripts are identical to the sine-curve scripts with the data-generating step changed.
* **Appendix C** (`appendices/appendixC/bellcurve_edf_diagnostic.R`): the three
  methods are fitted to 50 datasets of bell-curve setting 5 (n = 10, m = 10,
  σ = 0.4) and the effective degrees of freedom of their estimated shape
  function compared (Supplementary Figure 3). About 25 min on a laptop;
  `run_simulation <- FALSE` redraws from the manuscript results RDS.
* **Appendix D** (`appendices/appendixD/snmm_monotonesimulation.R`): a
  penalized spline fitted to y = x³ − x + ε under the monotonicity penalty
  weights λ_c ∈ {0, 0.1, 1, 10} (Supplementary Figure 4; n = m = 20, σ = 0.4, K = 10), using
  `src/snmmTMB_derivpenalty.cpp`. About 10 s, plus a
  one-off compilation of the model.

---


## References

* Cole, T. J., Donaldson, M. D. C. and Ben-Shlomo, Y. (2010). SITAR — a useful
  instrument for growth curve analysis. *International Journal of Epidemiology*.
* Elmi, A., Ratcliffe, S. J., Parry, S. and Guo, W. (2011). A B-spline based
  semiparametric nonlinear mixed effects model. *Journal of Computational and
  Graphical Statistics*.
* Ke, C. and Wang, Y. (2001). Semiparametric nonlinear mixed-effects models
  and their applications. *Journal of the American Statistical Association*.
* Kristensen, K., Nielsen, A., Berg, C. W., Skaug, H. and Bell, B. M. (2016).
  TMB: Automatic differentiation and Laplace approximation. *Journal of
  Statistical Software*.
