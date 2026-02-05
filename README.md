## Basic Bulk Monte Carlo (MATLAB) — Electron Transport Simulator

### Overview

This project is a learning-focused **MATLAB Bulk Monte Carlo (BMC) simulator** for modeling **electron transport in bulk semiconductors** under an applied electric field. The simulator evolves a large ensemble of electrons in time using the standard Monte Carlo cycle:

1. **Initialize** an electron ensemble and physical constants
2. **Compute scattering rates** and build a sampling table
3. Repeatedly perform **free-flight drift** + **stochastic scattering**
4. **Collect and average** transport metrics (e.g., velocity/energy trends)
5. **Plot** scattering rates and simulation results

The code is organized into small, readable modules to make the core physics and algorithm easy to follow.

---

## What this simulator models

* **Carrier drift** under an external electric field (free-flight motion in k-space / real space depending on implementation choices)
* **Random scattering events** drawn from precomputed total and mechanism-specific scattering rates
* **Ensemble averaging** to extract macroscopic transport behavior from microscopic particle trajectories

This is the standard approach used to estimate quantities like steady-state drift velocity trends vs. field and to visualize how scattering shapes transport.

---

## Repository structure (high level)

Key MATLAB files (as organized in the repo):

* `initializePhysicalConstants.m` — Defines material/physical constants and simulation parameters
* `initializeElectronEnsemble.m` — Creates the initial electron population (energies/valleys/momenta as applicable)
* `calculateScatteringRates.m` — Computes scattering rates as a function of energy (or other state variable)
* `buildScatteringTable.m` — Builds a sampling structure (CDF/table) to efficiently choose scattering events
* `drift.m` — Implements the free-flight update under the applied field
* `runMonteCarloSimulation.m` — Main Monte Carlo loop (drift → scatter → record)
* `averageResults.m` — Post-processing / averaging across time windows or particles
* `plotScatteringRates.m`, `plotSimulationResults.m` — Visualization helpers
* `runner.m` — “Main” script that runs a typical simulation configuration
* `smoke_test.m` — Quick sanity test to validate basic execution

---

## Method summary (how the Monte Carlo loop works)

A typical bulk Monte Carlo transport iteration follows:

1. **Precompute scattering data**

   * Calculate mechanism-specific rates (e.g., acoustic/optical/intervalley depending on what’s implemented)
   * Form the **total scattering rate** and a **mechanism selection table**

2. **Free-flight (drift)**

   * Advance electron state under the applied field for a sampled free-flight time
   * Update momentum/energy accordingly

3. **Scattering**

   * Randomly select if/when scattering occurs and which mechanism occurs (via the table)
   * Apply the post-scatter state update (energy/valley/momentum changes)

4. **Statistics**

   * Accumulate ensemble averages (velocity, energy, etc.)
   * Summarize results using `averageResults.m` and plot outputs

---

## How to run

Typical workflow:

1. Open the project folder in MATLAB
2. Run:

   * `runner.m` for a standard simulation run
3. Optional:

   * `smoke_test.m` to validate basic functionality quickly

Outputs are visualized through the provided plotting scripts.

---

## Why this project matters (what it demonstrates)

This repo shows practical skills that directly map to semiconductor modeling work:

* Turning **transport physics** into a working stochastic simulator
* Efficient organization of a Monte Carlo codebase (rates → tables → loop → post-process)
* Clean separation between **physics models** (rates) and **numerical engine** (drift + event sampling)
* Reproducible, testable simulation workflow (`runner`, `smoke_test`)

