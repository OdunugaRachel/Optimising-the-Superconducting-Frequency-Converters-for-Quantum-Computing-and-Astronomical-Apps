# TWPA Simulation Project

## 1. Project Overview

This project contains a set of MATLAB scripts to simulate the performance of a Travelling Wave Parametric Amplifier (TWPA). The primary goal of these simulations is to calculate the device's transmission characteristics (S-parameters) and its dispersion relation (how the phase velocity of a signal changes with frequency).

The simulation models the TWPA as a long superconducting microstrip transmission line periodically loaded with capacitor-like "fingers". The length of these fingers is modulated along the device, which is the key mechanism for achieving dispersion engineering and, ultimately, parametric gain.

The main outputs of the simulation are:
1.  **S21 vs. Frequency Plot**: Shows the transmission magnitude through the device.
2.  **Dispersion Plot (`k` vs. Frequency)**: Shows the propagation constant (wavenumber `k` or `β`) as a function of frequency. This is the most critical output for designing the TWPA, as it determines the phase-matching conditions required for amplification.
3.  **Saved Data (`.mat` files)**: Key results like the frequency array, S21 data, and calculated dispersion are saved for further analysis.

This document serves as a comprehensive guide to the project, explaining the simulation workflow, file structure, and how to use the scripts for analysis and device optimization.

---
## 2. Dispersion Simulaton Scripts

These are the foundational scripts. Their purpose is to calculate the dispersion relation `k(f)`, which is the input for all subsequent analysis.

*   **`network_calc_adaptive_freq.m` (Primary Script)**
    *   **Purpose**: Calculates the S21 transmission and dispersion (`k` vs. `f`) of the TWPA using an **adaptive frequency grid**.
    *   **How it Works**: It calls the `refine_freq_grid.m` function. This function performs an initial coarse simulation, identifies regions where the S21 phase changes rapidly (indicating interesting physics like a bandgap), and automatically inserts more frequency points in those regions. It repeats this until the phase is smoothly resolved. This provides high accuracy where needed without the massive computational cost of a dense linear grid.
    *   **Usage**: This is the **recommended script** for generating dispersion data. Modify the geometry (`g`), material, and modulation parameters at the top and run.
    *   **Output**: Saves a `.mat` file (e.g., `0714TWPaX.mat`) containing `f` (frequency array), `kperm` (dispersion data), and `S21`. Also saves plots.

*   **`network_calc.m`**
    *   **Purpose**: A simpler, non-adaptive version of the dispersion calculator.
    *   **How it Works**: Uses a `linspace` command to create a fixed, linearly spaced frequency grid. To resolve sharp features, you must manually increase the number of points to a very high value (e.g., 80,000+), making it much slower than the adaptive version.
    *   **Usage**: Good for educational purposes to understand the basic ABCD matrix calculation loop without the added complexity of adaptive meshing.
---

## 2.1. Simulation Workflow: From Physics to Code (Dispersion scipts)

The simulation follows a logical flow, translating the physical device structure into its high-frequency electromagnetic response.

**File Call Hierarchy:**
The main simulation script orchestrates calls to several core functions that model the underlying physics:
```
(User runs)
network_calc_adaptive_freq.m
 |
 +--> refine_freq_grid.m
       |
       +--> (Internal loop that calls the core calculation)
            |
            +--> mb2.m
            +--> mstrip_sc_Ls.m
            |    |
            |    +--> mstrip.m
            |
            +--> abcd2s.m
```

**Step-by-Step Simulation Process:**

1.  **Parameter Definition**: The user defines all physical parameters in the top section of a main script (e.g., `network_calc_adaptive_freq.m`). This includes:
    *   Finger geometry (`g.finger.l`, `g.finger.w`, `g.finger.p`).
    *   Microstrip geometry (`g.ms.h`, `g.ms.t1`, etc.).
    *   Modulation profile (`g.finger.modperiod`, `g.finger.modamp`).
    *   Material properties (`rhon`, `Tc`).

2.  **Frequency-Dependent Material Properties**: For each frequency in the sweep:
    *   The script calls `mb2.m` to calculate the Mattis-Bardeen conductivity.
    *   This is used to find the London penetration depth `lambda`, which quantifies the kinetic inductance of the superconductor.

3.  **Transmission Line Parameters**:
    *   The script calls `mstrip_sc_Ls.m` with the geometry and the calculated `lambda`. This function returns the per-unit-length inductance `Lperm(f)` and capacitance `Cperm(f)`.
    *   From these, the characteristic impedance `Z0(f)` and phase velocity `vph(f)` of the unloaded line are calculated.

4.  **Unit Cell ABCD Matrix Calculation**: The code models the entire TWPA as a long chain of identical "unit cells". A unit cell consists of the periodically loaded sections within one modulation period (`g.finger.modperiod`).
    *   The code enters a loop (`for jj = 1:n_unit_cell`) that builds the ABCD matrix for one complete modulation period.
    *   Inside this loop, it models each finger pair and the transmission line segment connecting them.
    *   A single finger-and-gap section is modeled as:
        *   A short piece of transmission line (`abcd_trl`).
        *   A shunt admittance representing the two capacitive fingers (`abcd_finger`).
        *   Another short piece of transmission line (`abcd_trl`).
    *   The `Lfinger` value is updated at each step `jj` using a `cos` function to model the sinusoidal modulation of finger lengths.
    *   These individual matrices are multiplied together to create the final ABCD matrix for one full modulation period (`abcd`).

5.  **Total Device ABCD Matrix**: The ABCD matrix of the single modulation period (`abcd`) is then multiplied by itself `ncell` times to get the ABCD matrix for the entire device (`abcd_tot`). This is a computationally efficient way to model a long, periodic structure.

6.  **S-Parameter Calculation**: The total device matrix `abcd_tot` is converted to S-parameters using the `abcd2s` function. The transmission parameter, `S21`, is stored.

7.  **Dispersion Calculation (Post-Processing)**:
    *   After the S21 is calculated for all frequencies, the phase of S21 is extracted (`angle(S21)`).
    *   The `unwrap` function is used to remove 2π jumps.
    *   A custom `while` loop is used to enforce that the phase is strictly monotonically increasing with frequency, which it must be for a physical passive device. This is a critical data-cleaning step.
    *   The final, clean phase `uth` is divided by the total device length to get the propagation constant `kperm = uth / len_meters`. This `kperm` is the `β` used in dispersion analysis.

8.  **Output and Analysis**: The script plots `abs(S21)` and `kperm` vs. frequency and saves the data.

---
## 3. Harmonic and Post-Processing Analysis

These scripts take the dispersion data generated by the main simulators as input to analyze the nonlinear performance of the TWPA.

*   **`simulate0714harmonics_frequency_optismation.m`**
    *   **Purpose**: Simulates nonlinear harmonic generation and **sweeps the pump frequency** to find the optimal operating point for frequency conversion.
    *   **How it Works**: It loads the dispersion data (`.mat` file) from `network_calc_adaptive_freq.m`. It then solves the coupled mode equations (CME) for a range of pump frequencies to model the power transfer between the pump and its odd harmonics (3rd, 5th, etc.).
    *   **Usage**: This is the primary script for analyzing and optimizing harmonic generation. Load a `.mat` file, set the pump parameters, and run.
    *   **Output**: Generates plots of harmonic output power versus pump frequency, helping to identify the best pump frequency for maximizing a specific harmonic.

*   **`simulate0714harmonics.m`**
    *   **Purpose**: Simulates nonlinear harmonic generation for a **single, fixed pump frequency**.
    *   **How it Works**: Similar to the frequency optimization script, it loads the dispersion data and solves the CMEs. However, it only does this for one specified pump frequency, showing how the power in each harmonic evolves along the length of the device.
    *   **Usage**: Useful for a detailed look at the power dynamics at a specific operating point you've already identified.
    *   **Output**: A plot showing the power of the pump and harmonics as a function of position along the TWPA.

*   **`simulateharmonics_withefficency_inclusion.m`**
    *   **Purpose**: Simulates harmonic generation while including efficiency calculations to determine the **optimal phase mismatch (Δβ)** for a given device.
    *   **How it Works**: For a working device configuration, this script solves the CMEs and calculates the conversion efficiency for each resulting harmonic. It then identifies the `Δβ` value that corresponds to the peak efficiency.
    *   **Usage**: Use this on a known-good device design to find the target `Δβ` for future optimization. This optimal `Δβ` is often slightly different from zero due to other nonlinear effects.
    *   **Output**: The optimal `Δβ` value that should be used as a target in `analyze_delta_beta.m` for new device designs.

*   **`analyze_delta_beta.m`**
    *   **Purpose**: Calculates and visualizes the **phase mismatch (Δβ)**, which is the critical parameter for achieving efficient parametric gain and harmonic generation.
    *   **How it Works**: It loads the dispersion data (`kperm` vs `f`). For a given nonlinear process (e.g., four-wave mixing: `Δβ = 2*k_p - k_s - k_i`), it calculates the phase mismatch `Δβ` across a range of frequencies.
    *   **Usage**: Run this after generating dispersion data to check if your device geometry supports the desired nonlinear process. While the ideal goal is to make `Δβ ≈ 0`, the true optimal value is often slightly non-zero. This optimal target `Δβ` can be found by running `simulateharmonics_withefficency_inclusion.m` on a working configuration.
    *   **Output**: A plot of `Δβ` vs. frequency. Regions where `Δβ` is close to the optimal target value are where the most efficient parametric processes can occur.

### 3.1. Harmonic and Post-Processing Analysis Workflow

The following workflow describes how the harmonic analysis scripts interact and how to use them to optimize nonlinear performance:

**File Call Hierarchy:**
The main harmonic analysis scripts use the dispersion data generated by the simulation scripts:
```
(User runs)
simulate0714harmonics_frequency_optismation.m
 |
 +--> solveCME.m
 |
 +--> analyze_5th_harmonic.m
 |
 +--> analyze_delta_beta_extra.m
 |
 +--> simulateharmonics_withefficency_inclusion.m
 |
 +--> analyze_delta_beta.m
```

**Step-by-Step Harmonic Analysis Process:**

1. **Load Dispersion Data:** Start by running `network_calc_adaptive_freq.m` (or `network_calc.m`) to generate and save the device's dispersion data (`.mat` file with `f`, `kperm`, `S21`).

2. **Sweep Pump Frequency for Harmonics:** Use `simulate0714harmonics_frequency_optismation.m` to sweep the pump frequency and solve the coupled mode equations for harmonic generation. This script helps identify the optimal pump frequency for maximizing a specific harmonic.

3. **Analyze Harmonic Power at Fixed Pump:** For a detailed look at harmonic evolution at a specific pump frequency, use `simulate0714harmonics.m`.

4. **Include Efficiency and Find Optimal Δβ:** Use `simulateharmonics_withefficency_inclusion.m` to calculate conversion efficiency and determine the optimal phase mismatch (Δβ) for your device.

5. **Phase Mismatch Analysis:** Run `analyze_delta_beta.m` (or `analyze_delta_beta_extra.m`) to visualize Δβ across frequencies and compare with the optimal target value. This guides device design for efficient parametric processes.

6. **Advanced Harmonic Analysis:** Use `analyze_5th_harmonic.m` for focused analysis of 5th harmonic generation, and `simulate_cascaded_triplers.m` for multi-stage/cascaded harmonic studies.

7. **Experimental Validation:** Use `harmonicPlotting.m` to compare simulation results with experimental data and validate device performance.

---

## 2.4. Extras & Bonuses

These scripts are central to understanding the nonlinear behavior and harmonic generation in the TWPA. They use the output of the network/diffraction simulation (e.g., `kperm`, `S21`) and model the coupled mode equations for parametric amplification and frequency conversion.

**`simulate_cascaded_triplers.m`**
*   **Purpose**: Simulates a two-stage device where the output of the first tripler (3rd harmonic generator) is used to pump a second tripler, enabling efficient generation of higher harmonics (e.g., 9th harmonic).
*   **Inputs**: Two `.mat` files for each stage, plus the output current from the first stage.
*   **Outputs**: Plots the evolution of harmonics in both stages, saves figures, and outputs the final harmonic currents.
*   **Usage**: Use this script to study cascaded harmonic generation and device optimization for high-order harmonics.

**`analyze_5th_harmonic.m`**
*   **Purpose**: Focused analysis of 5th harmonic generation, including frequency optimization and conversion efficiency vs. phase mismatch.
*   **Inputs**: Device `.mat` file, simulation parameters.
*   **Outputs**: Plots optimal pump frequency for 5th harmonic, conversion efficiency vs. phase mismatch, and saves results.
*   **Usage**: Use this script to optimize the device for strong 5th harmonic output.

**`harmonicPlotting.m`**
*   **Purpose**: Loads and visualizes experimental data (CSV files) from spectrum analyzer measurements, normalizes and overlays harmonic power traces for comparison with simulation.
*   **Inputs**: CSV files containing measured harmonic powers, reference S-parameter files.
*   **Outputs**: Plots normalized harmonic power vs. pump power, overlays multiple harmonics, and saves figures.
*   **Usage**: Use this script to analyze experimental results and validate simulation predictions.


# 4. Parameter Sweep Scripts

These scripts automate the process of sweeping device parameters (geometry, modulation amplitude, pump power, etc.) to optimize performance and explore the design space.

**`run_parameter_sweep_with_summary.m`, `run_parameter_sweep.m`**
*   **Purpose**: Systematically vary one or more device parameters, run the core simulation for each set, and summarize results (e.g., optimal gain, bandwidth, phase matching).
*   **Inputs**: Range of parameter values, device configuration.
*   **Outputs**: Plots and summary tables showing how performance metrics change with parameters.
*   **Usage**: Use these scripts to optimize device geometry, modulation, or operating conditions.

**`run_full_twpa_simulation.m`, `run_full_twpa_simulation_optimal_frequency.m`**
*   **Purpose**: Run the complete TWPA simulation for a given set of parameters, optionally finding the optimal pump frequency for gain or harmonic output.
*   **Inputs/Outputs**: Similar to above, but focused on full device simulation.

**`run_pump_power_sweep_with_summary.m`**
*   **Purpose**: Sweep pump power and analyze its effect on gain, harmonic generation, and efficiency.
*   **Usage**: Use to find optimal pump power for your device.

---

### Parameter Sweep Workflow & File Call Diagram

The following workflow describes how the parameter sweep scripts interact and how to use them to optimize device performance:

**File Call Hierarchy:**
The parameter sweep scripts orchestrate calls to the main simulation and analysis functions:
```
(User runs)
run_parameter_sweep_with_summary.m / run_parameter_sweep.m
 |
 +--> network_calc_adaptive_freq.m
       |
       +--> refine_freq_grid.m
             |
             +--> mb2.m
             +--> mstrip_sc_Ls.m
             |    |
             |    +--> mstrip.m
             |
             +--> abcd2s.m
       |
       +--> createTWPA.m
       +--> solveCME.m
       +--> removeHarmonics.m
       +--> ...other utility functions
```

**Step-by-Step Parameter Sweep Process:**

1. **Define Sweep Parameters:** Specify the range of device parameters to sweep (e.g., geometry, modulation amplitude, pump power) in the sweep script.

2. **Run Core Simulation:** For each parameter set, the sweep script calls the main simulation (`network_calc_adaptive_freq.m` or `network_calc.m`) to generate device data.

3. **Analyze Results:** The script collects outputs (e.g., gain, bandwidth, phase matching) and generates summary tables and plots for each parameter set.

4. **Optimize Device Design:** Use the summary results to identify optimal device configurations and operating conditions.

5. **Advanced Sweeps:** For full device simulations or optimal frequency searches, use `run_full_twpa_simulation.m`, `run_full_twpa_simulation_optimal_frequency.m`, or `run_pump_power_sweep_with_summary.m` as needed.

This workflow enables systematic exploration and optimization of the TWPA design space.

### 2.6. Core Functions Used Throughout

These functions are called by the main scripts and parameter sweeps. They encapsulate the physics and numerical methods:

*   **`createTWPA.m`**: Initializes a TWPA struct with all necessary fields for simulation (geometry, material, modes, currents, etc.).
*   **`solveCME.m`**: Solves the coupled mode equations for harmonic generation, given device parameters and pump frequency.
*   **`mb2.m`**: Calculates Mattis-Bardeen conductivity for superconductors.
*   **`mstrip_sc_Ls.m`**: Calculates inductance and capacitance for superconducting microstrip lines.
*   **`mstrip.m`**: Calculates standard microstrip properties.
*   **`abcd2s.m`**: Converts ABCD matrices to S-parameters.
*   **`removeHarmonics.m`**: Utility for filtering or post-processing harmonics.
*   **`refine_freq_grid.m`**: Implements adaptive frequency grid for efficient simulation.
*   **Other utility functions**: For plotting, data cleaning, and post-processing.
## 5. Suggestions for Code Improvement

The code is functional, but here are a few suggestions to improve its structure and reusability for future work.

1.  **Encapsulate Calculations into Functions**: The main calculation loop inside `network_calc.m` (and presumably inside `refine_freq_grid.m`) is very long. This core logic could be moved into its own function, for example:
    `[S21, Lperm, Cperm] = calculate_twpa_response(f, g, ms, ...)`
    This would make the main scripts cleaner and reduce code duplication.

2.  **Use a Configuration Struct/File**: Instead of hard-coding all the parameters at the top of the script, they could be defined in a separate configuration script or loaded from a `.mat` file. This makes it much easier to manage different device designs and run parameter sweeps without modifying the main calculation code.

3.  **Improve Variable Naming and Comments**:
    *   The variables `beta` and `betaf` are identical and can be consolidated.
    *   Add comments explaining the physical meaning of less obvious calculations, like the `xfactor` or the `vph` correction `vph = vph./sqrt(1+IcOverIstar^2+IcOverIstar^4)`. Explaining where these formulas come from (e.g., a reference to a paper or theory) would be invaluable.

4.  **Consolidate Constants**: Physical constants like `ev`, `mu0`, `kb` can be defined once in a separate script or function that is called at the beginning, rather than being redefined in every main script.

These changes would make the codebase more modular, easier to understand for new users, and more robust for future development.

---

## 2.7. Simulation Outputs & How to Use Them

**Simulation scripts output:**
- Plots of harmonic power vs. pump frequency, position, and pump power.
- Arrays of end currents for each harmonic (used for further analysis or as input to cascaded simulations).
- Conversion efficiency vs. phase mismatch (Δβ) plots.
- Saved `.mat` files with all relevant simulation data (frequency, S21, kperm, etc.).
- Figures and summary tables for parameter sweeps.

**How to use outputs:**
- Use harmonic power plots to identify optimal pump frequencies and device settings for desired harmonics.
- Use phase mismatch and efficiency plots to guide device design for strong harmonic generation.
- Use parameter sweep results to optimize geometry, modulation, and operating conditions.
- Compare simulation outputs to experimental data using `harmonicPlotting.m`.
- Use saved `.mat` files as input for further simulations or analysis scripts.

---

## 2.8. How Everything Fits Together

1. **Start with device/network simulation** (`network_calc_adaptive_freq.m` or `network_calc.m`) to generate `.mat` files with device parameters.
2. **Run harmonic simulation scripts** (`simulate0714harmonics_frequency_optismation.m`, `simulate_cascaded_triplers.m`, etc.) to model nonlinear behavior and harmonic generation.
3. **Analyze phase matching and efficiency** (`analyze_delta_beta_extra.m`, `analyze_5th_harmonic.m`) to optimize device design.
4. **Use parameter sweep scripts** to explore the design space and find optimal settings.
5. **Compare with experiment** using `harmonicPlotting.m` and overlay results for validation.
6. **Iterate and refine** device design and simulation parameters based on analysis and experimental feedback.

---

## 2.9. Summary Table of Key Files

| File Name                                      | Purpose/Functionality                                                      |
|------------------------------------------------|----------------------------------------------------------------------------|
| network_calc_adaptive_freq.m                   | Main device/network simulation (adaptive grid)                             |
| network_calc.m                                 | Main device/network simulation (fixed grid)                                |
| simulate0714harmonics_frequency_optismation.m  | Harmonic simulation, pump frequency sweep                                  |
| simulate_cascaded_triplers.m                   | Cascaded harmonic simulation (multi-stage)                                 |
| analyze_5th_harmonic.m                         | 5th harmonic optimization and efficiency analysis                          |
| analyze_delta_beta_extra.m                     | Phase mismatch (Δβ) analysis for harmonics                                 |
| harmonicPlotting.m                             | Experimental data analysis and visualization                               |
| run_parameter_sweep_with_summary.m             | Parameter sweep and summary analysis                                       |
| createTWPA.m                                   | TWPA struct initialization                                                 |
| solveCME.m                                     | Coupled mode equation solver                                               |
| mb2.m, mstrip_sc_Ls.m, mstrip.m, abcd2s.m      | Core physics/engineering functions                                         |
| removeHarmonics.m, refine_freq_grid.m          | Utility functions                                                          |
| *.mat files                                    | Saved device/network/simulation data                                       |

---

## 2.10. Best Practices for Future Users

- Always start with the device/network simulation to generate `.mat` files.
- Use harmonic simulation scripts to study nonlinear effects and optimize for desired harmonics.
- Use parameter sweep scripts to explore and optimize device design.
- Validate simulation results with experimental data using plotting scripts.
- Document any changes to device geometry, simulation parameters, or analysis methods.
- Save plots and results in organized output folders for reproducibility.
- Read comments in each script for detailed usage instructions.

---

# For questions or further development, start with the simulation and analysis scripts described above. Each script is documented and modular to support future research and device optimization.