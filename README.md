# TWPA Simulation and Analysis Suite

## 1. Project Overview

This project is a comprehensive MATLAB suite for the design, simulation, and analysis of a Travelling Wave Parametric Amplifier (TWPA). The core of the project is to model a superconducting microstrip transmission line whose dispersion is engineered by periodically modulating the length of capacitive "fingers" along its length.

The suite is designed to be a complete workflow, enabling a user to:
1.  **Simulate the fundamental device physics**: Calculate the frequency-dependent transmission (S21) and, most importantly, the dispersion relation (`k` vs. `f` or `β` vs. `f`).
2.  **Optimize the device geometry**: Systematically sweep device parameters (like modulation amplitude) to tune the dispersion for desired phase-matching characteristics.
3.  **Analyze phase-matching conditions**: Use the simulated dispersion to calculate the phase mismatch (`Δβ`) for four-wave mixing, which is the basis for parametric gain.
4.  **Simulate amplifier performance**: Model the gain, bandwidth, and generation of unwanted harmonics under a strong pump tone by solving the coupled-mode equations.

The ultimate goal is to provide a powerful in-silico tool to design a high-performance TWPA before fabrication.

---

## 2. Directory Structure

To keep the project organized, the files should be structured as follows. This separates the core simulators from analysis scripts, functions, and results.

```
TWPAv2/
├── 📜 README.md                 (You are here)
│
├── 🚀 Main_Simulators/
│   ├── network_calc.m            (Basic dispersion simulator with linear frequency grid)
│   └── network_calc_adaptive_freq.m (RECOMMENDED: Advanced simulator with adaptive grid)
│
├── ⚙️ Functions/
│   ├── mstrip.m                  (Calculates properties of a standard microstrip)
│   ├── mstrip_sc_Ls.m            (Calculates properties of a SUPERCONDUCTING microstrip)
│   ├── mb2.m                     (Calculates Mattis-Bardeen conductivity)
│   ├── abcd2s.m                  (Converts ABCD matrices to S-parameters)
│   ├── refine_freq_grid.m        (Implements the adaptive frequency meshing logic)
│   └── coupled_mode_odes.m       (Inferred: Defines the ODEs for harmonic simulation)
│
├──  Sweeps_And_Optimization/
│   └── run_parameter_sweep_with_summary.m (Sweeps a parameter and plots comparative results)
│
├── 📈 Analysis/
│   ├── analyze_delta_beta.m      (Calculates and plots phase mismatch Δβ)
│   └── simulate_harmonics_and_gain.m (Represents files like 'simulate0714harmonics_frequency_optismation.m')
│
├── 🗃️ Original/
│   ├── network_calc.m            (Legacy versions for historical reference)
│   └── network_calc_working.m
│
└── 💾 Results/
    ├── 0714TWPaX.mat             (Example output data file)
    ├── S21.png                   (Example output plot)
    └── dispersion.png            (Example output plot)
```

---

## 3. Detailed File Descriptions & Interactions

This section details the purpose of each script and how it connects to others.

### 3.1. Main Dispersion Simulators (`Main_Simulators/`)

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

### 3.2. Parameter Sweeps (`Sweeps_And_Optimization/`)

These scripts automate the process of running simulations for different device parameters to find an optimal design.

*   **`run_parameter_sweep_with_summary.m`**
    *   **Purpose**: To investigate how a specific geometric parameter affects the device dispersion.
    *   **How it Works**: This script contains a `for` loop that iterates over a range of values for a chosen parameter (e.g., `g.finger.modamp = [4, 5, 6, 7]`). In each iteration, it runs the core dispersion simulation logic (likely copied from `network_calc_adaptive_freq.m`). It stores the results of each run.
    *   **Usage**: Define the parameter and the range you want to sweep at the top of the script. Run it.
    *   **Output**: Produces a "summary" plot showing all the calculated dispersion curves on a single figure, color-coded or labeled by the parameter value. This is extremely useful for visualizing how to "engineer" the dispersion curve. It may also save a `.mat` file containing all the results.

### 3.3. Analysis and Performance Simulation (`Analysis/`)

These scripts take the dispersion data from the main simulators and use it to predict the actual performance of the TWPA as an amplifier.

*   **`analyze_delta_beta.m`**
    *   **Purpose**: To calculate the phase mismatch `Δβ` for four-wave mixing, which determines the potential for parametric gain.
    *   **Interaction**: It **loads** a `.mat` file generated by a dispersion simulator.
    *   **How it Works**:
        1.  Loads the `f` and `kperm` arrays from the results file.
        2.  Asks the user to define a pump frequency (`fp`).
        3.  It then iterates over a range of signal frequencies (`fs`). For each `fs`, it calculates the idler frequency `fi = fp - fs`.
        4.  Using interpolation (`interp1`), it finds the corresponding wavenumbers `kp`, `ks`, and `ki` from the loaded dispersion data.
        5.  It calculates the phase mismatch: `Δβ = kp - ks - ki`.
    *   **Output**: A plot of `Δβ` vs. signal frequency. The goal of dispersion engineering is to make `Δβ` as close to zero as possible over the widest possible bandwidth.

*   **`simulate_harmonics_and_gain.m`** (Represents files like `simulate0714harmonics_frequency_optismation.m`)
    *   **Purpose**: To simulate the full, nonlinear behavior of the amplifier, including signal gain, bandwidth, and the generation of unwanted higher-order harmonics.
    *   **Interaction**: This is the most complex script. It loads the dispersion data (`.mat` file) and likely takes pump power and device length as key inputs.
    *   **How it Works**: This script solves a system of coupled ordinary differential equations (ODEs).
        1.  It defines the set of interacting frequencies (pump, signal, idler, and various harmonics like `2*fp - fs`).
        2.  It uses the loaded dispersion data to find the `k` value for each of these frequencies.
        3.  It defines a function (e.g., `coupled_mode_odes.m`) that describes the derivatives `dA_i/dz` for the complex amplitude `A_i` of each wave. These equations include terms for the phase mismatch `Δβ` and a nonlinear coupling coefficient `γ` (which depends on pump power and device nonlinearity).
        4.  It uses a MATLAB ODE solver (like `ode45`) to integrate these equations along the length of the device (`z=0` to `z=L`).
    *   **Output**: Plots showing signal amplitude (gain) vs. device length or frequency, and the power in various harmonic "spurs". The "optimization" part of the filename suggests it may loop through different pump frequencies to find the one that gives the best gain-bandwidth product.

### 3.4. Core Physics & Utility Functions (`Functions/`)

These are the building blocks called by the main scripts. You do not run them directly.

*   **`mstrip.m`**: Calculates impedance (`Z0`) and effective dielectric constant (`epeff`) for a *standard* (non-superconducting) microstrip using the Hammerstad/Jensen formulas.
*   **`mstrip_sc_Ls.m`**: The superconducting version. It accounts for kinetic inductance by taking the London penetration depth (`lambda`) as an input. It is the core of the physics model for the transmission line itself.
*   **`mb2.m`**: Implements the Mattis-Bardeen theory to find the complex conductivity of the superconductor, which is needed to calculate the kinetic inductance.
*   **`abcd2s.m`**: A standard microwave utility to convert a 2x2 ABCD (chain) matrix into a 2x2 S-parameter matrix.
*   **`refine_freq_grid.m`**: The engine for the adaptive simulation. It wraps the core ABCD calculation in a loop that intelligently adds frequency points where the S21 phase changes rapidly.
*   **`coupled_mode_odes.m` (Inferred)**: A function that would be passed to an ODE solver. It takes the current position `z` and a vector of wave amplitudes `A` and returns the derivatives `dA/dz` based on the coupled-mode theory.

### 3.5. Legacy Files (`Original/`)

*   **`Original/network_calc.m`**, **`Original/network_calc_working.m`**: These are previous versions of the simulation scripts. They should **not** be used for new work but are kept for archival purposes. They can be useful for understanding the history of the project or reproducing an old result.

---

## 4. The Complete Workflow: A Step-by-Step Guide

Here is how you would use this suite to design and evaluate a TWPA from start to finish.

**Step 1: Initial Dispersion Simulation**
*   **Goal**: Get a first look at the dispersion of a new device design.
*   **Action**: Open `Main_Simulators/network_calc_adaptive_freq.m`. Enter the initial geometric and material parameters for your design. Run the script.
*   **Result**: A `.mat` file with the baseline dispersion `k(f)` and plots of `S21` and `k` vs `f`.

**Step 2: Design Optimization via Parameter Sweep**
*   **Goal**: Tune the dispersion curve to have the desired shape for broadband phase matching.
*   **Action**: Open `Sweeps_And_Optimization/run_parameter_sweep_with_summary.m`. Choose a parameter to tune (e.g., `g.finger.modamp`). Set a range of values for it. Run the script.
*   **Result**: A plot comparing the dispersion curves for each parameter value. From this, you can choose the optimal value. Update your design in `network_calc_adaptive_freq.m` with this new value and re-run it once to generate the final, optimized dispersion file.

**Step 3: Phase-Matching Analysis**
*   **Goal**: Quantitatively check how well your optimized design achieves phase matching.
*   **Action**: Open `Analysis/analyze_delta_beta.m`. Make sure it's set to load the `.mat` file from Step 2. Specify your intended pump frequency. Run the script.
*   **Result**: A plot of `Δβ` vs. signal frequency. This shows you the theoretical bandwidth over which you can expect gain.

**Step 4: Full Performance Simulation (Gain & Harmonics)**
*   **Goal**: Predict the real-world gain, bandwidth, and nonlinear performance of the final design.
*   **Action**: Open `Analysis/simulate_harmonics_and_gain.m`. Ensure it loads the correct dispersion file. Input the device length and expected pump power. Run the script.
*   **Result**: Plots of gain vs. frequency and the power of unwanted harmonic products. This is the final validation of your design.

---

## 5. Suggestions for Code Improvement & Future Work

This codebase is powerful and functional. To enhance its robustness and ease of use for future developers, consider the following improvements.

### 5.1. Refactor Core Calculation into a Function
The main calculation loop inside `network_calc.m` and `refine_freq_grid.m` is long and repeated. This logic should be extracted into a single, well-documented function to reduce code duplication and make the main scripts cleaner.

**Example Function Signature:**
```matlab
% [S21, kperm, f_out] = calculate_dispersion(g, material_params, sim_params)
```
This would dramatically simplify the main scripts and prevent bugs when changes are made in one place but not the other.

### 5.2. Centralize Parameter Management
Instead of hard-coding geometry and material parameters in every script, define them once in a configuration struct or a separate `config.m` file. This makes managing different designs and running sweeps much cleaner and less error-prone.

**Example Usage:**
```matlab
% In a main script:
device_config = get_twpa_design('Design_X_v1'); % A new function to load designs
[S21, kperm] = calculate_dispersion(device_config);
```

### 5.3. Create a Master Control Script
Write a single script `run_full_analysis.m` that orchestrates the entire workflow. It would call the dispersion calculation, then the delta-beta analysis, and finally the harmonic simulation, passing the results from one stage to the next automatically. This would make the entire suite much more user-friendly.

### 5.4. Improve Comments and Documentation
Add more comments explaining the *why* behind certain formulas, not just the *what*. For example, the source of the `vph = vph./sqrt(1+IcOverIstar^2+IcOverIstar^4)` correction should be cited (e.g., a paper or theory). This context is invaluable for new team members.

By implementing these software engineering best practices, the project will become even more powerful and maintainable for years to come.
