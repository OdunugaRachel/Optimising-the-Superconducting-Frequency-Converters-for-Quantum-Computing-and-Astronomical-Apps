# TWPA Simulation Project

## 1. Project Overview

This project contains a set of MATLAB scripts to simulate the performance of a Travelling Wave Parametric Amplifier (TWPA). The primary goal of these simulations is to calculate the device's transmission characteristics (S-parameters) and its dispersion relation (how the phase velocity of a signal changes with frequency).

The simulation models the TWPA as a long superconducting microstrip transmission line periodically loaded with capacitor-like "fingers". The length of these fingers is modulated along the device, which is the key mechanism for achieving dispersion engineering and, ultimately, parametric gain.

The main outputs of the simulation are:
1.  **S21 vs. Frequency Plot**: Shows the transmission magnitude through the device.
2.  **Dispersion Plot (`k` vs. Frequency)**: Shows the propagation constant (wavenumber `k` or `β`) as a function of frequency. This is the most critical output for designing the TWPA, as it determines the phase-matching conditions required for amplification.
3.  **Saved Data (`.mat` files)**: Key results like the frequency array, S21 data, and calculated dispersion are saved for further analysis.

---

## 2. File Structure & Purpose

The project files can be grouped into three main categories: Main Scripts, Function Files, and Legacy/Archival Files.

### 2.1. Main Simulation Scripts

These are the scripts you run to perform a simulation. You should only need to run one of these.

**`network_calc_adaptive_freq.m` (Recommended)**
*   **Purpose**: This is the most advanced and recommended script for running simulations. It calculates the S21 and dispersion of the TWPA.
*   **Key Feature**: It uses an *adaptive frequency grid*. It starts with a coarse frequency sweep and automatically adds more frequency points in regions where the phase of S21 is changing rapidly. This is crucial for accurately capturing sharp features in the dispersion curve (like near a bandgap) without needing an excessively large number of points across the whole spectrum, saving significant computation time.
*   **Usage**:
    1.  Open the file.
    2.  Modify the device parameters (geometry, materials, modulation) in the top sections.
    3.  Run the script.
    4.  The script will generate plots and save output `.mat` files.

**`network_calc.m`**
*   **Purpose**: A simpler version of the main simulation script.
*   **Key Difference**: This script uses a *fixed, linearly spaced frequency grid* (`linspace`). To accurately resolve sharp features, you may need to increase the number of points significantly (e.g., to 80,000+), which can make the simulation very slow.
*   **Usage**: It is functionally similar to the adaptive script but less efficient. It's useful for understanding the basic simulation flow without the complexity of the adaptive meshing.

### 2.2. Function Files (The "Engine")

These files contain the core physics and engineering models. They are called by the main scripts and you typically do not need to edit them.

**`mstrip.m`**
*   **Purpose**: Calculates the basic properties (characteristic impedance `Z0`, effective dielectric constant `epeff`) of a standard (non-superconducting) microstrip line.
*   **Method**: It implements the well-known Hammerstad and Jensen empirical formulas for microstrip lines.
*   **Called By**: This function is likely called by `mstrip_sc_Ls.m`.

**Other Important (but not provided) Functions:**
Based on the main scripts, the simulation relies on these other critical functions:
*   **`refine_freq_grid.m`**: This is the core of the adaptive frequency simulation. It likely takes an initial frequency grid and parameters, runs a coarse simulation, identifies large phase jumps, and inserts new frequency points in those regions, then re-runs the calculation until the phase is smoothly resolved. It is called by `network_calc_adaptive_freq.m`.
*   **`mstrip_sc_Ls.m`**: This function is central to the simulation. It calculates the per-unit-length inductance (`Lperm`) and capacitance (`Cperm`) for a *superconducting* microstrip. It accounts for the kinetic inductance of the superconductor, which is frequency-dependent and calculated from the London penetration depth (`lambda`). It likely calls `mstrip.m` for the geometric parts of the calculation.
*   **`mb2.m`**: Calculates the complex conductivity of the superconductor based on the Mattis-Bardeen theory. It takes the superconducting gap, temperature, and photon energy as inputs. This is essential for finding the kinetic inductance.
*   **`abcd2s.m`**: A standard microwave utility function that converts a 2x2 ABCD matrix (or "chain matrix") into a 2x2 S-parameter matrix.

### 2.3. Legacy/Archival Files (`Original/` directory)

The files in the `Original/` directory, such as `Original/network_calc.m` and `Original/network_calc_working.m`, appear to be earlier versions of the simulation scripts.

*   **Purpose**: They are kept for historical reference. You can see how the parameters and code have evolved.
*   **Usage**: You should not use these for new simulations. Refer to them only if you need to reproduce or understand a very specific old result. The scripts in the parent `TWPAv2/` directory are the current, active versions.

---

## 3. How the System Works: Simulation Workflow

The simulation follows a logical flow from defining the physical device to calculating its high-frequency behavior.

**File Call Hierarchy:**

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

## 4. Key Concepts for Analysis

### Delta Beta (`Δβ`) Analysis

The term "delta beta analysis" is not explicitly in the code, but the **output of this code is the input for that analysis**. Parametric amplification requires satisfying a phase-matching condition:

`Δβ = k_p - k_s - k_i ≈ 0`

where `k_p`, `k_s`, and `k_i` are the propagation constants (the `kperm` values from your simulation) for the pump, signal, and idler frequencies, respectively.

To perform this analysis, you would:
1.  Run the simulation to generate the `kperm` vs. `f` curve.
2.  Choose a pump frequency (`f_p`).
3.  For a range of signal frequencies (`f_s`), calculate the idler frequency (`f_i = f_p - f_s`).
4.  Look up the corresponding `k` values from your simulation results (`k_p = kperm(f_p)`, etc.).
5.  Calculate `Δβ` and plot it. The regions where `Δβ` is close to zero are where you can expect to get high parametric gain. The modulation (`g.finger.modamp`) is the primary tool used to shape the `kperm` curve to achieve this condition over a broad bandwidth.

### Parameter Sweeps

The current scripts are set up to run a single simulation. To perform a parameter sweep (e.g., to see how `g.finger.modamp` affects the dispersion), you would need to wrap the core logic of `network_calc_adaptive_freq.m` inside a `for` loop.

**Example of a sweep structure:**

```matlab
% Define a range of modulation amplitudes to test
mod_amps_to_sweep = [4, 5, 6, 7];

for i = 1:length(mod_amps_to_sweep)
    % Set the current modulation amplitude
    g.finger.modamp = mod_amps_to_sweep(i);

    % --- Paste the core simulation logic here ---
    % (From parameter setup down to plotting/saving)

    % It's a good idea to save the results with a unique name for each iteration
    % e.g., save(['results_modamp_' num2str(g.finger.modamp) '.mat'], ...);

    % You might also want to plot all dispersion curves on the same figure
    hold on;
    plot(f./1e9, kperm);
end
```

By doing this, you can systematically explore the design space and optimize the device geometry for the desired performance.

---

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