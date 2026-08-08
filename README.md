# Summary

This project contains the estimator and controller for the canards that flew on the 2025 rocket, Aurora, and further development towards the 2026 Waterloo Rocketry flight. 
Additionally, it contains a 6DOF rocket plant model designed to enable closed loop simulation, including sensor dynamics. 
There are a number of setup and support scripts for the main model (such as evaluating the Barrowman equations from input geometry). 

Ideally, setting the rocket to simulate should be as easy as changing the first line in `configure_plant_model.m` to run the appropriate script. 
It probably isn't, so at that point you should ask someone for help.

# How to use

## Running the Sim

0. Open `Simulinkcanards.prj` in MATLAB so that the repository folders are added to the path

#### Option A: `run_save_plot` (preferred)

1. Run `run_save_plot` from the MATLAB command window. This configures and runs the SIL simulation, saves the log to `results/sil/result.mat`, and displays both the live and saved-log plots.

Common variants:

```matlab
run_save_plot                                  % Run SIL for 180 seconds
run_save_plot("hil", 120)                      % Run HIL for 120 seconds
run_save_plot("sil", 120, "PlotMode", "log")   % Only plot the saved log
run_save_plot("sil", 120, "PlotMode", "none")  % Run and save without plotting
```

Arguments:

| Argument | Kind | Default | Description |
| --- | --- | --- | --- |
| `mode` | Positional | `"sil"` | Simulation mode: `"sil"` or `"hil"`. |
| `StopTime` | Positional | `180` | Simulation duration in seconds. |
| `ResultPath` | Name-value | Value of `mode` | Subdirectory beneath `results/` in which to save the output. |
| `FileName` | Name-value | `"result.mat"` | Name of the saved log file. |
| `PlotMode` | Name-value | `"both"` | Plot behavior: `"live"`, `"log"`, `"both"`, or `"none"`. |
| `ParameterOverrides` | Name-value | `struct()` | Scalar struct containing per-run overrides for configured model parameters. |

#### Option B: Run from Simulink

1. Run `configure_plant_model` in MATLAB.
2. Open `plant-model/CC_Flight_Simulation.slx`.
3. After the model has loaded, click **Run** in the Simulink toolbar.
4. View the scope blocks in `/visualization_estimator` or `/plant_combined/visualization_sim`.

## Setup
1. Clone the repo `git clone https://github.com/waterloo-rocketry/closedrocket.git`
2. Make sure you have MATLAB **2025b** installed (The specific version matters cause Simulink ;-;)
3. Install (do this with Matlab install when you can select multiple at once, if possible):
    - Aerospace Blockset
    - Aerospace Toolbox
    - Control System Toolbox
    - DSP System Toolbox
    - Instrument Control Toolbox
    - MATLAB Support for MinGW-w64 C/C++/Fortran Compiler
    - Signal Processing Toolbox
    - Simulink
4. In Matlab run `mex -setup C` and `mex -setup C++`
