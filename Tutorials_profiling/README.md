In addition to the Tutorial, we provide four shell scripts to help you easily modify parameters and evaluate the performance of AMR (Adaptive Mesh Refinement) technology on different problems and platforms.

Take **jobiamr2dcpu** as an example:

This script supports several test cases, including: *Bubble, ConvectedVortex, DoubleShearLayer, FlowPastCylinder, LidDrivenCavity, RayleighTaylor, TaylorGreen.*

To run a specific case, simply uncomment the corresponding line within the script and you can execute multiple test cases at once.

For each case, the script allows you to configure different  combinations of the following parameters:

```
skip_level_projector = [0,1]
cycling = [Auto, None]
max_level = [1, 2, ...]
max_grid_size = [8, 16, ...]
regrid_int = [4, 8, ...]
```
You can configure a series of parameter combinations, execute them at once and collect the results efficiently.

The results are organized into a structured directory format as follows:
```
Bubble
├── case_results_cpu2d
│   ├── cpu2d_skip0_Auto_mgs16_1_regrid4
│   │   ├── inputs.2d.lid_driven_cavity
│   │   └── log.txt
│   ├── cpu2d_skip0_Auto_mgs16_1_regrid8
│   │   ├── inputs.2d.lid_driven_cavity
│   │   └── log.txt
│   ├── cpu2d_skip0_Auto_mgs16_2_regrid4
│   │   ...
├── case_results_gpu2
│   │   ...
├── case_results_cpu3d
│   │   ...
├── case_results_gpu3d
│   │   ...

```


* bubble is the case name,
* cpu/gpu is the running platform,
* 2d/3d is the case dimension.
* The inputs file contains the configuration and parameter combinations of the simulation, which is of course indicated by the name of the directory (cpu2d_skip0_Auto_mgs16_1_regrid4).
* The log contains the runtime and function calls, memory and other information, which is the main object of our analysis.

The other scripts correspond to different situation, such as 2D/3D and CPU/GPU, which can be easily identified by their names.

we present a small study we conducted, [
IAMReX_Tutorial_Profiling_Results
](https://github.com/hswind4/IAMR_Tutorial_Profiling_Results), which can help you better understand the tutorial profiling and its results. 



