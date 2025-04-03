<div align="center">
<img src="./README_figures/IAMReX.png" alt="tittle" width="300">
<p align="center">
  <a href="https://arxiv.org/abs/2408.14140">
  <img src="https://img.shields.io/badge/arXiv-2408.14140-blue" alt="arxiv">
  </a>
</p>


[Overview](#Overview) -
[Features](#Features) -
[Examples](#Examples) -
[Install](#Install) -
[Tools](#Tools) -
[Acknowledgements](#Acknowledgements) -
[Contact](#Contact)

</div>


## Overview

This IAMReX repo extends the capability of original [IAMR](https://github.com/AMReX-Fluids/IAMR) codes, aiming at simulating the multiphase incompressible flows and fluid structure interaction problems on both CPUs and GPUs with/without subcycling. The Navier-Stokes euqations are solved on an adaptive semi-staggered grid using the projection method. The gas-liquid interface is captured using the level set (LS) method. The fluid-solid interface is resolved using the multidirect forcing immersed boundary method (IBM). The particle-wall as well as the particle-particle collisions are also captured by the adaptive collision time model (ACTM).



## Features

- LS method and reinitialization schemes
- Multidirect forcing Immersed Boundary Method
- Particle Collision Algorithms with Discrete Element Method



## Examples

- [Reversed Single Vortex (RSV)](./Tutorials/RSV/)

<div align="center">
    <img src="./README_figures/RSV.jpeg" alt="Profiles of drop interface in the RSV problem" width="500">
    <br>
    <figcaption style="text-align:center;">Profiles of drop interface in the RSV problem at t/T=1 after one rotation. Black line: Analytical Solution; Red line: 64*64; Blue line: 128*128; Green line: 256*256</figcaption>
    <br>
    <br>
</div>

- [Rayleigh-Taylor (RT) instability](./Tutorials/RayleighTaylor_LS/)

<div align="center">
    <!-- First Image -->
    <div style="display:inline-block; margin-right:10px; vertical-align:top; width:600px; height:400px; overflow:hidden;">
        <img src="./README_figures/IAMR_LSAMR.png" alt="Short Description 1" width="600">
        <br>
        <figcaption style="text-align:center;">(a) Density profile at t/T=2.42 using LS method. (b) Density profile at t/T=2.42 using IAMR convective scheme.</figcaption>
        <br>       
    </div>
    <!-- Second Image -->
	<br>
    <div style="display:inline-block; margin-left:10px; vertical-align:top; width:600px; height:400px; overflow:hidden;">
        <img src="./README_figures/RT_tip.png" alt="Short Description 2" width="600">
        <br>
        <figcaption style="text-align:center;">Comparison of the tip locations of the falling fluid and the rising fluid.</figcaption>
        <br>
        <br>
    </div>
</div>


- [Cluster of monodisperse particles](./Tutorials/Monodisperse/)

<div align="center">
    <img src="./README_figures/Monodisperse.png" alt="Cluster of monodisperse particles" width="600">
    <br>
    <figcaption style="text-align:center;"> Contours of velocity magnitude in y − z plane </figcaption>
    <br>
    <br>
</div>



## Install

### Download

Our codes rely on the following packages: 

1. AMReX:`git clone https://github.com/ruohai0925/amrex`
2. AMReX-Hydro:`git clone https://github.com/ruohai0925/AMReX-Hydro`
3. This repo:`git clone https://github.com/ruohai0925/IAMReX/tree/development`

After cloning, one will have three folder in your current directory:`AMReX`,`AMReX-Hydro`,`IAMReX`.

### Compile

We recommend using the GNU compiler to compile the program on the Linux platform. The compilation process requires preparing a `make` file, which you can find in the example folder under Tutorials. It is strongly recommended to use the `GNUmakefile` prepared in advance by the example. If you want to know more about the configuration information of the compilation parameters, you can check it in the [AMReX building](https://amrex-codes.github.io/amrex/docs_html/BuildingAMReX.html).

For example, if we want to compile in the `FlowPastSphere` , refer to the following steps:

1. `cd` to the FlowPastSphere directory

    ```shell
    cd IAMReX/Tutorials/FlowPastSphere
    ```

2. Modify compilation parameters in `GNUmakefile`.

   The compilation parameters depend on your computing platform. If you use `MPI` to run your program, then set `USE_MPI = TRUE`. If you are running the program with `Nvidia GPU(CUDA)`, then set `USE_CUDA = TRUE`. When using GPU runtime, please make sure your CUDA environment is ok.
   
3. Compile

   With the above settings, you can compile the program:

   ```shell
   make
   ```

   You can add parameters after `make`, such as `-f your_make_file` to customize your parameter file, and `-j 8` to specify `8` threads to participate in the compilation.

   If the compilation is successful, an executable file will be generated. Usually the executable file is named in the following format:`amr[DIM].[Compiler manufacturers].[computing platform].[is Debug].ex`. The executable file name of the Release version of the three-dimensional compiled using the GNU compiler in the MPI environment is: `amr3d.GNU.MPI.ex`.
   

### Usage

You should takes an input file as its first command-line argument.  The file may contain a set of parameter definitions that will overrides defaults set in the code. In the `FlowPastSphere`, you can find a file named `inputs.3d.flow_past_sphere`, run:

```shell
./amr3d.GNU.MPI.ex inputs.3d.flow_past_sphere
```

If you use `MPI` to run your program, you can type:

```shell
mpirun -np how_many_threads amr3d.GNU.MPI.ex inputs.3d.flow_past_sphere
```

This code typically generates subfolders in the current folder that are named `plt00000`, `plt00010`, etc, and `chk00000`, `chk00010`, etc. These are called plotfiles and checkpoint files. The plotfiles are used for visualization of derived fields; the checkpoint files are used for restarting the code. 



## Tools

In the `Tools` folder, there are some post-processing script for DIBM. The `postProcess` submodule is used to process the generated particle data file. The specific usage can be seen in the `README`. 

In addition, a Fortran code for generating a particle bed is provided, You can customize the domain size and particle size you want. 

```fortran
! Line 10 - 15
! domain scheme
LX_DOMAIN=3.0D0*Pi
LY_DOMAIN=6.0D0*Pi
LZ_DOMAIN=3.0D0*Pi
! D
SPD=1.D0
```

Compile the code and run it, and a file named `SPC.DAT` will be generated, in which each line represents the position of a particle.

```shell
gfortran -o fixBed DomainFill.F90
./fixBed
```

## Citation

To cite IAMReX, please use

```
@article{10.1063/5.0236509,
    author = {Li, Xuzhu (李虚竹) and Li, Chun (李春) and Li, Xiaokai (李晓凯) and Li, Wenzhuo (李文卓) and Tang, Mingze (唐铭泽) and Zeng, Yadong (曾亚东) and Zhu, Zhengping (朱正平)},
    title = {An open-source, adaptive solver for particle-resolved simulations with both subcycling and non-subcycling methods},
    journal = {Physics of Fluids},
    volume = {36},
    number = {11},
    pages = {113335},
    year = {2024},
    month = {11},
    issn = {1070-6631},
    doi = {10.1063/5.0236509},
    url = {https://doi.org/10.1063/5.0236509},
    eprint = {https://pubs.aip.org/aip/pof/article-pdf/doi/10.1063/5.0236509/20247663/113335\_1\_5.0236509.pdf},
}

@inbook{doi:10.2514/6.2025-1865,
    author = {Dewen Liu and Shuai He and Haoran Cheng and Yadong Zeng},
    title = {Investigate the Efficiency of Incompressible Flow Simulations on CPUs and GPUs With BSAMR},
    booktitle = {AIAA SCITECH 2025 Forum},
    chapter = {},
    pages = {},
    doi = {10.2514/6.2025-1865},
    URL = {https://arc.aiaa.org/doi/abs/10.2514/6.2025-1865},
    eprint = {https://arc.aiaa.org/doi/pdf/10.2514/6.2025-1865},
}

```

## Acknowledgements

We are grateful to Ann Almgren, Andy Nonaka, Andrew Myers, Axel Huebl, and Weiqun Zhang in the Lawrence Berkeley National Laboratory (LBNL) for their discussions related to [AMReX](https://github.com/AMReX-Codes/amrex) and [IAMR](https://github.com/AMReX-Fluids/IAMR). Y.Z. and Z.Z. also thank Prof. Lian Shen, Prof. Ruifeng Hu, and Prof. Xiaojing Zheng during their Ph.D. studies.



## Contact

If you have any question or wanna contribute to the code, please don't hesitate to contact us.
