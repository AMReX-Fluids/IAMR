---
title: 'IAMReX: an adaptive framework for the multiphase flow and fluid-particle interaction problems'
tags:
  - C++
  - Computational Fluid Dynamics
  - Adaptive Mesh Refinement
  - Immersed Boundary Method
  - Multiphase Flow
authors:
  - name: Chun Li
    orcid: 0009-0001-6081-5450
    affiliation: "1, 2"
  - name: Xuzhu Li
    orcid: 0009-0007-5348-4159
    affiliation: "1, 3"
  - name: Yiliang Wang
    orcid: 0000-0002-0156-2203
    affiliation: 4
  - name: Dewen Liu
    orcid: 0009-0006-2728-421X
    affiliation: 5
  - name: Shuai He
    orcid: 0009-0004-7754-8346
    affiliation: 6
  - name: Haoran Cheng
    orcid: 0009-0002-7636-038X
    affiliation: 7
  - name: Xiaokai Li
    orcid:  0009-0001-6639-0473
    affiliation: 8
  - name: Wenzhuo Li
    orcid: 0009-0001-7608-2992
    affiliation: 9
  - name: Mingze Tang
    orcid: 0009-0007-4194-9908
    affiliation: 10
  - name: Zhengping Zhu
    orcid: 0000-0002-1315-3554
    affiliation: 1
  - name: Yadong Zeng
    corresponding: true
    orcid: 0009-0001-7944-3597
    affiliation: 11
affiliations:
 - name: Research Center for Astronomical Computing, Zhejiang Laboratory, Hangzhou, 311100, China
   index: 1
 - name: School of Energy and Power Engineering, Lanzhou University of Technology, Lanzhou, 730050, China
   index: 2
 - name: School of Mechanical Engineering, Hefei University of Technology, Hefei, 230009, China
   index: 3
 - name: Independent Researcher, Loveland, 45140, USA
   index: 4
 - name: School of Civil Engineering, Southwest Jiaotong University, Chengdu, 611756, China 
   index: 5
 - name: School of Power and Energy, Northwestern Polytechnical University, Xi'an, 710129, China
   index: 6
 - name: Department of Electrical Engineering and Computer Science, University of Michigan, Ann Arbor, 48104, USA
   index: 7
 - name: School of Physical Science and Technology, ShanghaiTech University, Shanghai, 201210, China
   index: 8
 - name: Advanced Propulsion Laboratory, Department of Modern Mechanics, University of Science and Technology of China, Hefei, 230026, China
   index: 9
 - name: School of Aeronautics, Northwestern Polytechnical University, Xi'an, 710072, China
   index: 10
 - name: Department of Computer Science, University of Texas at Austin, Austin, 78712, USA
   index: 11

date: 17 January 2025
bibliography: paper.bib

# # Optional fields if submitting to a AAS journal too, see this blog post:
# # https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
# aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
# aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

IAMReX is an adaptive C++ solver designed for multiphase flow and fluid-particle interaction problems. It is built in a objected-oriented style and capable of high-performance massively parallel computing complex systems (e.g., gas-fluid interaction, cluster of particles).

The original goal of IAMReX is to extend the capability of IAMR codes [@almgren1998conservative], which only uses a density-based solver to capture the diffused interface of the two-phase flow. IAMReX offers the Level Set (LS) method and the reinitialization techniques for accurately capturing the two-phase interface [@zeng2022parallel], which increases the robustness of simulations with high Reynolds number [@zeng2023consistent]. For fluid-particle interaction problems, IAMReX employs the multidirect forcing immersed boundary method [@li2024open]. The associated Lagrangian markers used to resolve fluid-particle interface only exist on the finest-level grid, which greatly reduces memory cost. Both the subcycling and non-subcycling time advancement methods are implemented, and these methods help to decouple the time advancement at different levels. In addition, IAMReX is a publicly accessible platform designed specifically for developing massively parallel block-structured adaptive mesh refinement (BSAMR) applications. The code now supports hybrid parallelization using either pure MPI or MPI & OpenMP for multicore machines with the help of the AMReX framework [@zhang2019amrex].

The IAMReX code has undergone considerable development since 2023 and gained a few new contributors in the past two years. Although the projection-based flow solver is inherited from IAMR, IAMReX has added over 3,000 lines of new code, introduced 10 more new test cases, and contributed approximately 60 new commits on GitHub. The versatility, accuracy, and efficiency of the present IAMReX framework are demonstrated by simulating two-phase flow and fluid-particle interaction problems with various types of kinematic constraints. We carefully designed the document such that users can easily compile and run cases. Input files, profiling scripts, and raw postprocessing data are also available for reproducing all results.

# Statement of need

IAMReX is suitable for modeling multiphase flow problems and fluid-structure interaction problems. Its Level Set-based interface capturing technique can be beneficial for researchers studying phenomena such as wind over waves, breaking waves, and simulating the formation and disappearance of bubbles and droplets. Additionally, the immersed boundary method along with the collision models can parallelly resolve large-scale particles and capture their motions. Researchers working on studies of biological particle aggregation, sandstorms, wind erosion of ground surfaces, and seawater erosion of riverbeds are also among the target audience for this software.

# Acknowledgements

C.L., X.L., Y.Z., and Z.Z. are grateful to Ann Almgren, John Bell, Andy Nonaka, Candace Gilet, Andrew Myers, Axel Huebl, and Weiqun Zhang in the Lawrence Berkeley National Laboratory (LBNL) for their discussions related to AMReX and IAMR. We sincerely thank all the developers who have contributed to the original IAMR repository, although this paper focus on the newly added features and does not include the names of all contributors. Y.Z. and Z.Z. also thank Prof. Lian Shen, Prof. Ruifeng Hu, and Prof. Xiaojing Zheng during their Ph.D. studies.

# References
