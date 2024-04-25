# W-PIR-HT

This repository contains open-source code for the paper [Weakly Private Information Retrieval from Heterogeneously Trusted Servers](https://arxiv.org/abs/2402.17940).

## Usage
- auto_create_PIRtable ------------ Create PIR table and WPIR table
- auto_create_PIRtable_RK ------------ Create PIR table and WPIR table with random key output
- auto_compute_maxL_adjustp ------- Compute optimal MaxL with different p allocation
- auto_compute_MI_adjustp --------- Compute optimal MI with different p allocation
- plot_maxL_vs_D ------------------ Plot MaxL (leakage, download) region
- plot_MI_vs_D -------------------- Plot MI (leakage, download) region
- leakyPIR_LP_general ------------- linear programming in solving Leaky PIR 

For example, to reproduce experimental results of Max-L, run 'plot_maxL_vs_D.m'.

## Citing

If you find our work useful in your research, please cite the following paper:

> Zhao, W., Huang, Y. S., Zhou, R., & Tian, C. (2024). Weakly Private Information Retrieval from Heterogeneously Trusted Servers. arXiv preprint arXiv:2402.17940.

    @article{zhao2024weakly,
    title={Weakly Private Information Retrieval from Heterogeneously Trusted Servers},
    author={Zhao, Wenyuan and Huang, Yu Shin and Zhou, Ruida and Tian, Chao},
    journal={arXiv preprint arXiv:2402.17940},
    year={2024}
    }