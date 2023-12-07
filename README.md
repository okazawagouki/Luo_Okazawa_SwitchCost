# Data and code for Luo et al. (2023)

This repository contains data and test code for the following paper:

Luo T, Xu M, Zheng Z, Okazawa G (2023) Limitation of switching sensory information flow in flexible perceptual decision making.  bioRxiv 
[https://doi.org/10.1101/2023.12.03.569827](https://doi.org/10.1101/2023.12.03.569827)

## Requirement

The data are saved as MATLAB v7 mat file format.

The code is tested under MATLAB R2019b. Requires Statistics and Machine Learning Toolbox.

## Main dataset

[data.mat](./data/preproc/data.mat): Behavioral data collected from eight human participants in context-dependent face categorization task.

The dataset contains followings:

* cond (trial $\times$ 1): Task rule (either 1 or 2).
* cond_switch (trial $\times$ 1): 1 for non-switch trial, 2 for switch trial.
* csi (trial $\times$ 1): Cue-stimulus interval (s).
* fluc (trial $\times$ 1): Fluctuation of stimulus strength within each trial. In each cell, there is a matrix (number of frames $\times$ 6) which contains fluctuation of three facial features along two task axes.
* morph_level (trial $\times$ 2): Stimulus strength along two task axes.
* resp (trial $\times$ 1): Participants' choice (either 1 or 2).
* result (trial $\times$ 1): Result of a trial (either CORRECT, WRONG, NOFIX, FIXBREAK, or NOCHOICE).
* rt (trial $\times$ 1): Participants' reaction time (s).
* subj (trial $\times$ 1): Participants' ID.
* targ_cor (trial $\times$ 1): Correct target (either 1 or 2).

## Test code

[code.m](./code.m): The code will generate key figures of the paper. It will take a few seconds to generate the figures. Expected outcome of the figures is in the "figures" folder.

* Generate psychometric and chronometric functions averaged across participants.
* Generate psychophysical kernels.

## Contact

For questions or further inquiry about the dataset and code, please contact the first or lead author (luotl2021@ion.ac.cn, okazawa@ion.ac.cn).

## License

This material, including all data and code associated with this research, is protected under copyright. All rights are reserved by the authors. Upon the publication of the paper, the authors will distribute this material under a CC-BY 4.0 license.
