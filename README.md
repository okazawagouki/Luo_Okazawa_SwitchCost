# Data and code for Luo et al.

This repository contains data and test code for the following paper:

Luo T, Xu M, Zheng Z, Okazawa G (2023) Limitation of switching sensory information flow in flexible perceptual decision making.  bioRxiv 
[https://doi.org/10.1101/2023.12.03.569827](https://doi.org/10.1101/2023.12.03.569827)

## Requirement

The data are saved as MATLAB v7 mat file format.

The code is tested under MATLAB R2019b on Windows. Requires Statistics and Machine Learning Toolbox.

## Dataset

[context_dependent_face_categorization_task.mat](./data/preproc/context_dependent_face_categorization_task.mat): Behavioral data collected from eight human participants in context-dependent face categorization task. Reported in Fig. 2-4. This dataset contains the following basic task parameters:

* cond: Task rule (either 1 or 2).
* cond_switch: 1 for non-switch trial, 2 for switch trial.
* csi: Cue-stimulus interval (s).
* fluc: Fluctuation of stimulus strength within each trial in the main context dependent face categorization task (identity/age vs. expression). The matrix (number of frames $\times$ 6) contains fluctuation of three facial features along two task axes.
* morph_level: Stimulus strength along two task axes.
* resp: Participants' choice (either 1 or 2).
* result: Result of a trial.
* rt: Participants' reaction time (s).
* subj: Participants' ID.
* targ_cor: Correct target (either 1 or 2).

[fixed_stimulus_duration_task.mat](./data/preproc/fixed_stimulus_duration_task.mat): Seven human participants in total. Reported in Fig. 5. This dataset contains the basic task parameters above.

[identity_versus_color.mat](./data/preproc/identity_versus_color.mat): Seven human participants in total. Reported in Fig. 6. Besides the basic task paramters, this dataset contains:

* fluc_face: Fluctuation of stimulus strength along identity axis within each trial. The matrix (number of frames $\times$ 3) contains fluctuation of three facial features.
* fluc_color: Fluctuation along color axis.

[motion_versus_color.mat](./data/preproc/motion_versus_color.mat): Seven human participants in total. Reported in Fig. 6. Besides the basic task paramters, this dataset contains:

* motion_energy: The quantification of the fluctuating motion coherence of the random dots.
* fluc_color: The quantification of the fluctuating color coherence of the random dots.


## Code

[fig2.m](./fig2.mat): The code will generate figures in Fig. 2. It will take a few seconds to generate the figures.

Similarly, fig*.m will generate the corresponding figures reported in the paper. Figures related to model fitting result are more time-consuming, which may take a few minutes on a regular personal computer.


## Contact

For questions or further inquiry about the dataset and code, please contact the first or lead author (luotl2021@ion.ac.cn, okazawa@ion.ac.cn).

## License

This material, including all data and code associated with this research, is protected under copyright. All rights are reserved by the authors. Upon the publication of the paper, the authors will distribute this material under a CC-BY 4.0 license.
