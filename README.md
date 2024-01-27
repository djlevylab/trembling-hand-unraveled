# The trembling hand unraveled: the motor dynamics and neural sources of choice inconsistency
This is the code package for the manuscript “The trembling hand unraveled: the motor dynamics and neural sources of choice inconsistency” by Vered Kurtz-David, Asaf Madar, Adam Hakim, Noa Palmon & Dino J Levy.
<br>The manuscript presents the results from three studies (behavioral, MRI and replication studies). Only the behavioral study codes are currently shared in the repository. 
## System requirements
The code was run and tested on MATLAB 2021b edition, Windows 11.
## Installation guide
Install MATLAB 2020b, 2021b, or 2022b editions. 
<br>Install the Deep Learning Toolbox and the Symbolic Math Toolbox.
## Demo
To extract the mouse features from mouse tracking data, run the script located at extract_mouse_features\main.m using MATLAB’s interface.
- Input: see directories extract_mouse_features\data\Main task\Trial_level_data_and_parameters, extract_mouse_features\data\Main task\Mouse_tracking
Each subject has 6 xls files, which contain their trial level mouse tracking data, and trial level inconsistency data, for all three tasks: main task, motor task and numerical task (3 tasks * 2 files per task).
- Output: AllResults-Main task.mat, AllResults-Motor_localizer.mat, AllResults-Numerical_localizer.mat. These are mouse features structs stored in three files corresponding to the three tasks used in the manuscript, and are required for the next sections.
<br>Runtime: 35-40 minutes.
### Inconsistency indices
To generate the inconsistency indices download the code package provided by Halevey et al. (2018) from https://github.com/persitzd/RP-Toolkit 
<br>The output of the trial-level inconsistency indices is located in extract_mouse_features\data\Main task\Trial_level_data_and_parameters.
<br>The output of the subject-level inconsistency indices is located in figures_analyses\fig2\ subject_level_inconsistency_indices.xls
<br>The raw data for subjects’ choices is in directory raw_data_for_inconsistency_indices. All raw data files have the following structure:
1) Column A - subject ID
2) Column B - trial number
3) Column C - choice X-axis
4) Column D - choice Y-axis
5) Column E - max. X-axis (intersection of the budget line with the X-axis)
6) Column F - max. Y-axis (intersection of the budget line with the Y-axis)
7) Coulumn E+F provide the budget line parameters

### Figure 2
1) For plotting subjects’ inconsistency indices, run the code in directory fig2: Fig2a_plottingInconsistencyIndices.m.
Input: The code uses the subject-level inconsistency indices (output from the Halevy et al. code package) for all three studies. 
(*) This code also outputs Supplementary Fig. 7.
2) For plotting the share of total expenditure that subjects allocated to the Y-axis as a function of the budget line’s slope, run Fig2b_PlottingSubjectsChoicesLogPrices_aggregateAllSubs.m. 
Input: The code uses subjects’ choices raw data from all three studies.
3) For plotting the distributions of Euclidean Distances from targets in the non-value tasks, run Fig2c_EucDist_distributions.m
Input: The code uses as input the design matrices from the cross-tasks regressions. See ### Figures 4 & 5 below for details about generating those matrices (for the Behavioral study). 
### Figure 3
For plotting subjects’ mouse features distributions run the code in directory fig3: Fig3e_features_dist_across_trials_and_subjects_BehavioralStudy.m. 
- Input: AllResults-Main task.mat (output of extract_mouse_features\main.m)
### Figure 4 & 5
For regression analyses, run the following codes:
1) Create regression design matrices for each task (which also z-score features values), using: writeDesignMatrix_MainTask_BehavioralStudy.m, writeDesignMatrix_MotorTask_BehavioralStudy.m, writeDesignMatrix_NumericalTask_BehavioralStudy.m
- Inputs: AllResults-Main task.mat, AllResults-Motor_localizer.mat, AllResults-Numerical_localizer.mat
- Output: DesignMatrix_ZScoredAcrossSubject_MotorTask_BehavioralStudy.xls, DesignMatrix_ZScoredAcrossSubject_NumericalTask_BehavioralStudy.xls, DesignMatrix_ZScoredAcrossSubjectWithNans_MainTaskBehavioralStudy.xls
2) Combine those design matrices into one design matrix for all the cross-tasks regressions, using create_trial_by_trial_matrix.m
- Output: DesignMatrix_AcrossTasks_BehavioralStudy.xls
3) Run regressions that only include the main task, using Regressions_onlyMainTask_SelectedFeatures_BehavioralStudy.m
- Input: DesignMatrix_ZScoredAcrossSubjectWithNans_MainTaskBehavioralStudy.xls
- Apart from the regression tables, this code also output Fig. 4e (Behavioral study) and Supplementary Fig. 9.
4) Run the cross-tasks regressions, using Regressions_AllTasksCombined_SelectedFeatures_BehavioralStudy.m
- Input: DesignMatrix_AcrossTasks_BehavioralStudy.xls
- Apart from the regression tables, this code also output Fig. 4d (Behavioral study) and Fig. 5b.
