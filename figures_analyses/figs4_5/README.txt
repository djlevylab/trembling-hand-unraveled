
Inputs:
----- AllResults-Main task.mat (for the main task)
----- AllResults-Motor_localizer.mat (for the motor task)
----- AllResults-Numerical_localizer.mat (for the numerical task)

These structs are the output from the mouse features code package.

To generate Figures 4+5:

1) Create regression design matrices for each task (which also z-score features values), using:
----- writeDesignMatrix_MainTask_BehavioralStudy.m
----- writeDesignMatrix_MotorTask_BehavioralStudy.m
----- writeDesignMatrix_NumericalTask_BehavioralStudy.m

2) Combine those design matrices into one design matrix for all the cross-tasks regressions, using:
------ create_trial_by_trial_matrix.m

3) Run regressions that only include the main task, using:
------ Regressions_onlyMainTask_SelectedFeatures_BehavioralStudy.m

4) Run the cross-tasks regressions, using:
Regressions_AllTasksCombined_SelectedFeatures_BehavioralStudy.m


