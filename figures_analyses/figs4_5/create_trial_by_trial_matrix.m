clear
clc
close all

%% read main task files

fileName1 = 'DesignMatrix_ZScoredAcrossSubjectWithNans_MainTaskBehavioralStudy.xls';
[data1,titles1,raw1] = xlsread(fileName1);
Subs1=data1(:,1);
Varian_OutSample=data1(:,17); % Varian OUT sample inconsistency index (alt criterion)
MMI_OutSample=data1(:,13); % MMI OUT sample inconsistency index (alt criterion)
Varian_InSample=data1(:,15); % Varian IN sample inconsistency index
MMI_InSample=data1(:,11); % MMI IN sample inconsistency index
Varian_OutSample_diff=data1(:,18); % Varian OUT sample inconsistency index (diff from full index)
MMI_OutSample_diff=data1(:,14); % MMI OUT sample inconsistency index (diff from full index)
slope1 = data1(:,4); % budget set's slope
intersection1 = data1(:,5); % budget set's intersection
X_allocation1 = data1(:,6); %X-axis choice
Y_allocation1 = data1(:,7); %Y-axis choice
RT1 = data1(:,10); %RT
SV = data1(:,19); 
ChoiceDifficulty = data1(:,20);
subjects1 = unique(data1(:,1));

Features1=data1(:,21:54);
features_titles = raw1(1,21:54);
num_features = size(Features1,2);


%% read motor task files
fileName2 = 'DesignMatrix_ZScoredAcrossSubject_MotorTask_BehavioralStudy.xls';
[data2,titles2,raw2] = xlsread(fileName2);
Subs2=data2(:,1);
slope2 = data2(:,4); % budget set's slope
intersection2 = data2(:,5); % budget set's intersection
X_allocation2 = data2(:,6); %X-axis choice
Y_allocation2 = data2(:,7); %Y-axis choice
RT2 = data2(:,10); %RT
subjects2 = unique(data2(:,1));

Features2=data2(:,11:44);

%% read numerical cognition task files
fileName3 = 'DesignMatrix_ZScoredAcrossSubject_NumericalTask_BehavioralStudy.xls';
[data3,titles3,raw3] = xlsread(fileName3);
Subs3=data3(:,1);
Features3=data3(:,11:44);
slope3 = data3(:,4); % budget set's slope
intersection3 = data3(:,5); % budget set's intersection
X_allocation3 = data3(:,6); %X-axis choice
Y_allocation3 = data3(:,7); %Y-axis choice
RT3 = data3(:,10); %RT
subjects3 = unique(data3(:,1));

all_subjects_across_tasks = [];

%% Find missing subjects - to relate trials from numerical and motor tasks to their trial number in the main task
missing_motor = ismember(subjects1,subjects2);
missing_numeric = ismember(subjects1,subjects3);

path_budget_sets_localizers =  [pwd '\BudgetSetsForLocalizers\'];
subjectlist = dir([path_budget_sets_localizers '*.xls']);


%% Subjects loop
% Set up the subject 
for i=1:length(subjects1)
 
    sid = subjects1(i);
    subject_missing_motor=0;
    subject_missing_numeric=0;
    if missing_motor(i)==0
        subject_missing_motor=1;
    end
    if missing_numeric(i)==0
        subject_missing_numeric=1;
    end
    
% upload the subject-specific file with the budget sets for the localizers    
    subject_budget_sets_localizers = xlsread([path_budget_sets_localizers subjectlist(i).name]);
    original_trial_number = subject_budget_sets_localizers(:,1);
    original_slope = subject_budget_sets_localizers(:,3);
        
% Extract subject specific rows in all tasks
    subject_rows1 = find(data1(:,1)==sid); %extract subject's own data in the main task
    subject_mouse_features1 = Features1(subject_rows1,:); %extract the subject's mouse features in the main task
    subject_Varian_OutSample=Varian_OutSample(subject_rows1); 
    subject_MMI_OutSample=MMI_OutSample(subject_rows1); 
    subject_Varian_InSample=Varian_InSample(subject_rows1); 
    subject_MMI_InSample=MMI_InSample(subject_rows1); 
    subject_Varian_OutSample_diff=Varian_OutSample_diff(subject_rows1); 
    subject_MMI_OutSample_diff=MMI_OutSample_diff(subject_rows1);
    subject_slope1 = slope1(subject_rows1); 
    subject_intersection1 = intersection1(subject_rows1); 
    subject_X_allocation1 = X_allocation1(subject_rows1); 
    subject_Y_allocation1 = Y_allocation1(subject_rows1);
    subject_RT_main_task = RT1(subject_rows1);
    subject_SV = SV(subject_rows1);
    subject_ChoiceDifficulty = ChoiceDifficulty(subject_rows1);
    
    if subject_missing_motor==0
        subject_rows2 = find(data2(:,1)==sid); %extract subject's own data in the main task
        subject_mouse_features2 = Features2(subject_rows2,:); %extract the subject's mouse features in the main task
        subject_slope2 = -slope2(subject_rows2); 
        subject_intersection2 = intersection2(subject_rows2); 
        subject_X_allocation2 = X_allocation2(subject_rows2); 
        subject_Y_allocation2 = Y_allocation2(subject_rows2);
        subject_RT_motor_task = RT2(subject_rows2);
    end
        
    if subject_missing_numeric==0
        subject_rows3 = find(data3(:,1)==sid); %extract subject's own data in the main task
        subject_mouse_features3 = Features3(subject_rows3,:); %extract the subject's mouse features in the main task
        subject_slope3 = -slope3(subject_rows3); 
        subject_intersection3 = intersection3(subject_rows3);  
        subject_X_allocation3 = X_allocation3(subject_rows3); 
        subject_Y_allocation3 = Y_allocation3(subject_rows3);
        subject_RT_numeric_task = RT3(subject_rows3);
    end
    
    
% Identify trial locations in each task
    %first case - we have both localizers
    if subject_missing_numeric==0 && subject_missing_motor==0
        for a=1:length(original_trial_number)
             localizer_slopes_locations_main_task(a) = original_trial_number(a);
             localizer_slopes_locations_motor_localizer(a) = find(subject_slope2(:)==original_slope(a));
             localizer_slopes_locations_numeric_localizer(a)= find(subject_slope3(:)==original_slope(a));
        end
        
   %second case - we don't have the numeric localizer
    elseif subject_missing_numeric==1 && subject_missing_motor==0
        for a=1:length(original_trial_number)
            localizer_slopes_locations_main_task(a) = original_trial_number(a);
            localizer_slopes_locations_motor_localizer(a)= find(subject_slope2(:)==original_slope(a));
        end
        
    %third case - we don't have the motor localizer
    elseif subject_missing_numeric==0 && subject_missing_motor==1
        for a=1:length(original_trial_number)
            localizer_slopes_locations_main_task(a) = original_trial_number(a);
            localizer_slopes_locations_numeric_localizer(a)= find(subject_slope3(:)==original_slope(a));
        end
        
    %forth case - we don't have both localizers (sub 202) 
    elseif subject_missing_numeric==1 && subject_missing_motor==1
        for a=1:length(original_trial_number)
            localizer_slopes_locations_main_task(a) = original_trial_number(a);
        end
    end
    
% Create the matrix 
    % basic trial parameters
    subject_across_tasks = zeros(length(original_trial_number),125);
    subject_across_tasks(:,1)  = sid; %subject number
    subject_across_tasks(:,2)  = original_trial_number; %original main task trial number
    if subject_missing_motor==0
        subject_across_tasks(:,3)  = localizer_slopes_locations_motor_localizer; %trial number in motor localizer
    else
        subject_across_tasks(:,3)  = nan;
    end
    if subject_missing_numeric==0
        subject_across_tasks(:,4)  = localizer_slopes_locations_numeric_localizer; %trial number in numeric localizer
    else
        subject_across_tasks(:,4)  = nan;
    end    
    subject_across_tasks(:,5)  = original_slope; %slope
    subject_across_tasks(:,6)  = subject_budget_sets_localizers(:,4); %intersection
    
    % choices (targets) and choices in localizers
    subject_across_tasks(:,7)  = subject_X_allocation1(localizer_slopes_locations_main_task); %X allocation/target X
    subject_across_tasks(:,8)  = subject_Y_allocation1(localizer_slopes_locations_main_task); %Y allocation/target Y
    if subject_missing_motor==0 
        subject_across_tasks(:,9)  = subject_X_allocation2(localizer_slopes_locations_motor_localizer);
        subject_across_tasks(:,10) = subject_Y_allocation2(localizer_slopes_locations_motor_localizer);
    else
        subject_across_tasks(:,9)  = nan;
        subject_across_tasks(:,10) = nan;
    end
    if subject_missing_numeric==0 
        subject_across_tasks(:,11) = subject_X_allocation3(localizer_slopes_locations_numeric_localizer);
        subject_across_tasks(:,12) = subject_Y_allocation3(localizer_slopes_locations_numeric_localizer);
    else
        subject_across_tasks(:,11) = nan;
        subject_across_tasks(:,12) = nan;
    end
    
    % inconsistency indices
    subject_across_tasks(:,13)  = subject_MMI_OutSample(localizer_slopes_locations_main_task); %MMI out sample (alt criterion)
    subject_across_tasks(:,14)  = subject_MMI_InSample(localizer_slopes_locations_main_task); %MMI in sample (alt criterion)
    subject_across_tasks(:,15)  = subject_MMI_OutSample_diff(localizer_slopes_locations_main_task); %MMI out sample (diff from full index)
    subject_across_tasks(:,16)  = subject_Varian_OutSample(localizer_slopes_locations_main_task); %Varian out sample
    subject_across_tasks(:,17)  = subject_Varian_InSample(localizer_slopes_locations_main_task); %Varian in sample
    subject_across_tasks(:,18)  = subject_Varian_OutSample_diff(localizer_slopes_locations_main_task); %Varian out sample (diff from full index)
    
        
    % RTs
    subject_across_tasks(:,19) = subject_RT_main_task(localizer_slopes_locations_main_task);
    if subject_missing_motor==0
        subject_across_tasks(:,20) = subject_RT_motor_task(localizer_slopes_locations_motor_localizer);
    else
        subject_across_tasks(:,20) = nan;
    end
    if subject_missing_numeric==0
        subject_across_tasks(:,21) = subject_RT_numeric_task(localizer_slopes_locations_numeric_localizer);
    else
        subject_across_tasks(:,21) = nan;
    end
    
    % SV and Choice difficulty
    subject_across_tasks(:,22) = subject_SV(localizer_slopes_locations_main_task);
    subject_across_tasks(:,23) = subject_ChoiceDifficulty(localizer_slopes_locations_main_task);
    
    % mouse features 
    subject_across_tasks(:,24:57)= subject_mouse_features1(localizer_slopes_locations_main_task,:);
    % find nan trials
    nan_main = find(isnan(subject_across_tasks(:,7)));
    subject_across_tasks(nan_main,24:57) = nan;
    if subject_missing_motor==0
        subject_across_tasks(:,58:91) = subject_mouse_features2(localizer_slopes_locations_motor_localizer,:);
        % find nan trials
        nan_motor = find(isnan(subject_across_tasks(:,9)));
        subject_across_tasks(nan_motor,58:91) = nan;
    else 
        subject_across_tasks(:,58:91) = nan;
    end
    if subject_missing_numeric==0
        subject_across_tasks(:,92:125) = subject_mouse_features3(localizer_slopes_locations_numeric_localizer,:); 
        % find nan trials
        nan_numerical = find(isnan(subject_across_tasks(:,11)));
        subject_across_tasks(nan_numerical,92:125) = nan;
    else
        subject_across_tasks(:,92:125) = nan;
    end
    
    all_subjects_across_tasks = [all_subjects_across_tasks; subject_across_tasks];
    
clear subject_rows1 subject_mouse_features1 nan_numerical nan_motor nan_main
clear subject_Varian_OutSample subject_Varian_InSample subject_MMI_OutSample subject_MMI_InSample
clear subject_slope1 subject_intersextion1 subject_intersextion1 subject_X_allocation subject_Y_allocation
clear subject_RT_main_task subject_SV subject_ChoiceDifficulty
clear subject_rows2 subject_mouse_features2 motor subject_slope2 
clear subject_intersextion2 subject_X_allocation2 subject_Y_allocation2 subject_RT_motor_task
clear subject_rows3 subject_mouse_features3 subject_slope3 
clear subject_intersextion3 subject_X_allocation3 subject_Y_allocation3 subject_RT_numeric_task
clear localizer_slopes_locations_main_task localizer_slopes_locations_numeric_localizer localizer_slopes_locations_motor_localizer
clear subject_across_tasks subject_budget_sets_localizers original_trial_number

end

%%

trial_by_trial_matrix = array2table(all_subjects_across_tasks, 'VariableNames',...
        {'SID', 'trial_number_main_task', 'trial_number_motor_task' ,'trial_number_numeric_task', ...
        'slope' ,'IntersectionY' ,'X_allocation', 'Y_allocation' ,'X_coordinate_motor' ,'Y_coordinate_motor' , ...
        'X_coordinate_numeric' ,'Y_coordinate_numeric', ...
        'MMI_Out_of_Sample_alt' ,'MMI_In_Sample_component', 'MMI_Out_of_Sample_diff', ...
        'Varian_Out_of_Sample_alt', 'Varian_In_Sample_component', 'Varian_Out_of_Sample_diff', ...
        'RT_main_task', 'RT_motor_task' ,'RT_numeric_task', 'SV', 'Choice_difficulty', ...   
        'xflips_main', 'yflips_main', 'MaxXflip_main', 'MaxYflip_main' ,'MaxProgress_main' ,'medianProgress_main', ...
        'Time2MaxProgress_main', 'meanVel_main', 'MaxVel_main', ...
        'medianVel_main', 'VelOver10_main', 'EndofTrialVel_main', 'MaxVelTime_main', ...
        'meanAcc_main', 'MaxAcc_main', 'N_Acc_main', ...
        'MaxArcLen_main,' 'Curvenes_main', 'Angle_STD_main', 'MD_main', 'AUC_main', ...
        'Layover_main',  'numFixations_main', 'MaxFixTime_main', 'MaxFixChoiceDist_main' ,... 
        'numBudgetFixations_main', 'AveFixTime_main', ...
        'AboveLine_main', 'XTimeOutofBounds_main', 'YTimeOutofBounds_main', ...
        'TimeNearPredicted_main' ,'TimeBudgetLine_main', ...
        'XSampEn4_main', 'YSampEn4_main', ... 
        'xflips_motor', 'yflips_motor', 'MaxXflip_motor', 'MaxYflip_motor', 'MaxProgress_motor', 'medianProgress_motor', ...
        'Time2MaxProgress_motor', 'meanVel_motor', 'MaxVel_motor', ...
        'medianVel_motor' , 'VelOver10_motor', 'EndofTrialVel_motor', 'MaxVelTime_motor', ...
        'meanAcc_motor', 'MaxAcc_motor' ,'N_Acc_motor' ,...
        'MaxArcLen_motor', 'Curvenes_motor', 'Angle_STD_motor', 'MD_motor', 'AUC_motor', ...
        'Layover_motor',  'numFixations_motor', 'MaxFixTime_motor', 'MaxFixChoiceDist_motor', ... 
        'numBudgetFixations_motor', 'AveFixTime_motor', ...
        'AboveLine_motor', 'XTimeOutofBounds_motor', 'YTimeOutofBounds_motor', ...
        'TimeNearPredicted_motor', 'TimeBudgetLine_motor', ...
        'XSampEn4_motor', 'YSampEn4_motor', ... 
        'xflips_numeric', 'yflips_numeric', 'MaxXflip_numeric', 'MaxYflip_numeric', 'MaxProgress_numeric', 'medianProgress_numeric', ...
        'Time2MaxProgress_numeric', 'meanVel_numeric', 'MaxVel_numeric', ...
        'medianVel_numeric', 'VelOver10_numeric', 'EndofTrialVel_numeric', 'MaxVelTime_numeric', ...
        'meanAcc_numeric', 'MaxAcc_numeric', 'N_Acc_numeric', ...
        'MaxArcLen_numeric', 'Curvenes_numeric', 'Angle_STD_numeric', 'MD_numeric', 'AUC_numeric', ...
        'Layover_numeric',  'numFixations_numeric', 'MaxFixTime_numeric', 'MaxFixChoiceDist_numeric', ... 
        'numBudgetFixations_numeric', 'AveFixTime_numeric', ...
        'AboveLine_numeric', 'XTimeOutofBounds_numeric', 'YTimeOutofBounds_numeric', ...
        'TimeNearPredicted_numeric', 'TimeBudgetLine_numeric', ...
        'XSampEn4_numeric', 'YSampEn _numeric'}); 
    
file_name = 'DesignMatrix_AcrossTasks_BehavioralStudy.xls';
writetable(trial_by_trial_matrix, file_name, 'WriteVariableNames', true);