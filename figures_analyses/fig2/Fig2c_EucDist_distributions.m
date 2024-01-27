clear
clc
close all

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code creates a figure displaying the distributions of Euclidean
% distances from targets in the non-value additional tasks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% color palette
main_task_color = [255,197,197]./255;
motor_task_color = [189,215,238]./255;
numerical_task_color = [226,240,217]./255;

%% Behavioral study

fileName1 = 'DesignMatrix_AcrossTasks_BehavioralStudy.xls';
[data,titles,raw] = xlsread(fileName1);

% Calculate Euc Distance motor 
X_target_behavioral = data(:,7);
X_MousePressMotor_behavioral = data(:,9);
X_MousePressNumerical_behavioral = data(:,11);

Y_target_behavioral = data(:,8);
Y_MousePressMotor_behavioral = data(:,10);
Y_MousePressNumerical_behavioral = data(:,12);

EucDist_MotorTask_BehavioralStudy = sqrt((X_target_behavioral-X_MousePressMotor_behavioral).^2+(Y_target_behavioral-Y_MousePressMotor_behavioral).^2);
EucDist_NumericalTask_BehavioralStudy = sqrt((X_target_behavioral-X_MousePressNumerical_behavioral).^2+(Y_target_behavioral-Y_MousePressNumerical_behavioral).^2);

% Descriptive stats
mean_EucDist_Motor_Behavioral = nanmean(EucDist_MotorTask_BehavioralStudy);
std_EucDist_Motor_Behavioral = nanstd(EucDist_MotorTask_BehavioralStudy);
min__EucDist_Motor_Behavioral = min(EucDist_MotorTask_BehavioralStudy);
max__EucDist_Motor_Behavioral = max(EucDist_MotorTask_BehavioralStudy);
 
mean_EucDist_Numeric_Behavioral = nanmean(EucDist_NumericalTask_BehavioralStudy);
std_EucDist_Numeric_Behavioral = nanstd(EucDist_NumericalTask_BehavioralStudy);
min__EucDist_Numeric_Behavioral = min(EucDist_NumericalTask_BehavioralStudy);
max__EucDist_Numeric_Behavioral = max(EucDist_NumericalTask_BehavioralStudy);
 

%% Neuroimaging study

fileName2 = 'DesignMatrix_AcrossTasks_NeuroimagingStudy.xls';
[data2,titles2,raw2] = xlsread(fileName2);

% Calculate Euc Distance 
X_target = data2(:,6);
X_MousePress = data2(:,8);

Y_target = data2(:,7);
Y_MousePress = data2(:,9);

EucDist_MotorTask_NeuroimagingStudy = sqrt((X_target-X_MousePress).^2+(Y_target-Y_MousePress).^2);

% Descriptive stats
mean_EucDist_Motor_Neuroimaging = nanmean(EucDist_MotorTask_NeuroimagingStudy);
std_EucDist_Motor_Neuroimaging = nanstd(EucDist_MotorTask_NeuroimagingStudy);
min__EucDist_Motor_Neuroimaging = min(EucDist_MotorTask_NeuroimagingStudy);
max__EucDist_Motor_Neuroimaging = max(EucDist_MotorTask_NeuroimagingStudy);

%% Replication study
fileName3 = 'DesignMatrix_AcrossTasks_ReplicationStudy.xls';
[data3,titles3,raw3] = xlsread(fileName3);

% Calculate Euc Distance motor 
X_target_replication = data3(:,7);
X_MousePressMotor_replication = data3(:,9);
X_MousePressNumerical_replication = data3(:,11);

Y_targetreplication = data3(:,8);
Y_MousePressMotor_replication = data3(:,10);
Y_MousePressNumerical_replication = data3(:,12);

EucDist_MotorTask_ReplicationStudy = sqrt((X_target_replication-X_MousePressMotor_replication).^2+(Y_targetreplication-Y_MousePressMotor_replication).^2);
EucDist_NumericalTask_ReplicationStudy = sqrt((X_target_replication-X_MousePressNumerical_replication).^2+(Y_targetreplication-Y_MousePressNumerical_replication).^2);

% Descriptive stats
mean_EucDist_Motor_Replication = nanmean(EucDist_MotorTask_ReplicationStudy);
std_EucDist_Motor_Replication = nanstd(EucDist_MotorTask_ReplicationStudy);
min__EucDist_Motor_Replication = min(EucDist_MotorTask_ReplicationStudy);
max__EucDist_Motor_Replication = max(EucDist_MotorTask_ReplicationStudy);
 
mean_EucDist_Numeric_Replication = nanmean(EucDist_NumericalTask_ReplicationStudy);
std_EucDist_Numeric_Replication = nanstd(EucDist_NumericalTask_ReplicationStudy);
min__EucDist_Numeric_Replication = min(EucDist_NumericalTask_ReplicationStudy);
max__EucDist_Numeric_Replication = max(EucDist_NumericalTask_ReplicationStudy);

%% Figure
figure;
edges = 0:0.1:5;
h1 = histogram(EucDist_MotorTask_BehavioralStudy,edges,'Normalization', 'probability');
h1.FaceColor = motor_task_color; h1.FaceAlpha = 1.0; 
hold on 
h4 = histogram(EucDist_NumericalTask_BehavioralStudy,edges,'Normalization', 'probability');
h4.FaceColor = numerical_task_color; h4.FaceAlpha = 1.0; 

ylim([0 0.4]);
title('Behavioral study');

figure;
h2 = histogram(EucDist_MotorTask_NeuroimagingStudy,edges,'Normalization', 'probability');
h2.FaceColor = motor_task_color; h2.FaceAlpha = 1.0; 
ylim([0 0.4]);
% legend('Euc. dist. motor task');
 title('Neuroimaging study');

figure;
h3 = histogram(EucDist_MotorTask_ReplicationStudy,edges,'Normalization', 'probability');
h3.FaceColor = motor_task_color; h3.FaceAlpha = 1.0; 
hold on
h5 = histogram(EucDist_NumericalTask_ReplicationStudy,edges,'Normalization', 'probability');
h5.FaceColor = numerical_task_color; h5.FaceAlpha = 1.0; 
ylim([0 0.4]);
title('Replication study');

%% statistical tests (Wilcoxon rans-sum test, due to the skewed distributions)
p1 = ranksum(EucDist_MotorTask_BehavioralStudy,EucDist_MotorTask_NeuroimagingStudy);
p2 = ranksum(EucDist_MotorTask_BehavioralStudy,EucDist_MotorTask_ReplicationStudy);
p3 = ranksum(EucDist_NumericalTask_BehavioralStudy,EucDist_NumericalTask_ReplicationStudy);
