clear
clc
close all


%% Read mouse features

fileName1 = 'DesignMatrix_AcrossTasks_BehavioralStudy.xls';
[data,titles,raw] = xlsread(fileName1);


%% general control variables
SID = data(:,1);
main_task_trial_num = data(:,2);
motor_task_trial_num = data(:,3);
numeric_task_trial_num = data(:,4);
slope = data(:,5);
IntersectY = data(:,6);
RT_main_task = data(:,19);
RT_motor_task = data(:,20);
RT_numeric_task = data(:,21);
SV = data(:,22);
ChoiceDifficulty = data(:,23);
MMI = data(:,14); 

controls_main_task = [SV RT_main_task slope IntersectY ChoiceDifficulty main_task_trial_num];
controls_motor_task = [SV RT_motor_task slope IntersectY ChoiceDifficulty motor_task_trial_num ];
controls_numerical_task = [SV RT_numeric_task slope IntersectY ChoiceDifficulty numeric_task_trial_num];

mouse_features_main_task = data(:,24:57);
mouse_features_motor_task = data(:,58:91);
mouse_features_numerical_task = data(:,92:125);

titles_for_elastic_net_main_task = titles(:,24:57);
titles_for_elastic_net_motor_task = titles(:,58:91);
titles_for_elastic_net_numerical_task = titles(:,92:125);

%% Run linear models just for the MAIN task 
% Elastic net
[PredictorsForRegression_main, PredictorsNamesForRegression_main, winningModelLambda_main, winningModelAlpha_main, winningModelMSE_main] = elastic_net(mouse_features_main_task,MMI,titles_for_elastic_net_main_task);

%%
% Regression
% no controls no random slope
X = PredictorsForRegression_main;
varsNamesRegression = [ PredictorsNamesForRegression_main {'MMI'}];  
lme_noControls_noRandomSlope_MainTask = fitlm(X,MMI, ...
            'VarNames',varsNamesRegression);
adjRsqr(1,1) = lme_noControls_noRandomSlope_MainTask.Rsquared.Adjusted;
% find significant features
sig_features = find(table2array(lme_noControls_noRandomSlope_MainTask.Coefficients(:,4))<0.05);
sig_features_names{1,1} = lme_noControls_noRandomSlope_MainTask.CoefficientNames(sig_features);
clear X varsNamesRegression sig_features

% no controls with random slope
X = [ones(length(PredictorsForRegression_main),1) PredictorsForRegression_main];
predictorsNamesRegression = [{'Intercept'} PredictorsNamesForRegression_main ];  
Z = ones(length(PredictorsForRegression_main),1);
G = SID;
lme_noControls_withRandomSlope_MainTask = fitlmematrix(X,MMI,Z,G);
adjRsqr(2,1) = lme_noControls_withRandomSlope_MainTask.Rsquared.Adjusted; 
% find significant features
sig_features = find(table2array(dataset2table(lme_noControls_withRandomSlope_MainTask.Coefficients(:,6)))<0.05);
sig_features_names{2,1} = predictorsNamesRegression(sig_features);
clear X Z G predictorsNamesRegression sig_features        

% with controls with random slope
X = [ones(length(PredictorsForRegression_main),1) PredictorsForRegression_main controls_main_task];
predictorsNamesRegression = [ {'Intercept'}  PredictorsNamesForRegression_main {'SV'} {'RT'} {'slope'} {'IntersectY'} {'ChoiceDifficulty'} {'main_task_trial_num'} ];  
Z = ones(length(PredictorsForRegression_main),1);
G = SID;
lme_withControls_withRandomSlope_MainTask = fitlmematrix(X,MMI,Z,G);
adjRsqr(3,1) = lme_withControls_withRandomSlope_MainTask.Rsquared.Adjusted;        
% find significant features
sig_features = find(table2array(dataset2table(lme_withControls_withRandomSlope_MainTask.Coefficients(:,6)))<0.05);
sig_features_names{3,1} = predictorsNamesRegression(sig_features);
clear X Z G predictorsNamesRegression sig_features 


%% Run linear models just for the MOTOR task 
% Elastic net
[PredictorsForRegression_motor, PredictorsNamesForRegression_motor, winningModelLambda_motor, winningModelAlpha_motor, winningModelMSE_motor] = elastic_net(mouse_features_motor_task,MMI,titles_for_elastic_net_motor_task);

%%
% Regression
% no controls no random slope
X = PredictorsForRegression_motor;
varsNamesRegression = [ PredictorsNamesForRegression_motor {'MMI'}];  
lme_noControls_noRandomSlope_MotorTask = fitlm(X,MMI, ...
            'VarNames',varsNamesRegression);
adjRsqr(1,2) = lme_noControls_noRandomSlope_MotorTask.Rsquared.Adjusted; 
% find significant features
sig_features = find(table2array(lme_noControls_noRandomSlope_MotorTask.Coefficients(:,4))<0.05);
sig_features_names{1,2} = lme_noControls_noRandomSlope_MotorTask.CoefficientNames(sig_features);
clear X varsNamesRegression sig_features

% no controls with random slope
X = [ones(length(PredictorsForRegression_motor),1) PredictorsForRegression_motor];
predictorsNamesRegression = [{'Intercept'} PredictorsNamesForRegression_motor ];  
Z = ones(length(PredictorsForRegression_motor),1);
G = SID;
lme_noControls_withRandomSlope_MotorTask = fitlmematrix(X,MMI,Z,G);
adjRsqr(2,2) = lme_noControls_withRandomSlope_MotorTask.Rsquared.Adjusted;        
% find significant features
sig_features = find(table2array(dataset2table(lme_noControls_withRandomSlope_MotorTask.Coefficients(:,6)))<0.05);
sig_features_names{2,2} = predictorsNamesRegression(sig_features);
clear X Z G predictorsNamesRegression sig_features        

% with controls with random slope
X = [ones(length(PredictorsForRegression_motor),1) PredictorsForRegression_motor controls_motor_task];
predictorsNamesRegression = [ {'Intercept'}  PredictorsNamesForRegression_motor {'SV'} {'RT'} {'slope'} {'IntersectY'} {'ChoiceDifficulty'} {'motor_task_trial_num'} ];  
Z = ones(length(PredictorsForRegression_motor),1);
G = SID;
lme_withControls_withRandomSlope_MotorTask = fitlmematrix(X,MMI,Z,G);
adjRsqr(3,2) = lme_withControls_withRandomSlope_MotorTask.Rsquared.Adjusted;        
% find significant features
sig_features = find(table2array(dataset2table(lme_withControls_withRandomSlope_MotorTask.Coefficients(:,6)))<0.05);
sig_features_names{3,2} = predictorsNamesRegression(sig_features);
clear X Z G predictorsNamesRegression sig_features


%% Run linear models just for the NUMERICAL task 
% Elastic net
[PredictorsForRegression_numerical, PredictorsNamesForRegression_numerical, winningModelLambda_numerical, winningModelAlpha_numerical, winningModelMSE_numerical] = elastic_net(mouse_features_numerical_task,MMI,titles_for_elastic_net_numerical_task);

%%
% Regression
% no controls no random slope
X = PredictorsForRegression_numerical;
varsNamesRegression = [ PredictorsNamesForRegression_numerical {'MMI'}];  
lme_noControls_noRandomSlope_NumericalTask = fitlm(X,MMI, ...
            'VarNames',varsNamesRegression);
adjRsqr(1,3) = lme_noControls_noRandomSlope_NumericalTask.Rsquared.Adjusted; 
% find significant features
sig_features = find(table2array(lme_noControls_noRandomSlope_NumericalTask.Coefficients(:,4))<0.05);
sig_features_names{1,3} = lme_noControls_noRandomSlope_NumericalTask.CoefficientNames(sig_features);
clear X varsNamesRegression sig_features

% no controls with random slope
X = [ones(length(PredictorsForRegression_numerical),1) PredictorsForRegression_numerical];
predictorsNamesRegression = [{'Intercept'} PredictorsNamesForRegression_numerical ];  
Z = ones(length(PredictorsForRegression_numerical),1);
G = SID;
lme_noControls_withRandomSlope_NumericalTask = fitlmematrix(X,MMI,Z,G);
adjRsqr(2,3) = lme_noControls_withRandomSlope_NumericalTask.Rsquared.Adjusted;        
% find significant features
sig_features = find(table2array(dataset2table(lme_noControls_withRandomSlope_NumericalTask.Coefficients(:,6)))<0.05);
sig_features_names{2,3} = predictorsNamesRegression(sig_features);
clear X Z G predictorsNamesRegression sig_features n        

% with controls with random slope
X = [ones(length(PredictorsForRegression_numerical),1) PredictorsForRegression_numerical controls_numerical_task];
predictorsNamesRegression = [ {'Intercept'}  PredictorsNamesForRegression_numerical {'SV'} {'RT'} {'slope'} {'IntersectY'} {'ChoiceDifficulty'} {'numerical_task_trial_num'} ];  
Z = ones(length(PredictorsForRegression_numerical),1);
G = SID;
lme_withControls_withRandomSlope_NumericalTask = fitlmematrix(X,MMI,Z,G);
adjRsqr(3,3) = lme_withControls_withRandomSlope_NumericalTask.Rsquared.Adjusted;    
% find significant features
sig_features = find(table2array(dataset2table(lme_withControls_withRandomSlope_NumericalTask.Coefficients(:,6)))<0.05);
sig_features_names{3,3} = predictorsNamesRegression(sig_features);
clear X Z G predictorsNamesRegression sig_features

%% output for figure - AdjRsqr 
adj_Rsqr_table = array2table(adjRsqr,'VariableNames',...
            {'Main task' 'Motor task' 'Numerical task'} , ...
           'RowNames', {'no random slope, no controls', 'with random slope, no controls','with random slope, with controls'}); 

adjRsqr_forFigure(1,1:3) = adjRsqr(1,1:3);
adjRsqr_forFigure(2,1:3) = adjRsqr(2,1:3)-adjRsqr(1,1:3);
adjRsqr_forFigure(3,1:3) = adjRsqr(3,1:3)-adjRsqr(2,1:3);
figure;
b11=bar(1:3,adjRsqr_forFigure','stacked','LineWidth',1,'Facecolor','flat');
b11(1).FaceColor = [255 197 197]./255;
b11(2).FaceColor = 'none';
hatchfill2(b11(2),'single','HatchAngle',45,'HatchDensity',60,'hatchcolor',[255 197 197]./255,'HatchLineWidth',1.75);
b11(3).FaceColor = 'none';
hatchfill2(b11(3),'single','HatchAngle',0,'HatchDensity',60,'hatchcolor',[255 197 197]./255,'HatchLineWidth',1.75);

hold on 

b11=bar(2:3,adjRsqr_forFigure(:,2:3)','stacked','LineWidth',1,'Facecolor','flat');
b11(1).FaceColor = [189 215 238]./255;
b11(2).FaceColor = 'none';
hatchfill2(b11(2),'single','HatchAngle',45,'HatchDensity',60,'hatchcolor',[189 215 238]./255,'HatchLineWidth',1.75);
b11(3).FaceColor = 'none';
hatchfill2(b11(3),'single','HatchAngle',0,'HatchDensity',60,'hatchcolor',[189 215 238]./255,'HatchLineWidth',1.75);

hold on

b11=bar(3,adjRsqr_forFigure(:,3)','stacked','LineWidth',1,'Facecolor','flat');
b11(1).FaceColor = [226 240 217]./255;
b11(2).FaceColor = 'none';
hatchfill2(b11(2),'single','HatchAngle',45,'HatchDensity',60,'hatchcolor',[226 240 217]./255,'HatchLineWidth',1.75);
b11(3).FaceColor = 'none';
hatchfill2(b11(3),'single','HatchAngle',0,'HatchDensity',60,'hatchcolor',[226 240 217]./255,'HatchLineWidth',1.75);

xticklabels({'Main task', 'Motor task', 'Numerical task'}); ylabel('adj-Rsqr'); xlim([0.25 3.75]);

        
%% out-of-sample prediction
% I use the betas from the regression from the main task on mouse features from the motor and numerical tasks

locations_regressors = ismember(titles_for_elastic_net_main_task,PredictorsNamesForRegression_main);
index_regressors = find(locations_regressors==1);

% MOTOR TASK
X_new_motor = [ones(length(PredictorsForRegression_motor),1) mouse_features_motor_task(:,locations_regressors) controls_motor_task];
Z_new = ones(length(PredictorsForRegression_motor),1);
G_new = SID;
%  prediction
yfit_motor_test = predict(lme_withControls_withRandomSlope_MainTask,X_new_motor, Z_new, G_new);

% NUMERICAL TASK
X_new_numerical = [ones(length(PredictorsForRegression_numerical),1) mouse_features_numerical_task(:,locations_regressors) controls_numerical_task];
Z_new = ones(length(PredictorsForRegression_numerical),1);
G_new = SID;
%  prediction
yfit_numerical_test = predict(lme_withControls_withRandomSlope_MainTask,X_new_numerical, Z_new, G_new);

% correlations between fitted data and MMI
[corr_motor_main, p_motor_main, RL_motor_main, RU_motor_main] = corrcoef(MMI,yfit_motor_test,'Rows','complete'); 
[corr_numerical_main, p_numerical_main, RL_numerical_main, RU_numerical_main] = corrcoef(MMI,yfit_numerical_test,'Rows','complete'); 
% spearman
[corr_motor_main2, p_motor_main2] = corr(MMI,yfit_motor_test,'Rows','complete','Type','Spearman'); 
[corr_numerical_main2, p_numerical_main2] = corr(MMI,yfit_numerical_test,'Rows','complete','Type','Spearman'); 

% output for paper - cross-task prediction - scatterplot predicted and actual MMI
figure; 
sizeMarker = 10;
ax1=subplot(1,2,1);
scatter(yfit_motor_test,MMI,sizeMarker,'MarkerEdgeColor',[189, 215, 238]./255,'MarkerFaceColor',[189, 215, 238]./255 ,'LineWidth',1);
h = lsline(ax1); h.Color = 'k'; h.LineWidth = 0.5; h.LineStyle = '-'; xlim([min(yfit_motor_test) max(yfit_motor_test)]); ylim([0 1]);
        
ax2=subplot(1,2,2);
scatter(yfit_numerical_test,MMI,sizeMarker,'MarkerEdgeColor',[197,224,180]./255,'MarkerFaceColor',[197,224,180]./255,'LineWidth',1);
h = lsline(ax2); h.Color = 'k'; h.LineWidth = 0.5; h.LineStyle = '-'; xlim([min(yfit_numerical_test) max(yfit_numerical_test)]); ylim([0 1]);
        

%% Regressions that solely include MD and meanVel  - the "winning" features

% MAIN TASK
MD_main = data(:,43);
meanVel_main = data(:,31);

% no controls no random slope
X = [MD_main meanVel_main ];
varsNamesRegression = [ {'MD_main'} {'meanVel_main'} {'MMI'}];  
lme_noControls_noRandomSlope_MD_meanVel_MainTask = fitlm(X,MMI, ...
            'VarNames',varsNamesRegression);
adjRsqr_MD_meanVel(1,1) = lme_noControls_noRandomSlope_MD_meanVel_MainTask.Rsquared.Adjusted; 
% find significant features
sig_features = find(table2array(lme_noControls_noRandomSlope_MD_meanVel_MainTask.Coefficients(:,4))<0.05);
sig_features_MD_meanVel_names{1,1} = lme_noControls_noRandomSlope_MD_meanVel_MainTask.CoefficientNames(sig_features);
clear X varsNamesRegression sig_features

% no controls with random slope
X = [ones(length(MD_main),1) MD_main meanVel_main ];
predictorsNamesRegression = [{'Intercept'} {'MD_main'} {'meanVel_main'} ];  
Z = ones(length(MD_main),1);
G = SID;
lme_noControls_withRandomSlope_MD_meanVel_MainTask = fitlmematrix(X,MMI,Z,G, ...
            'FixedEffectPredictors',predictorsNamesRegression,...
            'RandomEffectPredictors',{'Intercept'},'RandomEffectGroups',{'SID'});
adjRsqr_MD_meanVel(2,1) = lme_noControls_withRandomSlope_MD_meanVel_MainTask.Rsquared.Adjusted;        
clear X Z G predictorsNamesRegression               

% with controls with random slope
X = [ones(length(MD_main),1) MD_main meanVel_main controls_main_task];
predictorsNamesRegression = [ {'Intercept'}  {'MD_main'} {'meanVel_main'} {'SV'} {'RT'} {'slope'} {'IntersectY'} {'ChoiceDifficulty'} {'main_task_trial_num'} ];  
Z = ones(length(MD_main),1);
G = SID;
lme_withControls_withRandomSlope_MD_meanVel_MainTask = fitlmematrix(X,MMI,Z,G);
adjRsqr_MD_meanVel(3,1) = lme_withControls_withRandomSlope_MD_meanVel_MainTask.Rsquared.Adjusted;        
clear X Z G predictorsNamesRegression 


% MOTOR TASK
% no controls no random slope
MD_motor = data(:,77);
meanVel_motor = data(:,65);

% no controls no random slope
X = [MD_motor meanVel_motor ];
varsNamesRegression = [ {'MD_motor'} {'meanVel_motor'} {'MMI'}];  
lme_noControls_noRandomSlope_MD_meanVel_MotorTask = fitlm(X,MMI, ...
            'VarNames',varsNamesRegression);
adjRsqr_MD_meanVel(1,2) = lme_noControls_noRandomSlope_MD_meanVel_MotorTask.Rsquared.Adjusted; 
% find significant features
sig_features = find(table2array(lme_noControls_noRandomSlope_MD_meanVel_MotorTask.Coefficients(:,4))<0.05);
sig_features_MD_meanVel_names{1,2} = lme_noControls_noRandomSlope_MD_meanVel_MotorTask.CoefficientNames(sig_features);
clear X varsNamesRegression sig_features

% no controls with random slope
X = [ones(length(MD_motor),1) MD_motor meanVel_motor ];
predictorsNamesRegression = [{'Intercept'} {'MD_motor'} {'meanVel_motor'} ];  
Z = ones(length(MD_motor),1);
G = SID;
lme_noControls_withRandomSlope_MD_meanVel_MotorTask = fitlmematrix(X,MMI,Z,G, ...
            'FixedEffectPredictors',predictorsNamesRegression,...
            'RandomEffectPredictors',{'Intercept'},'RandomEffectGroups',{'SID'});
adjRsqr_MD_meanVel(2,2) = lme_noControls_withRandomSlope_MD_meanVel_MotorTask.Rsquared.Adjusted;        
clear X Z G predictorsNamesRegression               

% with controls with random slope
X = [ones(length(MD_motor),1) MD_motor meanVel_motor controls_motor_task];
predictorsNamesRegression = [ {'Intercept'}  {'MD_motor'} {'meanVel_motor'} {'SV'} {'RT'} {'slope'} {'IntersectY'} {'ChoiceDifficulty'} {'motor_task_trial_num'} ];  
Z = ones(length(MD_motor),1);
G = SID;
lme_withControls_withRandomSlope_MD_meanVel_MotorTask = fitlmematrix(X,MMI,Z,G);
adjRsqr_MD_meanVel(3,2) = lme_withControls_withRandomSlope_MD_meanVel_MainTask.Rsquared.Adjusted;        
clear X Z G predictorsNamesRegression 


% NUMERICAL TASK
% no controls no random slope
MD_numerical = data(:,111);
meanVel_numerical = data(:,99);

% no controls no random slope
X = [MD_numerical meanVel_numerical ];
varsNamesRegression = [ {'MD_numerical'} {'meanVel_numerical'} {'MMI'}];  
lme_noControls_noRandomSlope_MD_meanVel_NumericalTask = fitlm(X,MMI, ...
            'VarNames',varsNamesRegression);
adjRsqr_MD_meanVel(1,3) = lme_noControls_noRandomSlope_MD_meanVel_NumericalTask.Rsquared.Adjusted; 
% find significant features
sig_features = find(table2array(lme_noControls_noRandomSlope_MD_meanVel_NumericalTask.Coefficients(:,4))<0.05);
sig_features_MD_meanVel_names{1,3} = lme_noControls_noRandomSlope_MD_meanVel_NumericalTask.CoefficientNames(sig_features);
clear X varsNamesRegression sig_features

% no controls with random slope
X = [ones(length(MD_numerical),1) MD_numerical meanVel_numerical ];
predictorsNamesRegression = [{'Intercept'} {'MD_numerical'} {'meanVel_numerical'} ];  
Z = ones(length(MD_numerical),1);
G = SID;
lme_noControls_withRandomSlope_MD_meanVel_NumericalTask = fitlmematrix(X,MMI,Z,G, ...
            'FixedEffectPredictors',predictorsNamesRegression,...
            'RandomEffectPredictors',{'Intercept'},'RandomEffectGroups',{'SID'});
adjRsqr_MD_meanVel(2,3) = lme_noControls_withRandomSlope_MD_meanVel_NumericalTask.Rsquared.Adjusted;        
clear X Z G predictorsNamesRegression               

% with controls with random slope
X = [ones(length(MD_numerical),1) MD_numerical meanVel_numerical  controls_numerical_task];
predictorsNamesRegression = [ {'Intercept'}  {'MD_numerical'} {'meanVel_numerical'} {'SV'} {'RT'} {'slope'} {'IntersectY'} {'ChoiceDifficulty'} {'numerical_task_trial_num'} ];  
Z = ones(length(MD_numerical),1);
G = SID;
lme_withControls_withRandomSlope_MD_meanVel_NumericalTask = fitlmematrix(X,MMI,Z,G);
adjRsqr_MD_meanVel(3,3) = lme_withControls_withRandomSlope_MD_meanVel_NumericalTask.Rsquared.Adjusted;        
clear X Z G predictorsNamesRegression 


