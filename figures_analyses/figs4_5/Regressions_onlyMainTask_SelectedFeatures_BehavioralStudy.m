clear
clc
close all

% BEHAVIORAL STUDY - Main task - selected features
% load design matrix
fileName1 = 'DesignMatrix_ZScoredAcrossSubjectWithNans_MainTaskBehavioralStudy.xls';
[data,titles,raw] = xlsread(fileName1);


%% variables for regressions
SID = data(:,1);
main_task_trial_num = data(:,2);
block_num = data(:,3);
slope = data(:,4);
IntersectY = data(:,5);
RT = data(:,10);
SV = data(:,19);
ChoiceDifficulty = data(:,20);
MMI = data(:,11);

mouse_features = data(:,21:54);
controls = [SV RT slope IntersectY ChoiceDifficulty main_task_trial_num];
titles_for_elastic_net =  titles(1,21:54); 

%% Elastic net
[PredictorsForRegression, PredictorsNamesForRegression, winningModelLambda, winningModelAlpha, winningModelMSE] = elastic_net(mouse_features,MMI,titles_for_elastic_net);

%% Regression
% no controls no random slope
X =  PredictorsForRegression;
varsNamesRegression = [ PredictorsNamesForRegression {'MMI'}];  
lme_noControls_noRandomSlope_MainTask = fitlm(X,MMI, ...
            'VarNames',varsNamesRegression);
adjRsqr(1,1) = lme_noControls_noRandomSlope_MainTask.Rsquared.Adjusted; 
% find significant features
sig_features = find(table2array(lme_noControls_noRandomSlope_MainTask.Coefficients(:,4))<0.05);
sig_features_names{1,1} = lme_noControls_noRandomSlope_MainTask.CoefficientNames(sig_features);
clear X varsNamesRegression sig_features

% no controls with random slope
X = [ones(length(PredictorsForRegression),1) PredictorsForRegression];
predictorsNamesRegression = [ {'Intercept'}  PredictorsNamesForRegression];  
Z = ones(length(PredictorsForRegression),1);
G = SID;
lme_noControls_withRandomSlope_MainTask = fitlmematrix(X,MMI,Z,G, ...
            'FixedEffectPredictors',predictorsNamesRegression,...
            'RandomEffectPredictors',{'Intercept'},'RandomEffectGroups',{'SID'});
adjRsqr(2,1) = lme_noControls_withRandomSlope_MainTask.Rsquared.Adjusted;  
% find significant features
sig_features = find(table2array(dataset2table(lme_noControls_withRandomSlope_MainTask.Coefficients(:,6)))<0.05);
sig_features_names{2,1} = predictorsNamesRegression(sig_features);
clear X Z G predictorsNamesRegression sig_features       

% with controls with random slope
X = [ones(length(PredictorsForRegression),1) PredictorsForRegression controls];
predictorsNamesRegression = [ {'Intercept'}  PredictorsNamesForRegression {'SV'} {'RT'} {'slope'} {'IntersectY'} {'ChoiceDifficulty'} {'main_task_trial_num'} ];  
Z = ones(length(PredictorsForRegression),1);
G = SID;
lme_withControls_withRandomSlope_MainTask = fitlmematrix(X,MMI,Z,G, ...
            'FixedEffectPredictors',predictorsNamesRegression,...
            'RandomEffectPredictors',{'Intercept'},'RandomEffectGroups',{'SID'});
adjRsqr(3,1) = lme_withControls_withRandomSlope_MainTask.Rsquared.Adjusted;  
% find significant features
sig_features = find(table2array(dataset2table(lme_withControls_withRandomSlope_MainTask.Coefficients(:,6)))<0.05);
sig_features_names{3,1} = predictorsNamesRegression(sig_features);
clear X Z G predictorsNamesRegression sig_features  

%% Regress SV on motor features from main task
% Regression
% no controls no random slope
X = PredictorsForRegression;
varsNamesRegression = [ PredictorsNamesForRegression {'SV'}];  
lme_noControls_noRandomSlope_MainTask_SV = fitlm(X,SV, ...
            'VarNames',varsNamesRegression);
adjRsqr_SV(1,1) = lme_noControls_noRandomSlope_MainTask_SV.Rsquared.Adjusted;
clear X varsNamesRegression 

% no controls with random slope
X = [ones(length(PredictorsForRegression),1) PredictorsForRegression];
predictorsNamesRegression = [{'Intercept'} PredictorsNamesForRegression ];  
Z = ones(length(PredictorsForRegression),1);
G = SID;
lme_noControls_withRandomSlope_MainTask_SV = fitlmematrix(X,SV,Z,G);
adjRsqr_SV(2,1) = lme_noControls_withRandomSlope_MainTask_SV.Rsquared.Adjusted; 
clear X Z G predictorsNamesRegression 

%% Leave one subject out - selected features
subjects = unique(SID);

for sub=1:length(subjects)
    
   % find sub_rows
   subNum = subjects(sub);
   sub_rows = find(SID==subNum);
   allOther_subjects_rows = find(SID~=subNum);
   current_subject_matrix = mouse_features(sub_rows,:);
   controls_subject_matrix = controls(sub_rows,:);
   MMI_subject = MMI(sub_rows,1);
   allOther_subjects_matrix = mouse_features(allOther_subjects_rows,:);
   controls_allOther_subjects_matrix = controls(allOther_subjects_rows,:);
   MMI_allOther_subjects = MMI(allOther_subjects_rows,1);
   
   % run elastic net on all other subjects
   [PredictorsForRegression, PredictorsNamesForRegression, ~, ~, ~] = elastic_net(allOther_subjects_matrix,MMI_allOther_subjects,titles_for_elastic_net);
   
   % run linear model based on all subjects except the current subject 
   X = [ones(length(PredictorsForRegression),1) PredictorsForRegression controls_allOther_subjects_matrix];
   predictorsNamesRegression = [ {'Intercept'}  PredictorsNamesForRegression {'SV'} {'RT'} {'slope'} {'IntersectY'} {'ChoiceDifficulty'} {'main_task_trial_num'} ];  
   Z = ones(length(PredictorsForRegression),1);
   G = SID(allOther_subjects_rows);
   lme_only_main_task_selectedFeatures_with_controls_leaveOneOut = fitlmematrix(X,MMI_allOther_subjects,Z,G, ...
                'FixedEffectPredictors',predictorsNamesRegression,...
                'RandomEffectPredictors',{{'Intercept'}},'RandomEffectGroups',{'SID'});       
   
   % test set - the current subject
   locations_regressors = ismember(titles_for_elastic_net,PredictorsNamesForRegression);
   index_regressors = find(locations_regressors==1);
   X_new = [ones(length(sub_rows),1) current_subject_matrix(:,index_regressors) controls_subject_matrix];
   Z_new = ones(length(sub_rows),1);
   G_new = SID(sub_rows);
   
   ypredict_current_sub = predict(lme_only_main_task_selectedFeatures_with_controls_leaveOneOut,X_new, Z_new, G_new);
   
   [corrSub, pSub] = corrcoef(MMI_subject,ypredict_current_sub,'Rows','complete'); 
   corr_leaveOneSubOut(sub,1) = corrSub(2,1);
   pval_leaveOneSubOut(sub,1) = pSub(2,1);
   
   clear sub_rows current_subject_matrix allOther_subjects_matrix allOther_subjects_rows
   clear controls_subject_matrix controls_allOther_subjects_matrix MMI_subject MMI_allOther_subjects
   clear ypredict_current_sub lme_only_main_task_selectedFeatures_with_controls_leaveOneOut corrSub pSub
   clear X Z G predictorsNamesRegression PredictorsForRegression PredictorsNamesForRegression
   clear locations_regressors index_regressors X_new Z_new G_new
   
end

% ttest
[h1,p1] = ttest(corr_leaveOneSubOut(:,1),0,'Tail','right');

% figure
group = [ones(length(corr_leaveOneSubOut))];

xCenter = size(corr_leaveOneSubOut,2); 
spread = 0.4; % 0=no spread; 0.5=random spread within box bounds (can be any value)

figure()
[h,~,MX,MED,bw]=violin(corr_leaveOneSubOut, 'mc',[], 'medc','--k', ...
    'facealpha',0.5,'facecolor', [255 197 197]./255,'edgecolor',[192 0 0]./255);
hold on   
scatter1 = plot(rand(size(corr_leaveOneSubOut))*spread -(spread/2) + xCenter, corr_leaveOneSubOut, ...
            'o','linewidth', 1,'MarkerEdgeColor',[192 0 0]./255);
ylabel('correlation (rho)'); 
ylim([-0.1 1]);

[h1,p1,ci1,stats1] = ttest(corr_leaveOneSubOut,0,'Tail','right');



%% Motor learning effects during the task 
% Compare results from regression from the fisrt half of the task to the second half

% first half of the task
First_half = find(main_task_trial_num<=75);

% Elastic net - first half
[PredictorsForRegression_first_half, PredictorsNamesForRegression_first_half, ~, ~, ~] = elastic_net(mouse_features(First_half,:),MMI(First_half,:),titles_for_elastic_net);

% Regression
% no controls no random slope
X =  PredictorsForRegression_first_half;
varsNamesRegression = [ PredictorsNamesForRegression_first_half {'MMI'}];  
lme_noControls_noRandomSlope_MainTask_first_half = fitlm(X,MMI(First_half), ...
            'VarNames',varsNamesRegression);
adjRsqr_first_second(1,1) = lme_noControls_noRandomSlope_MainTask_first_half.Rsquared.Adjusted; 
clear X varsNamesRegression 

% no controls with random slope
X = [ones(length(PredictorsForRegression_first_half),1) PredictorsForRegression_first_half];
predictorsNamesRegression = [ {'Intercept'}  PredictorsNamesForRegression_first_half];  
Z = ones(length(PredictorsForRegression_first_half),1);
G = SID(First_half);
lme_noControls_withRandomSlope_MainTask_first_half = fitlmematrix(X,MMI(First_half),Z,G);
adjRsqr_first_second(2,1) = lme_noControls_withRandomSlope_MainTask_first_half.Rsquared.Adjusted;  
clear X Z G predictorsNamesRegression        

% with controls with random slope
X = [ones(length(PredictorsForRegression_first_half),1) PredictorsForRegression_first_half controls(First_half,:)];
predictorsNamesRegression = [ {'Intercept'}  PredictorsNamesForRegression_first_half {'SV'} {'RT'} {'slope'} {'IntersectY'} {'ChoiceDifficulty'} {'main_task_trial_num'} ];  
Z = ones(length(PredictorsForRegression_first_half),1);
G = SID(First_half);
lme_withControls_withRandomSlope_MainTask_first_half = fitlmematrix(X,MMI(First_half),Z,G);
adjRsqr_first_second(3,1) = lme_withControls_withRandomSlope_MainTask_first_half.Rsquared.Adjusted;  
clear X Z G predictorsNamesRegression       

% second half of the task
Second_half = find(main_task_trial_num>75);

% Elastic net - second half
[PredictorsForRegression_second_half, PredictorsNamesForRegression_second_half, ~, ~, ~] = elastic_net(mouse_features(Second_half,:),MMI(Second_half,:),titles_for_elastic_net);

% Regression
% no controls no random slope
X =  PredictorsForRegression_second_half;
varsNamesRegression = [ PredictorsNamesForRegression_second_half {'MMI'}];  
lme_noControls_noRandomSlope_MainTask_second_half = fitlm(X,MMI(Second_half), ...
            'VarNames',varsNamesRegression);
adjRsqr_first_second(1,2) = lme_noControls_noRandomSlope_MainTask_second_half.Rsquared.Adjusted; 
clear X varsNamesRegression 

% no controls with random slope
X = [ones(length(PredictorsForRegression_second_half),1) PredictorsForRegression_second_half];
predictorsNamesRegression = [ {'Intercept'}  PredictorsNamesForRegression_second_half];  
Z = ones(length(PredictorsForRegression_second_half),1);
G = SID(Second_half);
lme_noControls_withRandomSlope_MainTask_second_half = fitlmematrix(X,MMI(Second_half),Z,G);
adjRsqr_first_second(2,2) = lme_noControls_withRandomSlope_MainTask_second_half.Rsquared.Adjusted;  
clear X Z G predictorsNamesRegression        

% with controls with random slope
X = [ones(length(PredictorsForRegression_second_half),1) PredictorsForRegression_second_half controls(Second_half,:)];
predictorsNamesRegression = [ {'Intercept'}  PredictorsNamesForRegression_second_half {'SV'} {'RT'} {'slope'} {'IntersectY'} {'ChoiceDifficulty'} {'main_task_trial_num'} ];  
Z = ones(length(PredictorsForRegression_second_half),1);
G = SID(Second_half);
lme_withControls_withRandomSlope_MainTask_second_half = fitlmematrix(X,MMI(Second_half),Z,G);
adjRsqr_first_second(3,2) = lme_withControls_withRandomSlope_MainTask_second_half.Rsquared.Adjusted;  
clear X Z G predictorsNamesRegression  

% figure
adj_Rsqr_first_second_table = array2table(adjRsqr_first_second,'VariableNames',...
            {'first half' 'second half'} , ...
           'RowNames', {'no random slope, no controls', 'with random slope, no controls','with random slope, with controls'}); 

adjRsqr_forFigure(1,1:2) = adjRsqr_first_second(1,1:2);
adjRsqr_forFigure(2,1:2) = adjRsqr_first_second(2,1:2)-adjRsqr_first_second(1,1:2);
adjRsqr_forFigure(3,1:2) = adjRsqr_first_second(3,1:2)-adjRsqr_first_second(2,1:2);
figure;
b11=bar(1:2,adjRsqr_forFigure','stacked','LineWidth',1,'Facecolor','flat');
b11(1).FaceColor = [255 197 197]./255;
b11(2).FaceColor = 'none';
hatchfill2(b11(2),'single','HatchAngle',45,'HatchDensity',60,'hatchcolor',[255 197 197]./255,'HatchLineWidth',1);
b11(3).FaceColor = 'none';
hatchfill2(b11(3),'single','HatchAngle',0,'HatchDensity',60,'hatchcolor',[255 197 197]./255,'HatchLineWidth',1);
xticklabels({'first half', 'second half'}); ylabel('adj-Rsqr'); xlim([0.25 2.75]);

% out-of-sample prediction
% I use the betas from the regression from the first half of the main task
% on mouse features from the second half of the main task

locations_regressors = ismember(titles_for_elastic_net,PredictorsNamesForRegression_first_half);
index_regressors = find(locations_regressors==1);

X_new_second_half = [ones(length(mouse_features(Second_half,1)),1) mouse_features(Second_half,locations_regressors) controls(Second_half,:)];
Z_new = ones(length(mouse_features(Second_half,1)),1);
G_new = SID(Second_half);
%  prediction
yfit_second_first = predict(lme_withControls_withRandomSlope_MainTask_first_half,X_new_second_half, Z_new, G_new);

% correlations between fitted data and MMI
[corr_second_first, p_second_first, RL_second_first, RU_second_first] = corrcoef(MMI(Second_half),yfit_second_first,'Rows','complete'); 

% output for paper - cross-task prediction - scatterplot - predicted and actual MMI
figure; 
ax1=gca;
sizeMarker = 10;
scatter(yfit_second_first,MMI(Second_half),sizeMarker,'MarkerEdgeColor',[255 197 197]./255,'MarkerFaceColor',[255 197 197]./255 ,'LineWidth',1);
h = lsline(ax1); h.Color = 'k'; h.LineWidth = 0.5; h.LineStyle = '-'; xlim([min(yfit_second_first) max(yfit_second_first)]); ylim([0 1]);
xlabel('predicted inconsistency score based on mouse features (a.u.)'); ylabel('MMI');
        
