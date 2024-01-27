clear
clc

close all

%%
load('AllResults-Main task.mat');

% Drop subjects 109 and 309, whose trajectories were faulty 
[mal_subjects,j] = find(cellfun(@isempty,AllResults));
AllResults(mal_subjects,:,:) = [];
Subjects = length(AllResults);

%%
z = 0;

xflips_allSubs = [];
meanVel_allSubs = [];
AUC_allSubs = [];
numFixations_allSubs = [];
subjectsNum = [];

for s=1:Subjects
   fprintf('Subject %d\n',s);
    %load subject-specific mouse trajectory mat
    TrajectoryResults=AllResults{s,2};
        
     %load subject-specific behavioral choice data mat
    ChoiceData=AllResults{s,3};

    %identify trials with no choice 
    NAtrials=isnan(cat(1,ChoiceData.Real_Y));

    % Representative features (according to Fig. 3 in the manuscript)
    % flips features
    xflips          = cat(1,TrajectoryResults.Xflips); %number of flips on the x-axis
    meanVel         = cat(1,TrajectoryResults.meanVel); %mean mouse velocity (NIS/sec) during the trial
    AUC             = cat(1,TrajectoryResults.AUC); % the area between the trajectory and the "choice line".
    numFixations    = cat(1,TrajectoryResults.numFixations); % num of mouse fixations (>0.25 seconds in the same spatial bin)

    % subjects number
    SID             =   cat(1,ChoiceData.SID); %subject number

    % contacanate 
    xflips_allSubs = [xflips_allSubs;xflips];
    meanVel_allSubs = [meanVel_allSubs;meanVel];
    AUC_allSubs = [AUC_allSubs;AUC];
    numFixations_allSubs = [numFixations_allSubs;numFixations];
    subjectsNum = [subjectsNum; SID];

    clear xflips meanVel AUC numFixations SID

end

%% color palette
color_map=[223,87,79; %red
          242,213,84; %yellow
          61,199,127; %green
            82,160,210; %light blue
          86,102,191; %dark blue
          255, 255,255]; % white
color_map = color_map./255;      

%% Heat map of subjects' distributions features values across trials
% meanVel
% create bins of meanVel values
[meanVel_bins, edges_meanVel] = discretize(meanVel_allSubs,10);

% create a matrix by subjects and trials and meanVel bins
for subject=1:Subjects
 sid = cell2mat(AllResults(subject,1));
 subject_rows = find(subjectsNum==sid);
    for t=1:length(subject_rows)
        subjects_across_trials_meanVel(subject,t) = meanVel_bins(subject_rows(t),1);
    end
  clear sid subject_rows
end

subjects_across_trials_meanVel = sort(subjects_across_trials_meanVel,2,'ascend');
[subjects_across_trials_meanVel, index] = sortrows(subjects_across_trials_meanVel,'descend');

nan_values = find(isnan(subjects_across_trials_meanVel));
subjects_across_trials_meanVel(nan_values)=11;

%
figure()
im = imagesc(subjects_across_trials_meanVel);
xlabel('trials');
ylabel('subjects');
ytick=1:1:Subjects;
subjects_id = unique(subjectsNum);
set(gca,'Ytick',ytick,'YTickLabel',subjects_id(index),'FontSize',10);

n = 10;               %// number of colors
R = linspace(243,255,n)./255;  %// Red from 243 to 255
G = linspace(39,255,n)./255;   %// Green from 39 to 255
B = linspace(72,255,n)./255;  %// Blue from 72 to 255

colormap_meanVel = ( [R(:), G(:), B(:)] );  %// create colormap
black = [0,0,0];
colormap_meanVel = [colormap_meanVel];

colormap(colormap_meanVel); colorbar; caxis([1 10]);
title('meanVel');

clear subjects_across_trials

%% xflips
% create bins of xflips values
[xflips_bins, edges_xflips] = discretize(xflips_allSubs,10);

% create a matrix by subjects and trials and meanVel bins
for subject=1:Subjects
 sid = cell2mat(AllResults(subject,1));
 subject_rows = find(subjectsNum==sid);
    for t=1:length(subject_rows)
        subjects_across_trials_xflips(subject,t) = xflips_bins(subject_rows(t),1);
    end
  clear sid subject_rows
end

subjects_across_trials_xflips = sort(subjects_across_trials_xflips,2,'ascend');
[subjects_across_trials_xflips, index] = sortrows(subjects_across_trials_xflips,'descend');

nan_values = find(isnan(subjects_across_trials_xflips));
subjects_across_trials_xflips(nan_values)=11;

% 
figure()
im = imagesc(subjects_across_trials_xflips);
xlabel('trials');
ylabel('subjects');
ytick=1:1:Subjects;
subjects_id = unique(subjectsNum);
set(gca,'Ytick',ytick,'YTickLabel',subjects_id(index),'FontSize',10);

n = 10;               %// number of colors
R = linspace(0,255,n)./255;  %// Red from 0 to 255
G = linspace(149,255,n)./255;   %// Green from 149 to 255
B = linspace(47,255,n)./255;  %// Blue from 47 to 255

colormap_xflips = ( [R(:), G(:), B(:)] );  %// create colormap
black = [0,0,0];
colormap_xflips = [colormap_xflips];

colormap(colormap_xflips); colorbar; caxis([1 10]);
title('xflips');

%% AUC
% create bins of AUC values
edges_AUC = [0 0.025 0.05 0.075 0.1 0.125 0.15 0.175 0.2 0.225 1];
AUC_bins = discretize(AUC_allSubs,edges_AUC);

% create a matrix by subjects and trials and meanVel bins
for subject=1:Subjects
 sid = cell2mat(AllResults(subject,1));
 subject_rows = find(subjectsNum==sid);
    for t=1:length(subject_rows)
        subjects_across_trials_AUC(subject,t) = AUC_bins(subject_rows(t),1);
    end
  clear sid subject_rows
end

subjects_across_trials_AUC = sort(subjects_across_trials_AUC,2,'ascend');
[subjects_across_trials_AUC, index] = sortrows(subjects_across_trials_AUC,'descend');

nan_values = find(isnan(subjects_across_trials_AUC));
subjects_across_trials_AUC(nan_values)=11;

% 
figure()
im = imagesc(subjects_across_trials_AUC);
xlabel('trials');
ylabel('subjects');
ytick=1:1:Subjects;
subjects_id = unique(subjectsNum);
set(gca,'Ytick',ytick,'YTickLabel',subjects_id(index),'FontSize',10);

n = 10;               %// number of colors
R = linspace(0,255,n)./255;  %// Red from 0 to 255
G = linspace(71,255,n)./255;   %// Green from 149 to 255
B = linspace(189,255,n)./255;  %// Blue from 67 to 255

colormap_AUC = ( [R(:), G(:), B(:)] );  %// create colormap
black = [0,0,0];
colormap_AUC = [colormap_AUC];

colormap(colormap_AUC); colorbar; caxis([1 10]);
title('AUC');

%% numFixations
% create bins of numFixations values
[numFixations_bins, edges_numFixations] = discretize(numFixations_allSubs,10);

% create a matrix by subjects and trials and meanVel bins
for subject=1:Subjects
 sid = cell2mat(AllResults(subject,1));
 subject_rows = find(subjectsNum==sid);
    for t=1:length(subject_rows)
        subjects_across_trials_numFixations(subject,t) = numFixations_bins(subject_rows(t),1);
    end
  clear sid subject_rows
end

subjects_across_trials_numFixations = sort(subjects_across_trials_numFixations,2,'ascend');
[subjects_across_trials_numFixations, index] = sortrows(subjects_across_trials_numFixations,'descend');

nan_values = find(isnan(subjects_across_trials_numFixations));
subjects_across_trials_numFixations(nan_values)=11;

% 
figure()
im = imagesc(subjects_across_trials_numFixations);
xlabel('trials');
ylabel('subjects');
ytick=1:1:Subjects;
subjects_id = unique(subjectsNum);
set(gca,'Ytick',ytick,'YTickLabel',subjects_id(index),'FontSize',10);

n = 10;               %// number of colors
R = linspace(130,255,n)./255;  %// Red from 130 to 255
G = linspace(0,255,n)./255;   %// Green from 0 to 255
B = linspace(172,255,n)./255;  %// Blue from 172 to 255

colormap_numFixations = ( [R(:), G(:), B(:)] );  %// create colormap
black = [0,0,0];
colormap_numFixations = [colormap_numFixations];

colormap(colormap_numFixations); colorbar; caxis([1 10]);
title('numFixations');

