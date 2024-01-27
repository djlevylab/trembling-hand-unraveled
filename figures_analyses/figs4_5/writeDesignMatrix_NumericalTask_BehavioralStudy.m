clear
clc


load('AllResults-Numerical_localizer.mat');

% Drop subjects 109 and 309, whose trajectories were faulty 
[mal_subjects,j] = find(cellfun(@isempty,AllResults));
AllResults(mal_subjects,:,:) = [];
Subjects = length(AllResults);

DesignMatrixNumerical = [];

%%

z = 0;
for s=1:Subjects
   fprintf('Subject %d\n',s);
    %load subject-specific mouse trajectory mat
    TrajectoryResults=AllResults{s,2};
        
     %load subject-specific behavioral choice data mat
    ChoiceData=AllResults{s,3};

    %identify trials with no choice 
    NAtrials=isnan(cat(1,ChoiceData.realY));

    % Trajectory features
    %flips features
    xflips          = cat(1,TrajectoryResults.Xflips); %number of flips on the x-axis
    yflips          = cat(1,TrajectoryResults.Yflips); %number of flips on y-axis
    MaxXflip        = cat(1,TrajectoryResults.MaxXflip); %maximal flip on the X-axis
    MaxYflip        = cat(1,TrajectoryResults.MaxYflip); %maximal flip on the Y-axis

    %velocity features
    MaxProgress     = cat(1,TrajectoryResults.MaxProgress); %maximal progress made within a single sample (beyond all other samples)
    medianProgress  = cat(1,TrajectoryResults.MedianProgress); %the median progress (in one sample) made in each tiral
    Time2MaxProgress= cat(1,TrajectoryResults.Time2MaxProgress); %time in sec it took the subject to perform the largest progess. Perhaps can indicate on the moment of epiphany.
    meanVel         = cat(1,TrajectoryResults.meanVel); %mean mouse velocity (NIS/sec) during the trial
    MaxVel          = cat(1,TrajectoryResults.MaxVel); %maximal mouse velocity (NIS/sec) during the trial
    medianVel       = cat(1,TrajectoryResults.MedianVel); %median mouse velocity (NIS/sec) during the trial
    VelOver10       = cat(1,TrajectoryResults.VelOver10); %number of samples in which the velocity was higher than 10 NIS/sec
    EndofTrialVel   = cat(1,TrajectoryResults.EndofTrialVel); %velocity at the end of the trial
    MaxVelTime      = cat(1,TrajectoryResults.MaxVelTime); %timing of maximal velocity
    meanAcc         = cat(1,TrajectoryResults.meanAcc); %mean mouse accelaration (NIS/sec^2) during the trial
    MaxAcc          = cat(1,TrajectoryResults.MaxAcc); %maximal mouse accelaration (NIS/sec^2) during the trial
    N_Accs          = cat(1,TrajectoryResults.N_Accs); %number of acceleration eriods within a trial

    %curvature features
    maxArcLen       = cat(1,TrajectoryResults.ArcLenMax); %length of the biggest arc in the trajectory
    Curveness       = cat(1,TrajectoryResults.Curveness); %sum of the normals of all the curvatures
    Angle_STD       = cat(1,TrajectoryResults.AngleSTD); %std of the average angle (tan) between the trajectory and the choice line
    MD              = cat(1,TrajectoryResults.MD); %maximal distance of mouse from "choice line" (straight line from origin)
    AUC             = cat(1,TrajectoryResults.AUC); % the area between the trajectory and the "choice line".

    %fixations features
    Layover         = cat(1,TrajectoryResults.Layover); % layover near the origin (measured in sec)
    numFixations    = cat(1,TrajectoryResults.numFixations); % num of mouse fixations (>0.25 seconds in the same spatial bin)
    MaxFixTime      = cat(1,TrajectoryResults.MaxFixTime); %the duration of maximal mouse fixation(in sec)
    MaxFixChoiceDist= cat(1,TrajectoryResults.MaxFixChoiceDist); %distance to choice coords from longest fixation coords
    numBudgetFixations = cat(1,TrajectoryResults.numBudgetFixations); %number of fixations near/on budget line (up to 2 NIS distance)
    AveFixTime      = cat(1,TrajectoryResults.AveFixTime); %average fixations durations 

    %deviation features
    AboveLine       = cat(1,TrajectoryResults.AboveLine); %time spent beyond the budget line (in sec)
    XTimeOutofBounds= cat(1,TrajectoryResults.XTimeOutofBounds); %time spent outside the grid - X_axis (in sec)
    YTimeOutofBounds= cat(1,TrajectoryResults.YTimeOutofBounds); %time spent outside the grid - Y_axis (in sec)

    %economic-theory related features
    TimeNearPredicted= cat(1,TrajectoryResults.TimeNearPredicted); %time spent near the predicted bundle (up to 2 NIS above and below)
    TimeBudgetLine  = cat(1,TrajectoryResults.TimeBudgetLine); %time spent near the budget line (up to 2 NIS above and below)

    %entropy features
    XSampEn4        = cat(1,TrajectoryResults.XSampEn4); %sample entropy for a given window size (4) and axis (X)
    YSampEn4        = cat(1,TrajectoryResults.YSampEn4); %sample entropy for a given window size (4) and axis (Y)

   
    subject_stata_matrix = [xflips yflips MaxXflip MaxYflip ...
                           MaxProgress medianProgress Time2MaxProgress meanVel MaxVel medianVel VelOver10 ...
                           EndofTrialVel MaxVelTime meanAcc MaxAcc N_Accs ...
                           maxArcLen Curveness ...
                           Angle_STD MD AUC ...
                           Layover numFixations MaxFixTime MaxFixChoiceDist ...
                           numBudgetFixations AveFixTime ...
                           AboveLine XTimeOutofBounds YTimeOutofBounds TimeNearPredicted TimeBudgetLine ...
                           XSampEn4 YSampEn4];
     
    % other parameters
    %Budget set and trial parameters
    a                           =   abs(cat(1,ChoiceData.a)); %budget set slope
    b                           =   cat(1,ChoiceData.b); %Y-axis intersection
    RT                          =   cat(1,ChoiceData.RT); %RT (sec)
    obs_num                     =   cat(1,ChoiceData.Trial); %trial number (1 to 150)
    block_num                   =   cat(1,ChoiceData.Block); %block number (1 to 3)
    X_allocation                =   cat(1,ChoiceData.lineX); %button-press in X-axis
    Y_allocation                =   cat(1,ChoiceData.lineY); %Choice in Y-axis
    X_target                    =   cat(1,ChoiceData.targetX); % X-coordinate in he target (choice in the main task)
    Y_target                    =   cat(1,ChoiceData.targetY); % Y-coordinate in he target (choice in the main task)
    subNum                      =   cell2mat(AllResults(s,1));
    SID                         =   subNum.*ones(length(Y_target),1);

    DesignMatrixNumerical = [DesignMatrixNumerical; SID obs_num block_num a b X_allocation Y_allocation ...
                        X_target Y_target RT subject_stata_matrix]; 

    clear SID obs_num block_num a b RT X_allocation Y_allocation 
    clear X_target Y_target MainTaskTrialNumber
    clear xflips yflips MaxXflip MaxYflip 
    clear MaxProgress medianProgress Time2MaxProgress meanVel MaxVel medianVel VelOver10 
    clear EndofTrialVel MaxVelTime meanAcc MaxAcc medianAcc N_Accs 
    clear maxArcLen RadiusofCurvMax RadiusofCurvMean RadiusofCurvSTD Curveness 
    clear Angle_STD MD AUC 
    clear Layover numFixations MaxFixTime MaxFixChoiceDist  
    clear numBudgetFixations AveFixTime 
    clear AboveLine XTimeOutofBounds YTimeOutofBounds TimeInFOSD TimeNearPredicted TimeBudgetLine 
    clear XSampEn4 YSampEn4 subject_stata_matrix_zscored subject_stata_matrix

end

 num_col = size(DesignMatrixNumerical,2);


%% Normalize - z-score mouse features
for col=12:num_col
    DesignMatrixNumerical(:,col)=(DesignMatrixNumerical(:,col)-nanmean(DesignMatrixNumerical(:,col)))./nanstd(DesignMatrixNumerical(:,col));
end
disp('Done Normalizing!');

%% Clean
num_to_round = (length(find(DesignMatrixNumerical(:,12:end)>5 | DesignMatrixNumerical(:,12:end)<-5))./numel(DesignMatrixNumerical(:,12:end)))*100;
disp(['Data rounded: ' num2str(num_to_round) '%']);
for col=12:num_col
%     [B,TF] = rmoutliers(data2(:,col),'quartiles','ThresholdFactor',4);
%     data2(TF,col)=NaN;
    DesignMatrixNumerical((DesignMatrixNumerical(:,col))>5,col)=5;
    DesignMatrixNumerical((DesignMatrixNumerical(:,col))<-5,col)=-5;

end
disp('Done Cleaning!');

% NaN Reports
samp_nans=[];
for col=12:num_col
    allnans(col)=mean(isnan(DesignMatrixNumerical(:,col)));
    samp_nans=[samp_nans isnan(DesignMatrixNumerical(:,col))];
end


%% Output
DesignMatrix_table = array2table(DesignMatrixNumerical,'VariableNames',...
            {'SID' 'obs_num' 'block_num' 'a' 'b' 'X_allocation' 'Y_allocation' ...
            'X_target' 'Y_target' 'RT' ...           
            'xflips' 'yflips' 'MaxXflip' 'MaxYflip' ...
            'MaxProgress' 'medianProgress' 'Time2MaxProgress' 'meanVel' 'MaxVel' 'medianVel' 'VelOver10' ...
            'EndofTrialVel' 'MaxVelTime' 'meanAcc' 'MaxAcc' 'N_Accs' ...
            'maxArcLen' 'Curveness' ...
            'Angle_STD' 'MD' 'AUC' ...
            'Layover' 'numFixations' 'MaxFixTime' 'MaxFixChoiceDist' ...
            'numBudgetFixations' 'AveFixTime' ...
            'AboveLine' 'XTimeOutofBounds' 'YTimeOutofBounds' 'TimeNearPredicted' 'TimeBudgetLine' ...
            'XSampEn4' 'YSampEn4'}); 
    
file_name = ['DesignMatrix_ZScoredAcrossSubject_NumericalTask_BehavioralStudy.xls'];
writetable(DesignMatrix_table, file_name, 'WriteVariableNames', true);
