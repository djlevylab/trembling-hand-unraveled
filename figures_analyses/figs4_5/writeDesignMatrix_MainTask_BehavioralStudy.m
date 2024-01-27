clear
clc

load('AllResults-Main task.mat');

% Drop subjects 109 and 309, whose trajectories were faulty 
[mal_subjects,j] = find(cellfun(@isempty,AllResults));
AllResults(mal_subjects,:,:) = [];
Subjects = length(AllResults);

DesignMatrix = [];

%%

z = 0;
for s=1:Subjects
   fprintf('Subject %d\n',s);
    %load subject-specific mouse trajectory mat
    TrajectoryResults=AllResults{s,2};
        
     %load subject-specific behavioral choice data mat
    ChoiceData=AllResults{s,3};

    %identify trials with no choice 
    NAtrials=isnan(cat(1,ChoiceData.Real_Y));

    % Trajectory features
    % flips features
    xflips          = cat(1,TrajectoryResults.Xflips); %number of flips on the x-axis
    yflips          = cat(1,TrajectoryResults.Yflips); %number of flips on y-axis
    MaxXflip        = cat(1,TrajectoryResults.MaxXflip); %maximal flip on the X-axis
    MaxYflip        = cat(1,TrajectoryResults.MaxYflip); %maximal flip on the Y-axis

    % velocity features
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

    %economic-theory related features                                                             % and up to 2 NIS around the y coordinate of the budget line. 
    TimeNearPredicted= cat(1,TrajectoryResults.TimeNearPredicted); %time spent near the predicted bundle (up to 2 NIS above and below)
    TimeBudgetLine  = cat(1,TrajectoryResults.TimeBudgetLine); %time spent near the budget line (up to 2 NIS above and below)

    %entropy features
    XSampEn4        = cat(1,TrajectoryResults.XSampEn4); %sample entropy for a given window size (4) and axis (X)
    YSampEn4        = cat(1,TrajectoryResults.YSampEn4); %sample entropy for a given window size (4) and axis (Y)

    % Behavioral measurements (inconsistency, SV and choice difficulty)
    % Inconsistency indices - MMI
    MMI_in_Sample_component     =   cat(1,ChoiceData.MMI_In_Sample_component); 
    MMI_in_Sample_diff          =   cat(1,ChoiceData.MMI_In_Sample_diff); 
    MMI_Out_of_Sample_alt       =   cat(1,ChoiceData.MMI_Out_of_Sample_alt); 
    MMI_Out_of_Sample_diff      =   cat(1,ChoiceData.MMI_Out_of_Sample_diff); 

    %Inconsistency indices - Varian
    Varian_in_Sample_component  =   cat(1,ChoiceData.Varian_In_Sample_component); 
    Varian_in_Sample_diff       =   cat(1,ChoiceData.Varian_In_Sample_diff); 
    Varian_Out_of_Sample_alt    =   cat(1,ChoiceData.Varian_Out_of_Sample_alt); 
    Varian_Out_of_Sample_diff   =   cat(1,ChoiceData.Varian_Out_of_Sample_diff);

    %Other behavioral measurements
    SV                          =   (cat(1,ChoiceData.SV)); %subjective value (arbitrary units)
    ChoiceDifficulty            =   (cat(1,ChoiceData.ChoiceDifficulty)); %choice difficulty (calculated using SV with arbitrary units)

    subject_state_matrix = [SV ChoiceDifficulty ...
                           xflips yflips MaxXflip MaxYflip ...
                           MaxProgress medianProgress Time2MaxProgress meanVel MaxVel medianVel VelOver10 ...
                           EndofTrialVel MaxVelTime meanAcc MaxAcc N_Accs ...
                           maxArcLen Curveness ...
                           Angle_STD MD AUC ...
                           Layover numFixations MaxFixTime MaxFixChoiceDist ...
                           numBudgetFixations AveFixTime ...
                           AboveLine XTimeOutofBounds YTimeOutofBounds TimeNearPredicted TimeBudgetLine ...
                           XSampEn4 YSampEn4];
    
    % z_score (only SV and choice difficulty) - ordinal measures 
    stata_matrix_zscored = subject_state_matrix;
     for col=1:2
         stata_matrix_zscored(:,col)=(subject_state_matrix(:,col)-nanmean(subject_state_matrix(:,col))) / nanstd(subject_state_matrix(:,col));
     end

    % other parameters
    %Budget set and trial parameters
    a                           =   abs(cat(1,ChoiceData.slope)); %budget set slope
    b                           =   cat(1,ChoiceData.IntersectionY); %Y-axis intersection
    RT                          =   cat(1,ChoiceData.RT); %RT (sec)
    obs_num                     =   cat(1,ChoiceData.obs); %trial number (1 to 150)
    block_num                   =   cat(1,ChoiceData.block); %block number (1 to 3)
    SID                         =   cat(1,ChoiceData.SID); %subject number
    X_allocation                =   cat(1,ChoiceData.X_allocation); %Choice in X-axis
    Y_allocation                =   cat(1,ChoiceData.Y_allocation); %Choice in Y-axis
    Predicted_Bundle_X          =   cat(1,ChoiceData.Predicted_Bundle_X); %Predicted bundle (X-axis) according to elicited parameters
    Predicted_Bundle_Y          =   cat(1,ChoiceData.Predicted_Bundle_Y); %Predicted bundle (Y-axis) according to elicited parameters
%     Endowment                   =   cat(1,ChoiceData.Endowment); % Endowment measured by the safe bundle

    DesignMatrix = [DesignMatrix; SID obs_num block_num a b X_allocation Y_allocation ...
                   Predicted_Bundle_X Predicted_Bundle_Y RT MMI_in_Sample_component MMI_in_Sample_diff ...
                   MMI_Out_of_Sample_alt MMI_Out_of_Sample_diff ...
                   Varian_in_Sample_component Varian_in_Sample_diff Varian_Out_of_Sample_alt ...
                   Varian_Out_of_Sample_diff stata_matrix_zscored]; 

    clear SID obs_num block_num a b RT X_allocation Y_allocation 
    clear MMI_in_Sample_component MMI_in_Sample_diff MMI_Out_of_Sample_alt MMI_Out_of_Sample_diff 
    clear Varian_in_Sample_component Varian_Out_of_Sample_alt Varian_in_Sample_diff Varian_Out_of_Sample_diff 
    clear SV ChoiceDifficulty Predicted_Bundle_X Predicted_Bundle_Y 
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

num_col = size(DesignMatrix,2);

%% Normalize
 for col=21:num_col % z_score mouse features
  DesignMatrix(:,col)=(DesignMatrix(:,col)-nanmean(DesignMatrix(:,col))) / nanstd(DesignMatrix(:,col));
end

%% Clean
num_to_round = (length(find(DesignMatrix(:,21:end)>5 | DesignMatrix(:,21:end)<-5))./numel(DesignMatrix(:,21:end)))*100;
disp(['Data rounded: ' num2str(num_to_round) '%']);
for col=21:num_col % just mouse features
    DesignMatrix((DesignMatrix(:,col))>5,col)=5;
    DesignMatrix((DesignMatrix(:,col))<-5,col)=-5;

end
disp('Done Cleaning!');

% NaN Reports
samp_nans=[];
for col=21:num_col % just mouse features
    allnans(col)=mean(isnan(DesignMatrix(:,col)));
    samp_nans=[samp_nans isnan(DesignMatrix(:,col))];
end
 disp(allnans)


%% Output
DesignMatrix_table = array2table(DesignMatrix,'VariableNames',...
            {'SID' 'obs_num' 'block_num' 'a' 'b' 'X_allocation' 'Y_allocation' ...
            'Predicted_Bundle_X' 'Predicted_Bundle_Y' 'RT' ...
            'MMI_in_Sample_component' 'MMI_in_Sample_diff' 'MMI_Out_of_Sample_alt' 'MMI_Out_of_Sample_diff' ...
            'Varian_in_Sample_component' 'Varian_in_Sample_diff' 'Varian_Out_of_Sample_alt' 'Varian_Out_of_Sample_diff' ...
            'SV' 'ChoiceDifficulty' ...
            'xflips' 'yflips' 'MaxXflip' 'MaxYflip' ...
            'MaxProgress' 'medianProgress' 'Time2MaxProgress' 'meanVel' 'MaxVel' 'medianVel' 'VelOver10' ...
            'EndofTrialVel' 'MaxVelTime' 'meanAcc' 'MaxAcc' 'N_Accs' ...
            'maxArcLen' 'Curveness' ...
            'Angle_STD' 'MD' 'AUC' ...
            'Layover' 'numFixations' 'MaxFixTime' 'MaxFixChoiceDist' ...
            'numBudgetFixations' 'AveFixTime' ...
            'AboveLine' 'XTimeOutofBounds' 'YTimeOutofBounds' 'TimeNearPredicted' 'TimeBudgetLine' ...
            'XSampEn4' 'YSampEn4'}); 
    
file_name = ['DesignMatrix_ZScoredAcrossSubjectWithNans_MainTaskBehavioralStudy.xls'];
writetable(DesignMatrix_table, file_name, 'WriteVariableNames', true);
