function Results=trajectory_analysis(trajectory,time,trialdata,params)

%% Clean Trajectory
[new_traj, new_time]=clean_traj(trajectory,time);
Results.RawTraj=trajectory;
Results.OrgTraj=new_traj;

%% Set Parameters
SegCount=size(new_traj,1);
try
    PredictedCoords=[trialdata{params.PredCols(1)} trialdata{params.PredCols(2)}];
catch
    PredictedCoords=[0 0];
end
a=trialdata{params.a_col}; bY=trialdata{params.b_col};
bX=-bY/a;

% %% Spatial Remapping: Rescale From Pixels to Standard Coord. Space
% %Convert all clicks to right-hand side - end is all (1,1)
% RescX=(Results.OrgTraj(:,1)-Results.OrgTraj(1,1))...
%     /(Results.OrgTraj(end,1)-Results.OrgTraj(1,1));
% RescY=(Results.OrgTraj(:,2)-Results.OrgTraj(1,2))...
%     /(Results.OrgTraj(end,2)-Results.OrgTraj(1,2));
% Results.RescaleTraj=[RescX RescY];

% %% Linear Temporal Interpolation on RescaleTraj (to set timepoints):
% TimeStamps = time;
% TimeSamples = linspace(TimeStamps(1),TimeStamps(end),100);
%     %Convert all clicks to right-hand side - end is all (1,1)
% Results.InterpRescale(:,1)=interp1(TimeStamps,... % timepoints (x)
%     Results.RescaleTraj(:,1),... % X coordinate (f(x))
%     TimeSamples,'linear');% timepoints to use for interploating - query points (xq)
% Results.InterpRescale(:,2)=interp1(TimeStamps,... % timepoints
%     Results.RescaleTraj(:,2),... % Y coordinate
%     TimeSamples,'linear'); % timepoints to use for interploating

%% Linear Temporal Interpolation for Original trajectory:
TimeStamps = new_time;
NewSegCount=100;
% TimeSamples = linspace(TimeStamps(1),TimeStamps(end),length(Results.OrgTraj(:,1)));
TimeSamples = linspace(TimeStamps(1),TimeStamps(end),NewSegCount);

    %Convert all clicks to right-hand side - end is all (1,1)
Results.InterpOrg(:,1)=interp1(TimeStamps,... % timepoints (x)
    Results.OrgTraj(:,1),... % X coordinate (f(x))
    TimeSamples,'linear');% timepoints to use for interploating - query points (xq)
Results.InterpOrg(:,2)=interp1(TimeStamps,... % timepoints
    Results.OrgTraj(:,2),... % Y coordinate
    TimeSamples,'linear'); % timepoints to use for interploating
    %Retain choice - end is (-1,1) or (1,1)
Results.interpTime=interp1(TimeStamps,... % timepoints
    TimeStamps,... % Y coordinate
    TimeSamples,'linear')'; 
Results.SR=mean(diff(Results.interpTime));

%% Linear Temporal Interpolation for Original trajectory - Multiple Samples:
TimeStamps = new_time;
% TimeSamples = linspace(TimeStamps(1),TimeStamps(end),length(Results.OrgTraj(:,1)));
TimeSamples = linspace(TimeStamps(1),TimeStamps(end),1000);

    %Convert all clicks to right-hand side - end is all (1,1)
Results.InterpOrgHighRes(:,1)=interp1(TimeStamps,... % timepoints (x)
    Results.OrgTraj(:,1),... % X coordinate (f(x))
    TimeSamples,'linear');% timepoints to use for interploating - query points (xq)
Results.InterpOrgHighRes(:,2)=interp1(TimeStamps,... % timepoints
    Results.OrgTraj(:,2),... % Y coordinate
    TimeSamples,'linear'); % timepoints to use for interploating
    %Retain choice - end is (-1,1) or (1,1)
Results.interpTimeHighRes=interp1(TimeStamps,... % timepoints
    TimeStamps,... % Y coordinate
    TimeSamples,'linear')'; 
Results.SRHighRes=mean(diff(Results.interpTimeHighRes));


NormX=(Results.InterpOrgHighRes(:,1)-Results.InterpOrgHighRes(1,1))...
    /(Results.InterpOrgHighRes(end,1)-Results.InterpOrgHighRes(1,1));
NormY=(Results.InterpOrgHighRes(:,2)-Results.InterpOrgHighRes(1,2))...
    /(Results.InterpOrgHighRes(end,2)-Results.InterpOrgHighRes(1,2));
Results.NormTraj=[NormX NormY];

%% Layover time at start point
FirstX=Results.InterpOrgHighRes(1,1);
FirstY=Results.InterpOrgHighRes(1,2);
LayoverThresh=2;
LengthCheck=1:round(0.3*length(Results.InterpOrgHighRes(:,1)));
Results.Layover=sum((Results.InterpOrgHighRes(LengthCheck,1)<=(FirstX+LayoverThresh)) & (Results.InterpOrgHighRes(LengthCheck,2)<=(FirstY+LayoverThresh)))*Results.SRHighRes;
% Results.Layover=sum(Results.RescaleTraj(:,1)<=0.01 & Results.RescaleTraj(:,2)<=0.01)/SegCount;

%% Mouse Trajectory Angle
Results.Angles= ...
    abs(atand(Results.InterpOrg(:,1)./ ...
    Results.InterpOrg(:,2))).* ...
    sign(Results.InterpOrg(:,1));
%         % Optional: exclude points where mouse moving downward
%     ExcludeInd = [false; Results.TemporalDirect(2:end,2)< ...
%          Results.TemporalDirect(1:(end-1),2)];
%     Results.Angle(ExcludeInd)=NaN;
Results.Angles(Results.InterpOrg(:,1)==0)=0; % Convert NaNs to 0s
Results.Angles=deg2rad(Results.Angles);
Results.AngleSTD=nanstd(Results.Angles);

%% Mouse Segment & Comulative Progress
X_Progress=[]; Y_Progress=[];
X_Progress(1)=Results.InterpOrg(1,1);
Y_Progress(1)=Results.InterpOrg(1,2);
X_Progress(2:NewSegCount)=diff(Results.InterpOrg(:,1));
Y_Progress(2:NewSegCount)=diff(Results.InterpOrg(:,2));
Results.Seg_D_Progress=sqrt(X_Progress'.^2+Y_Progress'.^2);
Comul=[]; Comul(1:NewSegCount)=0;
for k=2:NewSegCount
    Comul(k)=Comul(k-1)+Results.Seg_D_Progress(k);
end
Results.Comul_Progress=Comul';
Results.MaxProgress=max(Results.Seg_D_Progress);
Results.Time2MaxProgress=find(Results.Seg_D_Progress==Results.MaxProgress,1,'first')*Results.SR;
Results.MedianProgress=median(Results.Seg_D_Progress);
% Mouse Trial Velocity per Segment
% Results.TrVeloc=Results.D_Start./[1:100]';
% raw_distance=sqrt(Results.RescaleTraj(:,1).^2+Results.RescaleTraj(:,2).^2);
% raw_velocity=raw_distance./[1:SegCount]';
% Results.InterpTrVeloc(:,1)=interp1(TimeStamps,raw_velocity,TimeSamples,'linear');
%%%% choose interpolation of raw or use interpolated?
Results.Velocity=[0; Results.Seg_D_Progress(2:end)./diff(Results.interpTime)];
Results.VelOver10=sum(Results.Velocity>10)*Results.SR;
Results.meanVel=mean(Results.Velocity);
Results.MedianVel=median(Results.Velocity);
[M,I]=max(Results.Velocity);
Results.MaxVel=M;
Results.MaxVelTime=I*Results.SR;
Results.Acceleration=[0; diff(Results.Velocity)./diff(Results.interpTime)];
Results.meanAcc=mean(abs(Results.Acceleration));
% Results.medianAcc=median(abs(Results.Acceleration)); Nonsense
Results.MaxAcc=max(Results.Acceleration);
Results.EndofTrialVel=mean(Results.Velocity(end-10:end));
[pks,~] = findpeaks(Results.Acceleration,'MinPeakDistance',3,'MinPeakHeight',50,'Threshold',1);
Results.N_Accs=length(pks);


%% Curvature
%Find trajectory after layover and before reaching budget line
x = 0:(1/10):100;
y = a*x + bY;
[~, closest_D]=dsearchn([x' y'],Results.InterpOrg);
startind=find(Results.InterpOrg(:,1)>1 | Results.InterpOrg(:,2)>1,1,'first');
for e=2:20
    endind=find(closest_D(startind+1:end)<e,1,'first')+startind;
    if ~isempty(endind)
        break
    end
end
if isempty(endind)
    endind=find(closest_D(startind+1:end)<100,1,'first')+startind;
end

Results.Traj2Line=Results.InterpOrg(startind:endind,:);
TotDist=sum(Results.Seg_D_Progress(startind:endind));
%Get equal distance trajectory
Results.EqualTraj=interparc(max([round(TotDist/2),3]),Results.Traj2Line(:,1)+rand(1,size(Results.Traj2Line,1))'/100,Results.Traj2Line(:,2));
[L,R,K] = curvature(Results.EqualTraj);
Results.ArcLenCum=L;
[L,TF] = rmoutliers(L,'quartiles');
Results.ArcLenMax=max(L);
Results.RadiusofCurv=R;
% [R,TF] = rmoutliers(R,'quartiles');
% Results.RadiusofCurvMed=nanmedian(R);
% Results.RadiusofCurvIQR=iqr(R);
% Results.RadiusofCurvSTD=nanstd(R);
Results.Curvatures=K;
Results.Curveness=sum(vecnorm(K(~isnan(K(:,1)),:)',2));

% %plot
% nexttile;
% plot(Results.InterpOrg(:,1),Results.InterpOrg(:,2),'.');
% title('InterpOrg');
% xlim([-10 70]); ylim([-10 70]);
% nexttile
% h = plot(Results.EqualTraj(:,1),Results.EqualTraj(:,2)); grid on;
% set(h,'marker','.');
% xlabel x ; ylabel y
% title(['EqualTraj, Curveness: ' num2str(Results.Curveness)])
% hold on
% quiver(Results.EqualTraj(:,1)',Results.EqualTraj(:,2)',K(:,1)',K(:,2)');
% xlim([-10 70]); ylim([-10 70]);
% nexttile
% plot(L,R)
% xlabel L ; ylabel R

%% Get MD
%the maximal distance from a trajectory point to the straight line to choice
Results.MD=get_MD(Results.InterpOrg)/100;

%% Calculate AUC
%Area between trajectory and straight line to choice
% Results.AUC=abs(trapz(Results.InterpRescale(:,1),Results.InterpRescale(:,2))-0.5);
Results.AUC=traj2mask(Results.InterpOrgHighRes);

%% Fixation Points
tolerance=3;
fixtimeThresh=0.2;
Edges={0:5:100,0:5:100};
[n, c] = hist3([Results.InterpOrgHighRes(:,1),Results.InterpOrgHighRes(:,2)],'edges',Edges);
FixationTimes=n*Results.SRHighRes;
FixationTimes(1:3,1:3)=0; %Clear Layover
%Set Choice space to zero in grid
allpoints=combvec(c{1}, c{2})';
[K1,D1]=dsearchn(Results.OrgTraj(end,:),allpoints);
[M1,I1] = min(D1); closest_grid=allpoints(I1,:); 
cg_xcoord=find(closest_grid(1)==c{1},1,'first');
cg_ycoord=find(closest_grid(2)==c{2},1,'first'); 
%Get Maximal fixation time without choice space
FixationTimes(cg_xcoord,cg_ycoord)=0;

Results.numFixations=sum(FixationTimes>fixtimeThresh,'all');
if Results.numFixations>0
    Results.AveFixTime=mean(FixationTimes(FixationTimes>fixtimeThresh));
    Results.MaxFixTime=max(FixationTimes(FixationTimes>fixtimeThresh),[],'all');
    [x,y]=find(FixationTimes==Results.MaxFixTime,1,'first');
    Results.MaxFixChoiceDist = pdist([c{1}(x) c{2}(y) ; Results.OrgTraj(end,:)],'euclidean');
    Results.MaxFixPredictDist = pdist([c{1}(x) c{2}(y) ; PredictedCoords],'euclidean');
    %Calculate distance from fixations coords to budget line
    [x,y]=find(FixationTimes>fixtimeThresh); FixCoords=[c{1}(x)' c{2}(y)'];
    for i=1:size(FixCoords,1)
        Q1=[bX 0]; Q2=[0 bY];
        P=FixCoords(i,:);
        ds(i)=abs(det([Q2-Q1,;P-Q1]))/norm(Q2-Q1);
    end
    Results.numBudgetFixations=sum(ds<tolerance);
else
    Results.AveFixTime=0;
    Results.MaxFixTime=0;
    Results.MaxFixChoiceDist=0;
    Results.MaxFixPredictDist=0;
    Results.numBudgetFixations=0;
end

%% Time Above Budget Line
Results.AboveLine=sum(Results.InterpOrg(:,2) > a*Results.InterpOrg(:,1)+bY)*Results.SR;

%% Time Out of Bounds
%How much time did trajectory appear below 0 and above 100
Results.XTimeOutofBounds=sum(Results.InterpOrg(:,1)>100 | Results.InterpOrg(:,1)<0)*Results.SR;
Results.YTimeOutofBounds=sum(Results.InterpOrg(:,2)>100 | Results.InterpOrg(:,2)<0)*Results.SR;

%% Time in dominated budget sets
%Calculate point of budget equality
syms x45 y45;
eqn1=a*x45+bY - y45 == 0;
eqn2=x45 -y45 == 0;
[A,B] = equationsToMatrix([eqn1, eqn2], [x45, y45]);
Solution = double(linsolve(A,B));

%Find samples that were around dominated budget sets
tolerance=4; %4?
if a>=-1
    RelevantCoords=Results.InterpOrg(Results.InterpOrg(:,1)<(Solution(1)-tolerance),:);
    FinalCoords=RelevantCoords(abs(RelevantCoords(:,2) - (a*RelevantCoords(:,1)+bY))<tolerance,:);
elseif a<-1
    RelevantCoords=Results.InterpOrg(Results.InterpOrg(:,1)>(Solution(1)+tolerance),:);
    FinalCoords=RelevantCoords(abs(RelevantCoords(:,2) - (a*RelevantCoords(:,1)+bY))<tolerance,:);  
end
Results.TimeInFOSD=size(FinalCoords,1)*Results.SR;

%% Time Near Predicted Budget Set
NearPredThresh=10;
Results.D_Predicted=sqrt((PredictedCoords(1,1)-Results.InterpOrgHighRes(:,1)).^2+(PredictedCoords(1,2)-Results.InterpOrgHighRes(:,2)).^2);
Results.TimeNearPredicted=sum((Results.D_Predicted<NearPredThresh))*Results.SRHighRes;

%% Time Around Budget Line
tolerance=2;
Results.TimeBudgetLine=sum(abs(Results.InterpOrgHighRes(:,2) - (a*Results.InterpOrgHighRes(:,1)+bY))<tolerance)*Results.SRHighRes;

%% Sample Entropy
X=diff(Results.InterpOrg(:,1));
tolerance=0.2*std(X); 
for window=4
    fieldname=sprintf('XSampEn%i',window);
    Results.(fieldname)=SampEn( window, tolerance, X);
end

Y=diff(Results.InterpOrg(:,2));
tolerance=0.2*std(Y); 
for window=4
    fieldname=sprintf('YSampEn%i',window);
    Results.(fieldname)=SampEn( window, tolerance, Y);
end

%% X & Y Flips
[Results.Xflips, Results.MaxXflip]=count_flips_amit(Results.InterpOrg(:,1),3,Results.SR);
[Results.Yflips, Results.MaxYflip]=count_flips_amit(Results.InterpOrg(:,2),3,Results.SR);

