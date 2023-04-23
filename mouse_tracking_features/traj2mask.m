function Area=traj2mask(trajectory)

%% Prepare Trajectory
x=round(trajectory(:,1));
y=round(trajectory(:,2));
x(x<=0)=1;
y(y<=0)=1;
x(x>100)=100;
y(y>100)=100;

%% Turn to Mask
mask=zeros(100,100);
for k=1:length(x)
    mask(x(k),y(k))=1;
end

%% Completed Trajectory (Line Start-End)
straightX=round(linspace(x(1),x(end),1000));
straighty=round(linspace(y(1),y(end),1000));
for k=1:length(straightX)
    mask(straightX(k),straighty(k))=1;
end

%% Close Area
% Fill Edge Gaps
edgedmask = filledgegaps(mask, 4);

% Close Object in Mask
se = strel('disk',7);
J = imclose(edgedmask,se);

% Close Object in Mask
J =imfill(J,'holes');

Area=mean2(J);
