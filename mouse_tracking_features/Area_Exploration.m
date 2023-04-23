close all;
figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf,'color',[0.9 0.9 0.9]);
sbpl=1;

subplot(2,4,sbpl)
x=round(Results.InterpOrg(:,1));
y=round(Results.InterpOrg(:,2));
scatter(x,y);
xlim([1 100]);
ylim([1 100]);
title('Original Trajectory');
sbpl=sbpl+1;

%% High Resolution Trajectory
subplot(2,4,sbpl)
x=round(Results.InterpOrgHighRes(:,1));
y=round(Results.InterpOrgHighRes(:,2));
x(x<=0)=1;
y(y<=0)=1;
x(x>100)=100;
y(y>100)=100;
scatter(x,y);
xlim([1 100]);
ylim([1 100]);
title('High Resolution Trajectory (No Zeros)');
sbpl=sbpl+1;
 
%% Trajectory as Mask Image
subplot(2,4,sbpl)
mask=zeros(100,100);
for k=1:length(x)
mask(x(k),y(k))=1;
end
imshow(imrotate(mask,90))
title('Trajectory as Mask Image');
sbpl=sbpl+1;
%% Completed Trajectory (Line Start-End)
subplot(2,4,4)
straightX=round(linspace(x(1),x(end),1000));
straighty=round(linspace(y(1),y(end),1000));
for k=1:length(straightX)
    mask(straightX(k),straighty(k))=1;
end
imshow(imrotate(mask,90))
title('Completed Trajectory (Line Start-End)');
sbpl=sbpl+1;

hold on
scatter([x(1) x(end)],[100-y(1) 100-y(end)],40,[0.4660 0.6740 0.1880; 0.6350 0.0780 0.1840],'filled')

%% Fill Edge Gaps
subplot(2,4,sbpl)
edgedmask = filledgegaps(mask, 4);
imshow(imrotate(edgedmask,90))
sbpl=sbpl+1;
title('Fill Edge Gaps');
hold on
scatter([x(1) x(end)],[100-y(1) 100-y(end)],40,[0.4660 0.6740 0.1880; 0.6350 0.0780 0.1840],'filled')

%% Close Object in Mask
subplot(2,4,sbpl)
se = strel('disk',7);
J = imclose(edgedmask,se);
imshow(imrotate(J,90))
title('Closed Trajectory Area)');
sbpl=sbpl+1;
hold on
scatter([x(1) x(end)],[100-y(1) 100-y(end)],40,[0.4660 0.6740 0.1880; 0.6350 0.0780 0.1840],'filled')

%% Close Object in Mask
subplot(2,4,sbpl)
J =imfill(J,'holes');
imshow(imrotate(J,90))
title('Fill Trajectory Area(holes)');
sbpl=sbpl+1;
hold on
scatter([x(1) x(end)],[100-y(1) 100-y(end)],40,[0.4660 0.6740 0.1880; 0.6350 0.0780 0.1840],'filled')

%% Complete Contour with Convex Hull
subplot(2,4,sbpl);
CH_objects = bwconvhull(J,'objects');
imshow(imrotate(CH_objects,90));
title('Objects Convex Hull - filled mask ');
sbpl=sbpl+1;
hold on
scatter([x(1) x(end)],[100-y(1) 100-y(end)],40,[0.4660 0.6740 0.1880; 0.6350 0.0780 0.1840],'filled')

%% Save
saveas(gcf,sprintf('Areas/Example-S%d-T%d.jpg',sub,trialnum));
