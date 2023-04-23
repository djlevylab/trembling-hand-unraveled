function [new_traj, newtime]=clean_traj(trajectory,time)

MaxSample=min([8 size(trajectory,1)-1]);
MinDiff=14.5;
FirstSampleThresh=20;

if abs(trajectory(1,1))<FirstSampleThresh && abs(trajectory(1,2))<FirstSampleThresh
   new_traj=trajectory;
   newtime=time;
   return
end  

TrajDiffs=[diff(trajectory(:,1)), diff(trajectory(:,2))];
I = [find(abs(TrajDiffs(1:MaxSample,1))>MinDiff,1,'first') find(abs(TrajDiffs(1:MaxSample,2))>MinDiff,1,'first')];

if isempty(I)
   new_traj=trajectory;
   newtime=time;
   return
end

new_traj=trajectory((min(I)+1):end,:);
newtime=time((min(I)+1):end,:);

end