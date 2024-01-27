function d=get_MD(Trajectory)
%the max of the distance from the straight line between 0,0 and 1,1 (0 if
%its on the line)
N=4;
% Q1=[0,0]; Q2=[1,1];
Q1=[Trajectory(1,1),Trajectory(1,2)];
Q2=[Trajectory(end,1),Trajectory(end,2)];

ds=[];
for i=N+1:length(Trajectory)-N
    x1=Trajectory(i,1);
    y1=Trajectory(i,2);
    P=[x1 y1]; 
    ds(i)=abs(det([Q2-Q1,;P-Q1]))/norm(Q2-Q1);
end

d=max(ds);
end