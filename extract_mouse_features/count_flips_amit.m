function [nFlips,maxFlip] = count_flips_amit( vec,parameter,SR )
%Returns the number of flips which are bigger than parameter (default
%parameter is 3)
if ~exist ('parameter','var')
    parameter = 3;
end
%vec2 will have 0 where there is no change, 1 where ther is positive change
%and -1 where there is negative change
vec2 = zeros(length(vec),1);
for i=1:length(vec)-1
    if vec(i)== vec(i+1)
        vec2(i)=0;
    end
    if vec(i) < vec(i+1)
       vec2(i) = 1;
    end
    if (vec(i) > vec(i+1))   
       vec2(i) = -1;
    end
end
direction = vec2(1);

%   1 in vec3 = flip (from positive to negative or from negative to
%   positive). if there is positive change and then no change (0 change) and then positive it
%   will not be considered a flip. if there is positive and then no change and then negative
%   it will be considered flip.
vec2(length(vec2)) = 0;
vec3 = zeros(length(vec2),1);
for i=1:length(vec2)-1
    if vec2(i)~=0
        if vec2(i)~=direction
            direction = vec2(i);
            vec3(i) = 1;
        end
    end
end
%the parameter is the minimum size of a flip that will be taken into account


indx = find (vec3);
maxFlip = 0;
nFlips = 0;
for i=1:length(indx)-1
    a = indx(i+1)-indx(i);
    if a>=parameter
        nFlips = nFlips+1;
        if a>maxFlip
            maxFlip = a;
        end
    end   
end

%Convert Number of Samples to Real Time
maxFlip=maxFlip*SR;
%returns the number of flips and the max flip.
%[res,maxflip];

end


