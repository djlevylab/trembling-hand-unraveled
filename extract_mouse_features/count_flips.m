function count=count_flips(vec)

%Create vector where all positive progress is 1 and negative progress is 0.
vecdiff=diff(vec)>=0;

%Create vector where changes in sign is marked with 1.
SignChanges=abs(diff(vecdiff));

%Finish if there are no changes in sign
if sum(SignChanges)==0
    count=0;
    return
end

%Create Index vector for locations of flips that are larger than 3 bins.
indx=find(SignChanges);
locs=diff(indx)<3;

%Clear Sign Changes that are shorter than 3 bins.
if sum(locs)>=1
    SignChanges([indx(locs); indx(find(locs)+1)])=0;
end

%Count Flips
count=sum(SignChanges);