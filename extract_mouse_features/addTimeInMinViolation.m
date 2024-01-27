function AllResults=addTimeInMinViolation(AllResults)

filename='Residual Per Interval - Results - Date 19-Aug-2019 - 9513930843.xls';
[status,sheets] = xlsfinfo(filename);
tolerance=3; badsheets={}; badsubs=[];
for f=1:length(sheets)
    [num,txt,raw]=xlsread(filename,sheets{f});
    sub=raw{2,1};
    ind=find(sub==cell2mat(AllResults(:,1)));
    ChoiceData=AllResults{ind,3};
    if length(ChoiceData)~=size(num,1)
        fprintf('Bad Sheet: %d\n',sub);
        badsheets=[badsheets sheets(f)];
        badsubs=[badsubs sub];
        continue
    end
    
    TimeInMinViolation=NaN(length(ChoiceData),1);
    for t=1:length(ChoiceData)
        if num(t,3)~=0 %Change NaN only in trials with violations
            trajectory=AllResults{ind,2}(t).InterpOrgHighRes;
            SR=AllResults{ind,2}(t).SRHighRes;
            a=ChoiceData(t).slope;
            bY=ChoiceData(t).IntersectionY;
            col=5;
            while col<size(num,2)
                x1=min([num(t,col) num(t,col+2)]);
                x2=max([num(t,col) num(t,col+2)]);
                if isnan(x1)
                    break;
                end
                col=col+4;
                %Find samples that were around minimal violation segment
                locs=trajectory(:,1)>x1 & trajectory(:,1)<x2;
                RelevantCoords=trajectory(locs,:);
                FinalCoords=RelevantCoords(abs(RelevantCoords(:,2) - (a*RelevantCoords(:,1)+bY))<tolerance,:);
                TimeInMinViolation(t)=nansum([size(FinalCoords,1)*SR TimeInMinViolation(t)]);
                fprintf('Subject: %d, Trial: %d, Time: %f\n',sub,t,TimeInMinViolation(t))
            end
        end
        AllResults{ind,2}(t).TimeInMinViolation=TimeInMinViolation(t);
    end
end