%% Process & Generate AllResults Files
clear all; close all; clc
addpath(genpath(pwd));
AllTasks={'Main task', 'Motor_localizer', 'Numerical_localizer'};

BadSubjects=[309, 109];

%%change subjects to full strings, so no multiple identifications
for tp=1:length(AllTasks)
    TaskType=AllTasks{tp}; 
    fprintf('Task: %s\n',TaskType);
    
    %Load Choice Files
    pathChoice=['data\' TaskType '\Trial_level_data_and_parameters\'];
    dirChoice=dir([pathChoice '*.xls']);
    c = struct2cell(dirChoice);
    ChoiceFiles=c(1,:)';
    temp=cat(1,ChoiceFiles{:});
    Subjects=str2num(temp(:,1:3));

    %Load Mouse Tracking Files
    pathMT=['data\' TaskType '\Mouse_tracking\'];
    dirMT=dir([pathMT '/*.xls']);
    c = struct2cell(dirMT);
    MTFiles=c(1,:)';

    %Set MT Parameters
    MTtrialsCol=6;
    MTtimesCol=5;
    MTrajCols=1:2;

    %Set Task Parameters
    switch TaskType
        case 'Main task'
            params.PredCols=[25 26];
            params.a_col=4;
            params.b_col=5;
        case 'Motor_localizer'
            params.PredCols=[12 13];
            params.a_col=3;
            params.b_col=4;
        case 'Numerical_localizer'
            params.PredCols=[12 13];
            params.a_col=3;
            params.b_col=4;
    end     
    clear c temp

    % Load & Process MTs
    AllResults=cell(length(Subjects),3);
    for s=1:length(Subjects)
        %Subject Choice Data
        sub=Subjects(s);
        fprintf('Iter %d, Subject %d\n',s,sub);
        subChoice=ChoiceFiles{contains(ChoiceFiles,num2str(sub))};

        %Check if subject has Mouse Data & Load or if Bad Subject
        if sum(contains(MTFiles,num2str(sub)))==0 || ismember(sub,BadSubjects)
            fprintf('Subject %d has no/bad Mouse Data\n',sub);
            AllResults{s,1}=sub;
            AllResults{s,3}=ChoiceData;
            continue
        end
        MTFilesSubs = cellfun(@(x)(x(1:4)), MTFiles, 'UniformOutput', false);
        subMT=MTFiles{contains(MTFilesSubs,num2str(sub))};

        %Load Datasets
        [~, ~, ChoiceData]=xlsread([pathChoice subChoice]);
        [~, ~, MTData]=xlsread([pathMT subMT]);
        trials=cell2mat(MTData(:,MTtrialsCol));

        % Delete Training
        traininds=[];
        if min(trials)==-9
            traininds=ismember(trials,-9:-1);
        elseif max(trials)==159
            trials=trials-9;
            traininds=ismember(trials,-9:-1);
        end
        MTData(traininds,:)=[];
        trials(traininds)=[];

        %Collect Trajetories & Times
        trajs=cell2mat(MTData(:,MTrajCols));
        times=cell2mat(MTData(:,MTtimesCol));

        % Analyze
        clear SubResults
        for t=1:max(trials)
            trajectory=trajs(t==trials,:);
            time=times(t==trials,:);
            trialdata=ChoiceData(t+1,:);

            %Skip if Traj too short
            if size(trajectory,1)<2
                allfields=fieldnames(SubResults);
                for f=1:length(allfields)
                    SubResults(t).(allfields{f})=NaN;
                end
                continue
            end

            %Extract Mouse Features
            try
                SubResults(t)=trajectory_analysis(trajectory,time,trialdata,params);
            catch
                allfields=fieldnames(SubResults);
                for f=1:length(allfields)
                    SubResults(t).(allfields{f})=NaN;
                end
                continue
            end
        end
        AllResults{s,1}=sub;
        AllResults{s,2}=SubResults;
        ChoiceData(1,:)=strrep(ChoiceData(1,:),' (sec)','');
        ChoiceData(1,:)=strrep(ChoiceData(1,:),' ','_');
        ChoiceData(1,:)=strrep(ChoiceData(1,:),'-','_');
        AllResults{s,3}=cell2struct(ChoiceData(2:end,:), ChoiceData(1,:), 2);
    end
    
    try
    AllResults=addTimeInMinViolation(AllResults);
    catch
    end
    %clear empty
    emptyinds=cellfun(@isempty,AllResults(:,2));
    AllResults(emptyinds,:)=[];
    
    save(['AllResults-' TaskType '.mat'],'AllResults');
end
disp('Done!');