% plot inconsistency indices, and compare with previous samples

clear
close all

fileName = 'subject_level_inconsistency_indices.xlsx';

KPWL_data = xlsread(fileName, 'KPWL', 'A2:F36');
KHL_data = xlsread(fileName, 'KHL', 'A2:F92');
KHL_MRI_data = xlsread(fileName, 'KHL MRI_final', 'A2:F43');
KHL_Replication_data = xlsread(fileName, 'replication', 'A2:D70');
RandomDM_data = xlsread(fileName, 'RandomDM', 'A2:E1001');
RandomDM_75_data = xlsread(fileName, 'RandomDM_75', 'A2:F1001');
RandomDM_108_data = xlsread(fileName, 'RandomDM_108', 'A2:D1001');

RGB1 = [150,150,150]/255;
RGB2 = [114,47,55]/255;
RGB3 = [181,96,106]/255;
RGB4 = [107,150,181]/255;
RGB5 = [209,143,85]/255;
RGB6 = [111,163,122]/255;



%% GARP violations
% GARP_viol_Choi = Choi_data(:,2);
GARP_viol_KPWL = KPWL_data(:,2);
GARP_viol_KHL = KHL_data(:,2);
GARP_viol_KHL_MRI = KHL_MRI_data(:,2);
GARP_viol_KHL_Replication = KHL_Replication_data(:,2);
GARP_viol_Random = RandomDM_data(:,2);
GARP_viol_Random_75 = RandomDM_75_data(:,2);
GARP_viol_Random_108 = RandomDM_108_data(:,2);

% group = [ones(size(GARP_viol_Choi));
group = [ones(size(GARP_viol_KPWL)); 
            2 * ones(size(GARP_viol_KHL));
            3 * ones(size(GARP_viol_KHL_MRI));
            4 * ones(size(GARP_viol_KHL_Replication));
            5 * ones(size(GARP_viol_Random));
            6 * ones(size(GARP_viol_Random_75))];

% GARP_violations = {GARP_viol_Choi; GARP_viol_KPWL; GARP_viol_KHL; GARP_viol_Random}; 
GARP_violations = {GARP_viol_KPWL; GARP_viol_KHL; GARP_viol_KHL_MRI;GARP_viol_KHL_Replication; GARP_viol_Random; GARP_viol_Random_75}; 

xCenter = 1:numel(GARP_violations); 
spread = 0.5; % 0=no spread; 0.5=random spread within box bounds (can be any value)
figure()
subplot(1,3,1);
for i = 1:numel(GARP_violations)
     plot(rand(size(GARP_violations{i}))*spread -(spread/2) + xCenter(i), GARP_violations{i}, 'mo','linewidth', 0.5,...    
        'MarkerEdgeColor',[0.8 0.8 0.8]);
    hold on
end
hold on          

GARP_violations_double = cell2mat(GARP_violations);
box = boxplotCsub(GARP_violations_double, group,[],[],1,1,'k',1,1,1,[1 1],[],[],[]);
    set(box(8,6),'FaceColor',RGB6);
     set(box(8,5),'FaceColor',RGB5);
    set(box(8,4),'FaceColor',RGB4);
     set(box(8,3),'FaceColor',RGB3);
    set(box(8,2),'FaceColor',RGB2);
    set(box(8,1),'FaceColor',RGB1);

title('GARP violations');
ylim([0 11000]);
[p1(1),h1(1),stats1(1)] = ranksum(GARP_viol_KHL, GARP_viol_Random,'tail','left');
[p2(1),h2(1),stats2(1)] = ranksum(GARP_viol_KHL_MRI, GARP_viol_Random_75,'tail','left');
[p3(1),h3(1),stats3(1)] = ranksum(GARP_viol_KHL_MRI, GARP_viol_KPWL);
[p4(1),h4(1),stats4(1)] = ranksum(GARP_viol_KHL, GARP_viol_KPWL);
[p5(1),h5(1),stats5(1)] = ranksum(GARP_viol_KHL, GARP_viol_KHL_MRI);
[p6(1),h6(1),stats6(1)] = ranksum(GARP_viol_Random, GARP_viol_Random_108);
[p7(1),h7(1),stats7(1)] = ranksum(GARP_viol_Random_75, GARP_viol_Random_108);

% ES = Z/sqrt(n1+n2)
es1(1) = abs(stats1(1).zval)./(sqrt(length(GARP_viol_KHL)+length(GARP_viol_Random)));
es2(1) = abs(stats2(1).zval)./(sqrt(length(GARP_viol_KHL_MRI)+length(GARP_viol_Random_75)));
es3(1) = abs(stats3(1).zval)./(sqrt(length(GARP_viol_KHL_MRI)+length(GARP_viol_KPWL)));
es4(1) = abs(stats4(1).zval)./(sqrt(length(GARP_viol_KHL)+length(GARP_viol_KPWL)));
es5(1) = abs(stats5(1).zval)./(sqrt(length(GARP_viol_KHL)+length(GARP_viol_KHL_MRI)));
es6(1) = abs(stats6(1).zval)./(sqrt(length(GARP_viol_Random)+length(GARP_viol_Random_108)));
es7(1) = abs(stats7(1).zval)./(sqrt(length(GARP_viol_Random_75)+length(GARP_viol_Random_108)));


%% Afriat Index
Afriat_KPWL = KPWL_data(:,3);
Afriat_KHL = KHL_data(:,3);
Afriat_KHL_MRI = KHL_MRI_data(:,3);
Afriat_KHL_Replication = KHL_Replication_data(:,3);
Afriat_Random = RandomDM_data(:,3);
Afriat_Random_75 = RandomDM_75_data(:,3);
Afriat_Random_108 = RandomDM_108_data(:,3);

subplot(1,3,2);
Afriat_Index = {Afriat_KPWL; Afriat_KHL; Afriat_KHL_MRI; Afriat_KHL_Replication; Afriat_Random; Afriat_Random_75};
for i = 1:numel(Afriat_Index)
    plot(rand(size(Afriat_Index{i}))*spread -(spread/2) + xCenter(i), Afriat_Index{i}, 'mo','linewidth', 0.5,...
        'MarkerEdgeColor',[0.8 0.8 0.8]);
    hold on
end

hold on          

Afriat_Index_double = cell2mat(Afriat_Index);
box = boxplotCsub(Afriat_Index_double, group,[],[],1,1,'k',1,1,1,[1 1],[],[],[]);
    set(box(8,6),'FaceColor',RGB6);
    set(box(8,5),'FaceColor',RGB5);
    set(box(8,4),'FaceColor',RGB4);
    set(box(8,3),'FaceColor',RGB3);
    set(box(8,2),'FaceColor',RGB2);
    set(box(8,1),'FaceColor',RGB1);
                 
title('Afriat Index');
ylim([0 1]);
[p1(2),h1(2),stats1(2)] = ranksum(Afriat_KHL, Afriat_Random,'tail','left');
[p2(2),h2(2),stats2(2)] = ranksum(Afriat_KHL_MRI, Afriat_Random_75,'tail','left');
[p3(2),h3(2),stats3(2)] = ranksum(Afriat_KHL_MRI, Afriat_KPWL);
[p4(2),h4(2),stats4(2)] = ranksum(Afriat_KHL, Afriat_KPWL);
[p5(2),h5(2),stats5(2)] = ranksum(Afriat_KHL, Afriat_KHL_MRI);
[p6(2),h6(2),stats6(2)] = ranksum(Afriat_Random, Afriat_Random_108);
[p7(2),h7(2),stats7(2)] = ranksum(Afriat_Random_75, Afriat_Random_108);

% ES = Z/sqrt(n1+n2)
es1(2) = abs(stats1(2).zval)./(sqrt(length(GARP_viol_KHL)+length(GARP_viol_Random)));
es2(2) = abs(stats2(2).zval)./(sqrt(length(GARP_viol_KHL_MRI)+length(GARP_viol_Random_75)));
es3(2) = abs(stats3(2).zval)./(sqrt(length(GARP_viol_KHL_MRI)+length(GARP_viol_KPWL)));
es4(2) = abs(stats4(2).zval)./(sqrt(length(GARP_viol_KHL)+length(GARP_viol_KPWL)));
es5(2) = abs(stats5(2).zval)./(sqrt(length(GARP_viol_KHL)+length(GARP_viol_KHL_MRI)));
es6(2) = abs(stats6(2).zval)./(sqrt(length(GARP_viol_Random)+length(GARP_viol_Random_108)));
es7(2) = abs(stats7(2).zval)./(sqrt(length(GARP_viol_Random_75)+length(GARP_viol_Random_108)));


%% MMI
MMI_KPWL = KPWL_data(:,6);
MMI_KHL = KHL_data(:,6);
MMI_KHL_MRI = KHL_MRI_data(:,6);
MMI_KHL_Replication = KHL_Replication_data(:,4);
MMI_Random = RandomDM_data(:,4);
MMI_Random_75 = RandomDM_75_data(:,6);
MMI_Random_108 = RandomDM_108_data(:,4);

MMI = {MMI_KPWL; MMI_KHL; MMI_KHL_MRI; MMI_KHL_Replication; MMI_Random; MMI_Random_75};
subplot(1,3,3);
for i = 1:numel(MMI)
    plot(rand(size(MMI{i}))*spread -(spread/2) + xCenter(i), MMI{i}, 'mo','linewidth', 0.5,...
        'MarkerEdgeColor',[0.8 0.8 0.8]);
    hold on
end

hold on  

MMI_double = cell2mat(MMI);
box = boxplotCsub(MMI_double, group,[],[],1,1,'k',1,1,1,[1 1],[],[],[]);
     set(box(8,6),'FaceColor',RGB6);
     set(box(8,5),'FaceColor',RGB5);
     set(box(8,4),'FaceColor',RGB4);
     set(box(8,3),'FaceColor',RGB3);
     set(box(8,2),'FaceColor',RGB2);
     set(box(8,1),'FaceColor',RGB1);
     
title('MMI');
ylim([0 0.35]);

[p1(3),h1(3),stats1(3)] = ranksum(MMI_KHL, MMI_Random,'tail','left');
[p2(3),h2(3),stats2(3)] = ranksum(MMI_KHL_MRI, MMI_Random_75,'tail','left');
[p3(3),h3(3),stats3(3)] = ranksum(MMI_KHL_MRI, MMI_KPWL);
[p4(3),h4(3),stats4(3)] = ranksum(MMI_KHL, MMI_KPWL);
[p5(3),h5(3),stats5(3)] = ranksum(MMI_KHL, MMI_KHL_MRI);
[p6(3),h6(3),stats6(3)] = ranksum(MMI_Random, MMI_Random_108);
[p7(3),h7(3),stats7(3)] = ranksum(MMI_Random_75, MMI_Random_108);

% ES = Z/sqrt(n1+n2)
es1(3) = abs(stats1(3).zval)./(sqrt(length(GARP_viol_KHL)+length(GARP_viol_Random)));
es2(3) = abs(stats2(3).zval)./(sqrt(length(GARP_viol_KHL_MRI)+length(GARP_viol_Random_75)));
es3(3) = abs(stats3(3).zval)./(sqrt(length(GARP_viol_KHL_MRI)+length(GARP_viol_KPWL)));
es4(3) = abs(stats4(3).zval)./(sqrt(length(GARP_viol_KHL)+length(GARP_viol_KPWL)));
es5(3) = abs(stats5(3).zval)./(sqrt(length(GARP_viol_KHL)+length(GARP_viol_KHL_MRI)));
es6(3) = abs(stats6(3).zval)./(sqrt(length(GARP_viol_Random)+length(GARP_viol_Random_108)));
es7(3) = abs(stats7(3).zval)./(sqrt(length(GARP_viol_Random_75)+length(GARP_viol_Random_108)));

%% FDR
[FDR, Q] = mafdr([p3 p4 p5]);

%% histograms of GARP violations (supplementary materials)
% first we need to log transform the distribution of GARP violations
% so we will replace 0 with epsilon, to avoid -inf after the transformation
figure;
subplot(2,3,1);
edges = 0:0.01:0.3;
histogram(Afriat_KHL,edges,'Normalization','probability','FaceColor',RGB2);
ylabel('% of subjects'); xlabel('AI'); title('behavioral study'); xticks(0:0.1:0.3); ylim([0 0.5]);

subplot(2,3,2);
edges2 = 0:0.007:0.2;
histogram(Afriat_KHL_MRI,edges2,'Normalization','probability','FaceColor',RGB3);
xlabel('AI'); title('neuroimaging study'); xticks(0:0.1:0.2);  ylim([0 0.5]);

subplot(2,3,3);
edges2 = 0:0.007:0.2;
histogram(Afriat_KHL_Replication,edges2,'Normalization','probability','FaceColor',RGB3);
xlabel('AI'); title('replication study'); xticks(0:0.1:0.2);  ylim([0 0.5]);

subplot(2,3,4);
histogram(MMI_KHL,edges,'Normalization','probability','FaceColor',RGB2);
ylabel('% of subjects'); xlabel('subject-level MMI');  xticks(0:0.1:0.3); ylim([0 0.5]);

subplot(2,3,5);
histogram(MMI_KHL_MRI,edges2,'Normalization','probability','FaceColor',RGB3);
xlabel('subject-level MMI');  xticks(0:0.1:0.2); ylim([0 0.5]);

subplot(2,3,6);
histogram(MMI_KHL_Replication,edges2,'Normalization','probability','FaceColor',RGB3);
xlabel('subject-level MMI');  xticks(0:0.1:0.2); ylim([0 0.5]);

%% scatter plots correlating MMI with AI and GARP violations
figure;
subplot(2,3,1); 
scatter(MMI_KHL,Afriat_KHL,'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',RGB2);
hold on
h1 = lsline;
h1.LineWidth = 1.5; h1.LineStyle = ':'; h1.Color = RGB2;
xlabel('MMI'); ylabel('AI');
title('Behavioral study');
corr(MMI_KHL,Afriat_KHL);

subplot(2,3,2); 
scatter(MMI_KHL_MRI,Afriat_KHL_MRI,'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',RGB3);
hold on
h3 = lsline;
h3.LineWidth = 1.5; h3.LineStyle = ':'; h3.Color = RGB3;
xlabel('MMI'); ylabel('AI');
title('Neuroimaging study');
corr(MMI_KHL_MRI,Afriat_KHL_MRI);

subplot(2,3,3); 
scatter(MMI_KHL_Replication,Afriat_KHL_Replication,'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',RGB3);
hold on
h3 = lsline;
h3.LineWidth = 1.5; h3.LineStyle = ':'; h3.Color = RGB3;
xlabel('MMI'); ylabel('AI');
title('Relication study');
corr(MMI_KHL_MRI,Afriat_KHL_MRI);


subplot(2,3,4); 
scatter(MMI_KHL,GARP_viol_KHL,'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',RGB2);
hold on
h2 = lsline;
h2.LineWidth = 1.5; h2.LineStyle = ':'; h2.Color = RGB2;
xlabel('MMI'); ylabel('GARP violations (#)');
corr(MMI_KHL,GARP_viol_KHL);

subplot(2,3,5); 
scatter(MMI_KHL_MRI,GARP_viol_KHL_MRI,'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',RGB3);
xlabel('MMI'); ylabel('GARP violations (#)'); ylim([0 40]);
hold on
axes('Position',[.53 .37 .08 .08]);
box on
scatter(MMI_KHL_MRI,GARP_viol_KHL_MRI,'MarkerEdgeColor',[1 1 1], ...
        'MarkerFaceColor',RGB3); ylim([40 2000]); xlim([0 0.3]);
        xlabel('MMI'); ylabel('GARP viol.'); yticks([40 1000 2000]);
corr(MMI_KHL_MRI,GARP_viol_KHL_MRI);


subplot(2,3,6); 
scatter(MMI_KHL_Replication,GARP_viol_KHL_Replication,'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',RGB3);
hold on
h3 = lsline;
h3.LineWidth = 1.5; h3.LineStyle = ':'; h3.Color = RGB3;
xlabel('MMI'); ylabel('GARP violations (#)'); ylim([0 8000]);
corr(MMI_KHL_MRI,GARP_viol_KHL_MRI);