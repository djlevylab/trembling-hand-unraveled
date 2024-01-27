%% plotting precentage of points to the Y-product as function of log price ratio (budget set's slope%%

clear
close all;

RGB2 = [114,47,55]/255;
RGB3 = [181,96,106]/255;
RGB4 = [107,150,181]/255;

%% behavioral study
fileName = 'BehavioralStudy_all_sessions_RawData.csv'; 
data = load(fileName); 

%% MRI study
fileNameMRI = 'MRIStudy_all_sessions_RawData.csv';
dataMRI = load(fileNameMRI);

%% Replication study
fileNameReplication = 'ReplicationStudy_all_sessions_RawData.csv';
dataReplication = load(fileNameReplication);

%% Share of tokens to the Y-account - behavioral study

priceRatio = data(:,9)./data(:,8);
logPriceRatio = log(priceRatio);
ShareTokensToY = data(:,7)./(data(:,6)+data(:,7)); %share of tokens to Y of total

% bin by slopes/price ratios
logPriceRatioQuantiles = discretize(logPriceRatio,[-inf -1.75 -1.5 -1.25 -1 -0.75 -0.5 -.025 0 0.25 0.5 0.75 1 1.25 1.5 1.75 inf]);

%% Share of tokens to the Y-account - MRI study

priceRatioMRI = dataMRI(:,6)./dataMRI(:,5);
logPriceRatioMRI = log(priceRatioMRI);
ShareTokensToY_MRI = dataMRI(:,4)./(dataMRI(:,3)+dataMRI(:,4)); %share of tokens to Y of total

% bin by slopes/price ratios
logPriceRatioQuantilesMRI = discretize(logPriceRatioMRI,[-inf -1.75 -1.5 -1.25 -1 -0.75 -0.5 -.025 0 0.25 0.5 0.75 1 1.25 1.5 1.75 inf]);

%% Share of tokens to the Y-account - Repliation study

priceRatioReplication = dataReplication(:,9)./dataReplication(:,8);
logPriceRatioReplication = log(priceRatioReplication);
ShareTokensToY_Replication = dataReplication(:,7)./(dataReplication(:,6)+dataReplication(:,7)); %share of tokens to Y of total

% bin by slopes/price ratios
logPriceRatioQuantilesReplication = discretize(logPriceRatioReplication,[-inf -1.75 -1.5 -1.25 -1 -0.75 -0.5 -.025 0 0.25 0.5 0.75 1 1.25 1.5 1.75 inf]);

%% %%%% figure - multiple boxplots %%%% %%

figure()
box = boxplotCsub(ShareTokensToY, logPriceRatioQuantiles,[],[],1,1,'k',1,0.5,1,[1 3],[],[],[]);
    set(box(8,:),'FaceColor',RGB2);
hold on
box2 = boxplotCsub(ShareTokensToY_MRI, logPriceRatioQuantilesMRI,[],[],1,1,'k',1,0.5,1,[2 3],[],[],[]);
    set(box2(8,:),'FaceColor',RGB3);
hold on
box3 = boxplotCsub(ShareTokensToY_Replication, logPriceRatioQuantilesReplication,[],[],1,1,'k',1,0.5,1,[3 3],[],[],[]);
    set(box3(8,:),'FaceColor',RGB4);

xlabel('log price ratio','FontName','Arial','FontSize',12);
ylabel('tokens share to Y','FontName','Arial','FontSize',12);
axis([0.5 14.5 0 1]);
xticks([0.5 7.25 14.5]);
xticklabels({'-2.5', '0', '2.5'});
hold on

