% Plotting average volume fraction at division for key1 mutant and
% wild-type cells for Figure 2

Title = 'Box and Whisker Volume Fraction';

%%
% KEY1 rescue

DataPath1 = 'E:\PrincetonData\20240410\KEY1_rescue';

DataSubPaths = {'\Cell_1_register','\Cell_2_register','\Cell_3_register','\Cell_4_register',...
    '\Cell_5_register','\Cell_6_register','\Cell_7_register','\Cell_8_register','\Cell_9_register','\Cell_10_register'};

Y = {};
TIME = {};
DIVISIONFRAMESMIN = {};

for k = 1:size(DataSubPaths,2)
    PlotDataPath = [DataPath1,DataSubPaths{k},'\plots_20240625'];
    CellDataPath = [DataPath1,DataSubPaths{k}];
    
    load([PlotDataPath,'\plot_data.mat']);
    load([CellDataPath,'\PhysicalParameters.mat']);
    y = TOTAL_VOLUME./CELL_VOLUME;
   % DivisionFramesMin = DivisionFrames.*FrameInterval;
    DivisionFramesMin = time(DivisionFrames);
    Y{k} = y;%/max(y);
     TIME{k} = time;
   % TIME{k} = [time(1:6), time(7:end)+5];
    DIVISIONFRAMESMIN{k} = DivisionFramesMin;
    clearvars -except Y TIME DIVISIONFRAMESMIN DataPath DataSubPaths DataPath1 
end

Rescue_VolFrac = [];

for k = 1:size(Y,2)
    idx = find(TIME{k}==DIVISIONFRAMESMIN{k}(3));
    vf = Y{k}(idx);
    Rescue_VolFrac = [Rescue_VolFrac, vf];
end

%%
DataPath1 = 'E:\PrincetonData\20220209\k4;RBCS1-Venus';

DataSubPaths = {'\Cell1','\Cell2','\Cell3','\Cell4','\Cell5','\Cell6','\Cell7'};


Y = {};
TIME = {};
DIVISIONFRAMESMIN = {};

for k = 1:size(DataSubPaths,2)
    PlotDataPath = [DataPath1,DataSubPaths{k},'\plots_20240419'];
    CellDataPath = [DataPath1,DataSubPaths{k}];
    
    load([PlotDataPath,'\plot_data.mat']);
    load([CellDataPath,'\PhysicalParameters.mat']);
    y = TOTAL_VOLUME./CELL_VOLUME;
   % DivisionFramesMin = DivisionFrames.*FrameInterval;
    DivisionFramesMin = time(DivisionFrames);
    Y{k} = y;%/max(y);
    % TIME{k} = time;
    TIME{k} = [time(1:6), time(7:end)+5];
    DIVISIONFRAMESMIN{k} = DivisionFramesMin;
    clearvars -except Y TIME DIVISIONFRAMESMIN DataPath DataSubPaths DataPath1 Rescue_VolFrac
end

k4_VolFrac = [];

for k = 1:size(Y,2)
    idx = find(TIME{k}==DIVISIONFRAMESMIN{k}(3));
    vf = Y{k}(idx);
    k4_VolFrac = [k4_VolFrac, vf];
end

%%
DataPath1 = 'E:\PrincetonData\20220209\WT;RBCS1-Venus';



DataSubPaths = {'\Cell1','\Cell2','\Cell3','\Cell4','\Cell5','\Cell6'};


Y = {};
TIME = {};
DIVISIONFRAMESMIN = {};

for k = 1:size(DataSubPaths,2)
    PlotDataPath = [DataPath1,DataSubPaths{k},'\plots_20240419'];
    CellDataPath = [DataPath1,DataSubPaths{k}];
    
    load([PlotDataPath,'\plot_data.mat']);
    load([CellDataPath,'\PhysicalParameters.mat']);
    y = TOTAL_VOLUME./CELL_VOLUME;
   % DivisionFramesMin = DivisionFrames.*FrameInterval;
    DivisionFramesMin = time(DivisionFrames);
    Y{k} = y;%/max(y);
    % TIME{k} = time;
    TIME{k} = [time(1:6), time(7:end)+5];
    DIVISIONFRAMESMIN{k} = DivisionFramesMin;
    clearvars -except Y TIME DIVISIONFRAMESMIN DataPath DataSubPaths DataPath1 k4_VolFrac Rescue_VolFrac
end

WT_VolFrac = [];

for k = 1:size(Y,2)
    idx = find(TIME{k}==0);
    vf = Y{k}(idx);
    WT_VolFrac = [WT_VolFrac, vf];
end
%%

best_rescue = [1 2 3 4 5 6 7 8 10];

best_k4 = [2 3 4 5 6 7];

best_WT = [1 2 3 4 5 6 ];

%%

%best_k4 = [2 3 4 5 6];

Y_ALL{1} = WT_VolFrac(best_WT);
Y_ALL{2} = k4_VolFrac(best_k4);
Y_ALL{3} = Rescue_VolFrac(best_rescue);

WTData = WT_VolFrac(best_WT);
key1Data = k4_VolFrac(best_k4);
rescueData = Rescue_VolFrac(best_rescue);;

best = [1 2 3];

FUNC_PlotBoxWhisker(Title,best,Y_ALL);
hold on
scatter(ones(size(WTData)),WTData,'filled','MarkerFaceAlpha',0.3,'MarkerFaceColor','blue')
scatter(ones(size(key1Data)).*2,key1Data,'filled','MarkerFaceAlpha',0.3,'MarkerFaceColor','blue')
scatter(ones(size(rescueData)).*3,rescueData,'filled','MarkerFaceAlpha',0.3,'MarkerFaceColor','blue')
scatter(1,mean(WTData),'filled','MarkerFaceColor','k')
scatter(2,mean(key1Data),'filled','MarkerFaceColor','k')
scatter(3,mean(rescueData),'filled','MarkerFaceColor','k')
%ylim([0 1.1])

%%
[~,p_value] = ttest2(key1Data,rescueData);
