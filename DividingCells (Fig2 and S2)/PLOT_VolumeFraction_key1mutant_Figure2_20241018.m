%% Zeiss data 
% key1;EPYC1-Venus
% Cell 1 is the one shown in the Figure and in the SI movie

DataPath1 = 'E:\PrincetonData\20220117';

CellNumbers = [1 2];

for i = 1:length(CellNumbers)
    k = CellNumbers(i);
    datasubpath = ['\k4;EPYC1_dividing_cell',num2str(k)];
    DataSubPaths{i} = datasubpath;
end

Y = {};
TIME = {};
DIVISIONFRAMESMIN = {};

for k = 1:size(DataSubPaths,2)
    PlotDataPath = [DataPath1,DataSubPaths{k},'\plots_20240708_exp2'];
    CellDataPath = [DataPath1,DataSubPaths{k}];
    
    load([PlotDataPath,'\plot_data.mat']);
    load([CellDataPath,'\PhysicalParameters.mat']);
    y = TOTAL_VOLUME./CELL_VOLUME;
   % DivisionFramesMin = DivisionFrames.*FrameInterval;
    DivisionFramesMin = time(DivisionFrames);
    Y{k} = y;%/max(y);
    % TIME{k} = time;
    TIME{k} = time;
    DIVISIONFRAMESMIN{k} = DivisionFramesMin;
    clearvars -except Y TIME DIVISIONFRAMESMIN DataPath DataSubPaths DataPath1
end


figure9 = figure('color',[1 1 1]);
hold on

best = [1 2];

%colors = redpeachblue(size(best,2));
colors = gray(size(best,2)+3);

k = best(1);
time = TIME{k};
y = Y{k}%./max(Y{k}(10:end));
plot(time(10:end-55),y(10:end-55),'linewidth',2);

hold on

k = best(2);
time = TIME{k};
y = Y{k}%./max(Y{k}(10:end));
plot(time(10:end),y(10:end),'linewidth',2);

title('Volume Fraction')
xlabel('Time (min)')
ylabel('V(dense phase)/V(cell)')
set(gca,'fontsize',14);
%legend

%ylim([0 1])
xlim([-200 200]) 
ylim([0 0.07])

%% Colton's data iSIM
DataPath1 = 'E:\PrincetonData\fromColton\20231126';

CellNumbers = [2 3 4 5 6 7 8 9 10 11 12];

for i = 1:length(CellNumbers)
    k = CellNumbers(i);
    datasubpath = ['\Cell',num2str(k),'_register'];
    DataSubPaths{i} = datasubpath;
end

Y = {};
TIME = {};
DIVISIONFRAMESMIN = {};

for k = 1:size(DataSubPaths,2)
    PlotDataPath = [DataPath1,DataSubPaths{k},'\plots_20240708_exp2'];
    CellDataPath = [DataPath1,DataSubPaths{k}];
    
    load([PlotDataPath,'\plot_data.mat']);
    load([CellDataPath,'\PhysicalParameters.mat']);
    y = TOTAL_VOLUME./CELL_VOLUME;
   % DivisionFramesMin = DivisionFrames.*FrameInterval;
    DivisionFramesMin = time(DivisionFrames);
    Y{k} = y;%/max(y);
    % TIME{k} = time;
    TIME{k} = time;
    DIVISIONFRAMESMIN{k} = DivisionFramesMin;
    clearvars -except Y TIME DIVISIONFRAMESMIN DataPath DataSubPaths DataPath1
end

k = 12;

PlotDataPath = ['E:\PrincetonData\fromColton\20231126\Cell1_register\plots_20240708_exp2'];
CellDataPath = ['E:\PrincetonData\fromColton\20231126\Cell1_register'];

load([PlotDataPath,'\plot_data.mat']);
load([CellDataPath,'\PhysicalParameters.mat']);
y = TOTAL_VOLUME./CELL_VOLUME;
% DivisionFramesMin = DivisionFrames.*FrameInterval;
DivisionFramesMin = time(DivisionFrames);
Y{k} = y;%/max(y);
% TIME{k} = time;
TIME{k} = time;
DIVISIONFRAMESMIN{k} = DivisionFramesMin;
clearvars -except Y TIME DIVISIONFRAMESMIN DataPath DataSubPaths DataPath1

% k = 13;
% 
% PlotDataPath = ['E:\PrincetonData\fromColton\20231126\Cell13_not-dividing\plots_20240708_exp2'];
% CellDataPath = ['E:\PrincetonData\fromColton\20231126\Cell13_not-dividing'];
% 
% load([PlotDataPath,'\plot_data.mat']);
% load([CellDataPath,'\PhysicalParameters.mat']);
% y = TOTAL_VOLUME./CELL_VOLUME;
% % DivisionFramesMin = DivisionFrames.*FrameInterval;
% DivisionFramesMin = time(DivisionFrames);
% Y{k} = y;%/max(y);
% % TIME{k} = time;
% TIME{k} = time;
% DIVISIONFRAMESMIN{k} = DivisionFramesMin;
% clearvars -except Y TIME DIVISIONFRAMESMIN DataPath DataSubPaths DataPath1


figure9 = figure('color',[1 1 1]);
hold on

best = [1 2 8];% 4 5  7 8 9 10 11 ];

best = [1 2 3 4 5 6 7 8 9 10  ];

%colors = redpeachblue(size(best,2));
colors = gray(size(best,2)+3);

for i = 1:size(best,2)
    k = best(i);
    time = TIME{k};
    y = Y{k};%./max(Y{k});
    plot(time,y, 'color',colors(i,:),'linewidth',2);
end

% No division in this cell
% hold on
% k=13;
% time = TIME{k};
% y = Y{k};
% plot(time+50,y,'color','r','linewidth',2)


title('Volume Fraction')
xlabel('Time (min)')
ylabel('V(dense phase)/V(cell)')
set(gca,'fontsize',14);
%legend

ylim([0 0.07])
xlim([-200 200])

%legend
