% Plotting Partition Coefficient during division for Wild Type EPYC1-Venus
% cells Figure 2
% 2024-10-18 LL

DataPath1 = 'E:\PrincetonData\20220103';

CellNumbers = [1 2 3 4 5 6 7 8];

for i = 1:length(CellNumbers)
    k = CellNumbers(i);
    datasubpath = ['\Cell',num2str(k)];
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
    Y{k} = y%/max(y);
    % TIME{k} = time;
    TIME{k} = time;
    DIVISIONFRAMESMIN{k} = DivisionFramesMin;
    clearvars -except Y TIME DIVISIONFRAMESMIN DataPath DataSubPaths DataPath1
end

k = 9;

PlotDataPath = ['E:\PrincetonData\20220103\Cell9-noDivision\plots_20240708_exp2'];
CellDataPath = ['E:\PrincetonData\20220103\Cell9-noDivision'];

load([PlotDataPath,'\plot_data.mat']);
load([CellDataPath,'\PhysicalParameters.mat']);
y = TOTAL_VOLUME./CELL_VOLUME;
% DivisionFramesMin = DivisionFrames.*FrameInterval;
DivisionFramesMin = time(DivisionFrames);
Y{k} = y;
% TIME{k} = time;
TIME{k} = time;
DIVISIONFRAMESMIN{k} = DivisionFramesMin;
clearvars -except Y TIME DIVISIONFRAMESMIN DataPath DataSubPaths DataPath1
%% three cells
figure9 = figure('color',[1 1 1]);
hold on

%best = [1 2 3 4 5 6 7 8 ];
best = [1 3 5];
%colors = redpeachblue(size(best,2));
colors = gray(size(best,2)+3);

for i = 1:size(best,2)
    k = best(i);
    time = TIME{k};
    y = Y{k};%/max(Y{k});
    plot(time(10:end),y(10:end), 'color',colors(i,:),'linewidth',2);
end

title('Volume Fraction')
xlabel('Time (min)')
ylabel('V(dense phase)/V(cell)')
set(gca,'fontsize',14);
%legend

xlim([-200 200])
%ylim([0 1])
%% five cells

figure9 = figure('color',[1 1 1]);
hold on

best = [1 3 4 5 6];
%colors = redpeachblue(size(best,2));
colors = gray(size(best,2)+3);

for i = 1:size(best,2)
    k = best(i);
    time = TIME{k};
    y = Y{k};%/max(Y{k});
    plot(time(10:end),y(10:end), 'color',colors(i,:),'linewidth',2);
end
title('Volume Fraction')
xlabel('Time (min)')
ylabel('V(dense phase)/V(cell)')
set(gca,'fontsize',14);
%legend

xlim([-200 200])
ylim([0 0.07])
legend

%% all cells

figure9 = figure('color',[1 1 1]);
hold on

best = [1 2 3 4 5 6 7 8 ];
%best = [1 3 4 5 6];
%colors = redpeachblue(size(best,2));
colors = gray(size(best,2)+3);

for i = 1:size(best,2)
    k = best(i);
    time = TIME{k};
    y = Y{k};%/max(Y{k});
    plot(time(10:end),y(10:end), 'color',colors(i,:),'linewidth',2);
end
title('Volume Fraction')
xlabel('Time (min)')
ylabel('V(dense phase)/V(cell)')
set(gca,'fontsize',14);
%legend

xlim([-200 200])
ylim([0 0.07])
legend

%% with a non-dividing cell

hold on
k=9;
time = TIME{k};
y = Y{k};
plot(time+100,y,'color','r','linewidth',2)
legend

title('Volume Fraction')
xlabel('Time (min)')
ylabel('V(dense phase)/V(cell)')
set(gca,'fontsize',14);
%legend




