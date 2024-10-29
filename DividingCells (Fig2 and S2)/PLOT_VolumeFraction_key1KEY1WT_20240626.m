% Plotting Volume Fraction in Dense phase during division for Figure 2

DataPath1 = 'E:\PrincetonData\20221005';

DataSubPaths = {'\key1_no-rescue_cell_1','\KEY1_rescue_cell1_registered','\WT_Cell1_registered'};

Y = {};
TIME = {};
DIVISIONFRAMESMIN = {};

for k = 1:size(DataSubPaths,2)
    PlotDataPath = [DataPath1,DataSubPaths{k},'\plots_20240626'];
    CellDataPath = [DataPath1,DataSubPaths{k}];
    load([PlotDataPath,'\plot_data.mat']);
    load([CellDataPath,'\PhysicalParameters.mat']);
    y = TOTAL_VOLUME./CELL_VOLUME;
   % DivisionFramesMin = DivisionFrames.*FrameInterval;
    DivisionFramesMin = time(DivisionFrames);
    Y{k} = y;%./max(y);
    TIME{k} = time;
    %TIME{k} = [time(1:6), time(7:end)+5];
    DIVISIONFRAMESMIN{k} = DivisionFramesMin;
    clearvars -except Y TIME DIVISIONFRAMESMIN DataPath DataSubPaths DataPath1
end

best = [1 2 3];

figure1 = figure('color',[1 1 1]);
hold on
for i = 1:size(best,2)
    k = best(i);
    time = TIME{k};
    y = Y{k};
    plot(time,y, 'color',[0.5 0.5 0.5],'linewidth',2);
end

title('Volume Fraction')
xlabel('Time (min)')
ylabel('V(dense phase)/V(cell)')
set(gca,'fontsize',14);
xlim([-100 100]);
ylim([0 0.014]);