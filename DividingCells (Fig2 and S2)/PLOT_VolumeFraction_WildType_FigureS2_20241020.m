% Plotting Volume Fraction in Dense phase during division for Figure 1 

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
    Y{k} = y;%/max(y); %If you want to normalize the data by the maximum condensed volume fraction to mitigate bleaching and diffraction effects
    % TIME{k} = time;
    TIME{k} = [time(1:6), time(7:end)+5];
    DIVISIONFRAMESMIN{k} = DivisionFramesMin;
    clearvars -except Y TIME DIVISIONFRAMESMIN DataPath DataSubPaths DataPath1
end

%%

figure9 = figure('color',[1 1 1]);
hold on

best = [1 2 3  5 6 ];

colors = redpeachblue(size(best,2));

for i = 1:size(best,2)
    k = best(i);
    time = TIME{k};
    y = Y{k};
    plot(time,y, 'color',colors(i,:),'linewidth',2);
end

title('Volume Fraction')
xlabel('Time (min)')
ylabel('V(dense phase)/V(cell)')
set(gca,'fontsize',14);
legend

xlim([-100 100])
ylim([0 0.04])