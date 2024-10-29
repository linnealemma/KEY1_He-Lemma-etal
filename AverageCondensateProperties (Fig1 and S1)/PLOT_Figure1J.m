% Code used to make Figure 5E
% 2024-01-22 LL

%% Volume of Largest
% From volume
Title = 'Volume of Largest';

DataPaths = {DataSubPath_wt1,DataSubPath_mutant1,DataSubPath_mutant2,DataSubPath_mutant3,DataSubPath_rescue1,DataSubPath_rescue2, DataSubPath_rescue3};

X = categorical({'Wild Type','key1 mutant1','key1 mutant2','key1 mutant3','KEY1 rescue1','KEY1 rescue2','KEY1 rescue3'});
X = reordercats(X,{'Wild Type','key1 mutant1','key1 mutant2','key1 mutant3','KEY1 rescue1','KEY1 rescue2','KEY1 rescue3'});

Y_ALL = cell(size(DataPaths));
Y = [];
Yerr = [];
for i = 1:length(DataPaths)
    load([DataPaths{i},'\PROPERTIES.mat'])
    load([DataPaths{i},'\CELL.mat'])
    load([DataPaths{i},'\PhysicalParameters.mat']);
    % Create a logical index for structures you want to keep
    indexToKeep = [PROPERTIES.number] ~= 0;
    
    % Use this index to create a new array
    PROPERTIES = PROPERTIES(indexToKeep);
    CELL = CELL(indexToKeep);
    y1 = [];
    a = [];
    for j = 1:size(PROPERTIES,2)
        y1 = [y1, max(PROPERTIES(j).volume)];
    end
    y = y1;
    y1(y1>21) = [];
    Y_ALL{i} = y1;
    Y = [Y,mean(y)];
    Yerr = [Yerr,std(y)];
end

Y_ALL_Largest = Y_ALL;
%% Volume of All Condensates 

Title = 'Volume of All Condensates'
DataPaths = {DataSubPath_wt1,DataSubPath_mutant1,DataSubPath_mutant2,DataSubPath_mutant3,DataSubPath_rescue1,DataSubPath_rescue2, DataSubPath_rescue3};

X = categorical({'Wild Type','key1 mutant1','key1 mutant2','key1 mutant3','KEY1 rescue1','KEY1 rescue2','KEY1 rescue3'});
X = reordercats(X,{'Wild Type','key1 mutant1','key1 mutant2','key1 mutant3','KEY1 rescue1','KEY1 rescue2','KEY1 rescue3'});

Y_ALL = cell(size(DataPaths));
Y = [];
Yerr = [];
for i = 1:length(DataPaths)
    load([DataPaths{i},'\PROPERTIES.mat'])
    load([DataPaths{i},'\CELL.mat'])
    load([DataPaths{i},'\PhysicalParameters.mat']);
    y1 = [];
    a = [];
    for j = 1:size(PROPERTIES,2)
        y1 = [y1, sum(PROPERTIES(j).volume)];
        a = [a, (CELL(j).bb(4)*CELL(j).bb(5)*CELL(j).bb(6))*voxel];
    end
    y = y1;
    y(y>21) = [];
    Y_ALL{i} = y;
    Y = [Y,mean(y)];%./mean(a);
    Yerr = [Yerr,std(y)];%./mean(a);
end

Y_ALL_all = Y_ALL;
%%
labels2 ={'Wild-type total','Wild-type largest','key1 mutant total','key1 mutant largest','KEY1 rescue total','KEY1 rescue largest'};

data = []; % This will hold the corresponding data
group=[];
c = 1;
for i = 1:length(best)
j = best(i);
y1 = reshape(Y_ALL_Largest{j},1,[]);
y2 = reshape(Y_ALL_all{j},1,[]);
data = [data, y1, y2]; % Append the data
group = [group, repmat(labels2(c+1), 1, length(y1))];
group = [group, repmat(labels2(c), 1, length(y2))];
c = c+2;
%group = [group, repmat(g2, 1, length(y2))]; % Append the label
end
figure3 = figure('color', [1 1 1]);
V = violinplot(data,group,'GroupOrder',labels2);
hold on
title(Title);
box

%% Percent Volume in the main condensate


Title = 'Percent Volume in the main condensate All Condensates'
DataPaths = {DataSubPath_wt1,DataSubPath_mutant1,DataSubPath_mutant2,DataSubPath_mutant3,DataSubPath_rescue1,DataSubPath_rescue2, DataSubPath_rescue3};

X = categorical({'Wild Type','key1 mutant1','key1 mutant2','key1 mutant3','KEY1 rescue1','KEY1 rescue2','KEY1 rescue3'});
X = reordercats(X,{'Wild Type','key1 mutant1','key1 mutant2','key1 mutant3','KEY1 rescue1','KEY1 rescue2','KEY1 rescue3'});

Y_ALL = cell(size(DataPaths));
Y = [];
Yerr = [];
for i = 1:length(DataPaths)
    load([DataPaths{i},'\PROPERTIES.mat'])
    load([DataPaths{i},'\CELL.mat'])
    load([DataPaths{i},'\PhysicalParameters.mat']);
    y1 = [];
    y2 = [];
    a = [];
    for j = 1:size(PROPERTIES,2)
        y1 = [y1, sum(PROPERTIES(j).volume)];
        y2 = [y2, max(PROPERTIES(j).volume)];
        if isempty(max(PROPERTIES(j).volume)) == 1
            y2 = [y2, 0];
        end
        a = [a, (CELL(j).bb(4)*CELL(j).bb(5)*CELL(j).bb(6))*voxel];
    end
    y = y2./y1*100;
    %y(y>21) = [];
    Y_ALL{i} = y;
    Y = [Y,mean(y)];%./mean(a);
    Yerr = [Yerr,std(y)];%./mean(a);
end

FUNC_PlotViolin(Title, best, Y_ALL, labels);
box

FUNC_PlotViolin(Title, [1 2 3 4 5 6 7], Y_ALL, labels);
box
