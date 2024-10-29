%%
DataPaths = {'E:\PrincetonData\20231103\wild-type\bin2\R0C0','E:\PrincetonData\20231103\wild-type\bin2\R0C1',...
    'E:\PrincetonData\20231103\KEY1-rescue\R0C0','E:\PrincetonData\20231103\KEY1-rescue\R1C1','E:\PrincetonData\20231103\KEY1-rescue\R0C1',...
    'E:\PrincetonData\20231103\woRBM\a5_bin2', 'E:\PrincetonData\20231103\woRBM\a5_2_bin2'};

%% Sum Intensity Inside condensate each

Title = 'Sum Intensity In Condensate';

X = categorical({'Wild-type 1', 'Wild-type 2', 'KEY1 rescue 1','KEY! rescue 2', 'KEY1 rescue 3','KEY1_woRBM rescue 1','KEY1_woRBM rescue 2'});
X = reordercats(X,{'Wild-type 1', 'Wild-type 2', 'KEY1 rescue 1','KEY! rescue 2', 'KEY1 rescue 3','KEY1_woRBM rescue 1','KEY1_woRBM rescue 2'});


Y_ALL = cell(size(DataPaths));
Y = [];
Yerr = [];
for i = 1:length(DataPaths)
    load([DataPaths{i},'\PROPERTIES_SNAP.mat']);
    load([DataPaths{i},'\CELL.mat']);
    load([DataPaths{i},'\PhysicalParameters.mat']);
    % Create a logical index for structures you want to keep
    indexToKeep = [PROPERTIES.number] ~= 0;
    
    % Use this index to create a new array
    PROPERTIES = PROPERTIES(indexToKeep);
    CELL = CELL(indexToKeep);
    y1 = [];
    a = [];
    for j = 1:size(PROPERTIES,2)
        y1 = [y1, sum(PROPERTIES(j).SumIntSNAPIn)];
    end
    y = y1;
    Y_ALL{i} = y;
    Y = [Y,mean(y)];
    Yerr = [Yerr,std(y)];
end

labels = {'Wild-type 1', 'Wild-type 2', 'KEY1 rescue 1','KEY! rescue 2', 'KEY1 rescue 3','KEY1_woRBM rescue 1','KEY1_woRBM rescue 2'};
best = [1 2 3 4 5 6 7];

FUNC_PlotBoxWhisker(Title,best,Y_ALL)
FUNC_PlotViolin(Title, best, Y_ALL, labels);
box

%% Sum Intensity In Cell each

Title = 'Sum Intensity';

X = categorical({'Wild-type 1', 'Wild-type 2', 'KEY1 rescue 1','KEY! rescue 2', 'KEY1 rescue 3','KEY1_woRBM rescue 1','KEY1_woRBM rescue 2'});
X = reordercats(X,{'Wild-type 1', 'Wild-type 2', 'KEY1 rescue 1','KEY! rescue 2', 'KEY1 rescue 3','KEY1_woRBM rescue 1','KEY1_woRBM rescue 2'});


Y_ALL = cell(size(DataPaths));
Y = [];
Yerr = [];
for i = 1:length(DataPaths)
    load([DataPaths{i},'\PROPERTIES_SNAP.mat'])
    % Create a logical index for structures you want to keep
    indexToKeep = [PROPERTIES.number] ~= 0;
    
    % Use this index to create a new array
    PROPERTIES = PROPERTIES(indexToKeep);

    y = cat(1,PROPERTIES.SumIntSNAPCell);
    Y_ALL{i} = y;
    Y = [Y,mean(y)];
    Yerr = [Yerr,std(y)];
end

labels = {'Wild-type 1', 'Wild-type 2', 'KEY1 rescue 1','KEY! rescue 2', 'KEY1 rescue 3','KEY1_woRBM rescue 1','KEY1_woRBM rescue 2'};
best = [1 2 3 4 5 6 7];

FUNC_PlotBoxWhisker(Title,best,Y_ALL)
FUNC_PlotViolin(Title, best, Y_ALL, labels);
box

%% Sum Intensity Inside condensate

Title = 'Sum Intensity In Condensate';

X = categorical({'Wild-type', 'KEY1 rescue','KEY1_woRBM rescue'});
X = reordercats(X,{'Wild-type','KEY1 rescue', 'KEY1_woRBM rescue'});

xIndex = [1 1 2 2 2 3 3];

%Y_ALL = cell(size(DataPaths));
Y_ALL = cell(1,3);
Y = [];
Yerr = [];
for i = 1:length(DataPaths)
    k = xIndex(i);
    load([DataPaths{i},'\PROPERTIES_SNAP.mat']);
    load([DataPaths{i},'\CELL.mat']);
    load([DataPaths{i},'\PhysicalParameters.mat']);
    % Create a logical index for structures you want to keep
    indexToKeep = [PROPERTIES.number] ~= 0;
    
    % Use this index to create a new array
    PROPERTIES = PROPERTIES(indexToKeep);
    CELL = CELL(indexToKeep);
    y1 = [];
    a = [];
    for j = 1:size(PROPERTIES,2)
        y1 = [y1, sum(PROPERTIES(j).SumIntSNAPIn)];
    end
    y = y1;
    Y_ALL{k} = [Y_ALL{k}, y];
    Y = [Y,mean(y)];
    Yerr = [Yerr,std(y)];
end

labels = {'Wild Type','KEY1 Rescue','KEY1_woRBM'};
best = [1 2 3];

FUNC_PlotBoxWhisker(Title,best,Y_ALL)
FUNC_PlotViolin(Title, best, Y_ALL, labels);
box


%% Sum Intensity In Cell

Title = 'Sum Intensity';

X = categorical({'Wild-type', 'KEY1 rescue','KEY1_woRBM rescue'});
X = reordercats(X,{'Wild-type','KEY1 rescue', 'KEY1_woRBM rescue'});

%Y_ALL = cell(size(DataPaths));
Y_ALL = cell(1,3);
Y = [];
Yerr = [];
for i = 1:length(DataPaths)
    k = xIndex(i);
    load([DataPaths{i},'\PROPERTIES_SNAP.mat'])
    % Create a logical index for structures you want to keep
    indexToKeep = [PROPERTIES.number] ~= 0;
    
    % Use this index to create a new array
    PROPERTIES = PROPERTIES(indexToKeep);

    y = cat(1,PROPERTIES.SumIntSNAPCell);
    Y_ALL{k} = vertcat(Y_ALL{k}, y);
    Y = [Y,mean(y)];
    Yerr = [Yerr,std(y)];
end

labels = {'Wild Type','KEY1 Rescue','KEY1_woRBM'};
best = [1 2 3];

FUNC_PlotBoxWhisker(Title,best,Y_ALL)
FUNC_PlotViolin(Title, best, Y_ALL, labels);
box


%% Enrichment in Condensate

Title = 'Enrichment of Intensity in Condensate';

X = categorical({'Wild-type', 'KEY1 rescue','KEY1_woRBM rescue'});
X = reordercats(X,{'Wild-type','KEY1 rescue', 'KEY1_woRBM rescue'});

%Y_ALL = cell(size(DataPaths));
Y_ALL = cell(1,3);
Y = [];
Yerr = [];
for i = 1:length(DataPaths)
    k = xIndex(i);
    load([DataPaths{i},'\PROPERTIES_SNAP.mat']);
    load([DataPaths{i},'\CELL.mat']);
    load([DataPaths{i},'\PhysicalParameters.mat']);
    % Create a logical index for structures you want to keep
    indexToKeep = [PROPERTIES.number] ~= 0;
    
    % Use this index to create a new array
    PROPERTIES = PROPERTIES(indexToKeep);
    CELL = CELL(indexToKeep);
    y1 = [];
    a = [];
    for j = 1:size(PROPERTIES,2)
        y1 = [y1, sum(PROPERTIES(j).SumIntSNAPIn)];
    end
    y2 = reshape(cat(1,PROPERTIES.SumIntSNAPCell),1,[]);
    y = y1./y2;
    Y_ALL{k} = [Y_ALL{k}, y];
    Y = [Y,mean(y)];
    Yerr = [Yerr,std(y)];
end

labels = {'Wild Type','KEY1 Rescue','KEY1_woRBM'};
best = [1 2 3];

FUNC_PlotBoxWhisker(Title,best,Y_ALL)
FUNC_PlotViolin(Title, best, Y_ALL, labels);
box

%% Average Intensity Inside condensate

Title = 'Average Intensity In Condensate';

X = categorical({'Wild-type', 'KEY1 rescue','KEY1_woRBM rescue'});
X = reordercats(X,{'Wild-type','KEY1 rescue', 'KEY1_woRBM rescue'});

xIndex = [1 1 2 2 2 3 3];

%Y_ALL = cell(size(DataPaths));
Y_ALL = cell(1,3);
Y = [];
Yerr = [];
for i = 1:length(DataPaths)
    k = xIndex(i);
    load([DataPaths{i},'\PROPERTIES_SNAP.mat']);
    load([DataPaths{i},'\CELL.mat']);
    load([DataPaths{i},'\PhysicalParameters.mat']);
    % Create a logical index for structures you want to keep
    indexToKeep = [PROPERTIES.number] ~= 0;
    
    % Use this index to create a new array
    PROPERTIES = PROPERTIES(indexToKeep);
    CELL = CELL(indexToKeep);
    y1 = [];
    a = [];
    for j = 1:size(PROPERTIES,2)
        y1 = [y1, sum(PROPERTIES(j).MeanIntSNAPIn)];
    end
    y = y1;
    Y_ALL{k} = [Y_ALL{k}, y];
    Y = [Y,mean(y)];
    Yerr = [Yerr,std(y)];
end

labels = {'Wild Type','KEY1 Rescue','KEY1_woRBM'};
best = [1 2 3];

FUNC_PlotBoxWhisker(Title,best,Y_ALL)
FUNC_PlotViolin(Title, best, Y_ALL, labels);
box


%% Average Intensity In Cell

Title = 'Average Intensity';

X = categorical({'Wild-type', 'KEY1 rescue','KEY1_woRBM rescue'});
X = reordercats(X,{'Wild-type','KEY1 rescue', 'KEY1_woRBM rescue'});

%Y_ALL = cell(size(DataPaths));
Y_ALL = cell(1,3);
Y = [];
Yerr = [];
for i = 1:length(DataPaths)
    k = xIndex(i);
    load([DataPaths{i},'\PROPERTIES_SNAP.mat'])
    % Create a logical index for structures you want to keep
    indexToKeep = [PROPERTIES.number] ~= 0;
    
    % Use this index to create a new array
    PROPERTIES = PROPERTIES(indexToKeep);

    y = cat(1,PROPERTIES.MeanIntSNAPCell);
    Y_ALL{k} = vertcat(Y_ALL{k}, y);
    Y = [Y,mean(y)];
    Yerr = [Yerr,std(y)];
end

labels = {'Wild Type','KEY1 Rescue','KEY1_woRBM'};
best = [1 2 3];

FUNC_PlotBoxWhisker(Title,best,Y_ALL)
FUNC_PlotViolin(Title, best, Y_ALL, labels);
box
%%
dataset1 = Y_ALL{2};
dataset2 = Y_ALL{3};
[~,p_value] = ttest2(dataset1,dataset2)

%% Enrichment in Condensate from Average

Title = 'Enrichment of Intensity in Condensate SumIn/(Total-SumIn)';

X = categorical({'Wild-type', 'KEY1 rescue','KEY1_woRBM rescue'});
X = reordercats(X,{'Wild-type','KEY1 rescue', 'KEY1_woRBM rescue'});

%Y_ALL = cell(size(DataPaths));
Y_ALL = cell(1,3);
Y = [];
Yerr = [];
for i = 1:length(DataPaths)
    k = xIndex(i);
    load([DataPaths{i},'\PROPERTIES_SNAP.mat']);
    load([DataPaths{i},'\CELL.mat']);
    load([DataPaths{i},'\PhysicalParameters.mat']);
    % Create a logical index for structures you want to keep
    indexToKeep = [PROPERTIES.number] ~= 0;
    
    % Use this index to create a new array
    PROPERTIES = PROPERTIES(indexToKeep);
    CELL = CELL(indexToKeep);
    y1 = [];
    a = [];
    for j = 1:size(PROPERTIES,2)
        y1 = [y1, sum(PROPERTIES(j).SumIntSNAPIn)];
    end
    y2 = reshape(cat(1,PROPERTIES.SumIntSNAPCell),1,[]);
    y = y1./(y2-y1);
    Y_ALL{k} = [Y_ALL{k}, y];
    Y = [Y,mean(y)];
    Yerr = [Yerr,std(y)];
end

labels = {'Wild Type','KEY1 Rescue','KEY1_woRBM'};
best = [1 2 3];

FUNC_PlotBoxWhisker(Title,best,Y_ALL)
FUNC_PlotViolin(Title, best, Y_ALL, labels);
box

%% Number of condensates

Title = 'Number of pyrenoids per cell';

X = categorical({'Wild-type', 'KEY1 rescue','KEY1_woRBM rescue','key1 mutant'});
X = reordercats(X,{'Wild-type','KEY1 rescue', 'KEY1_woRBM rescue','key1 mutant'});

%Y_ALL = cell(size(DataPaths));
Y_ALL = cell(1,4);
Y = [];
Yerr = [];
for i = 1:length(DataPaths)
    k = xIndex(i);
    load([DataPaths{i},'\PROPERTIES_SNAP.mat'])
    % Create a logical index for structures you want to keep
    indexToKeep = [PROPERTIES.number] ~= 0;
    
    % Use this index to create a new array
    PROPERTIES = PROPERTIES(indexToKeep);

    y = cat(1,PROPERTIES.number);
    Y_ALL{k} = vertcat(Y_ALL{k}, y);
    Y = [Y,mean(y)];
    Yerr = [Yerr,std(y)];
end

% Add in key1 mutant for pyrenoid number count

DataPath = 'E:\PrincetonData\20231103';

DataSubPath_mutant1 = [DataPath, '\key1\R0C0'];

DataSubPath_mutant2 = [DataPath, '\key1\R0C1'];

DataPaths2 = {DataSubPath_mutant1, DataSubPath_mutant2};

for i = 1:length(DataPaths2)
    k = 4;
    load([DataPaths2{i},'\PROPERTIES.mat'])
    % Create a logical index for structures you want to keep
    indexToKeep = [PROPERTIES.number] ~= 0;
    
    % Use this index to create a new array
    PROPERTIES = PROPERTIES(indexToKeep);

    y = cat(1,PROPERTIES.number);
    Y_ALL{k} = vertcat(Y_ALL{k}, y);
    Y = [Y,mean(y)];
    Yerr = [Yerr,std(y)];
end


labels = {'Wild Type','KEY1 Rescue','KEY1_woRBM','key1 mutant'};
best = [1 2 3 4];

FUNC_PlotBoxWhisker(Title,best,Y_ALL)
FUNC_PlotViolin(Title, best, Y_ALL, labels);
box

