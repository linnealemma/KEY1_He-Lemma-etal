% Script to Plot average condensate properties for key1 mutant, KEY1 rescue
% and wild-type (from 2D analysis)

%% RBCS1-Venus

DataPath = 'E:\PrincetonData\20231103';

DataSubPath_mutant1 = [DataPath, '\key1\R0C0'];

DataSubPath_mutant2 = [DataPath, '\key1\R0C1'];


DataSubPath_rescue1 = [DataPath, '\KEY1-rescue\R0C0'];


DataSubPath_rescue2 = [DataPath, '\KEY1-rescue\R1C1'];

DataSubPath_wt1 = [DataPath, '\wild-type']; %bin 1, all the rest of the data is bin 2, so we do not include this in the plots

DataSubPath_wt2 = [DataPath, '\wild-type\bin2\R0C0'];

DataSubPath_wt3 = [DataPath, '\wild-type\bin2\R0C1'];

%% Figure 1I (2024-10-18) Number of pyrenoids/cell for wild type, key1-1 mutant and key1-1;KEY1-SNAP rescue with RBCS1-Venus labeled

DataPaths = {DataSubPath_wt1, DataSubPath_wt2, DataSubPath_wt3, DataSubPath_mutant1, DataSubPath_mutant2, DataSubPath_rescue1, DataSubPath_rescue2};

X = categorical({'Wild Type1','Wild Type2','Wild Type3','key1 mutant1','key1 mutant2','KEY1 rescue1','KEY1 rescue2'});
X = reordercats(X,{'Wild Type1','Wild Type2','Wild Type3','key1 mutant1','key1 mutant2','KEY1 rescue1','KEY1 rescue2'});

Y_ALL = cell(size(DataPaths));
Y = [];
Yerr = [];
for i = 1:length(DataPaths)
    load([DataPaths{i},'\PROPERTIES.mat'])
    y = cat(1,PROPERTIES.number);
    Y_ALL{i} = y;
    Y = [Y,mean(y)];
    Yerr = [Yerr,std(y)];
end

figure1 = figure('color',[1 1 1]);

Y_ALL_compact = Y_ALL;

Y_ALL_compact{1} = [Y_ALL{2}; Y_ALL{3}]
Y_ALL_compact{2} = [Y_ALL{4}; Y_ALL{5}];
Y_ALL_compact{3} = [Y_ALL{6}; Y_ALL{7}];


Title = 'Number of pyrenoids per cell';
group = []; % This will hold the group labels
data = []; % This will hold the corresponding data
labels = {'Wild Type1', 'WildType2','Wild Type3','key1 mutant1', 'key1 mutant2', 'KEY1 rescue1', 'KEY1 rescue2'};
best = [1 2 3 4 5 6 7];
% labels = {'Wild Type', 'key1 mutant','KEY1 rescue'};
% best = [1 2 3]
%g1 = repmat(labels(best(1)),size(Y_ALL{1},2),1);
X=[];

Y_ALL=Y_ALL_compact;

for i = 1:length(best)
    j = best(i);
    y = reshape(Y_ALL{j},1,[]);
    data = [data, y]; % Append the data
    X = [X,  repmat(i, 1, length(Y_ALL{j}))];
    group = [group, repmat(labels(j),1, length(Y_ALL{j}))]; % Append the labels
end

Y_ALL{1}(Y_ALL{1}>3)=[];

FUNC_PlotBoxWhisker(Title,best,Y_ALL)
FUNC_PlotViolin(Title, best, Y_ALL, labels);
box


labels = {'Wild Type', 'key1 mutant','KEY1 rescue'};
best = [1 2 3]

FUNC_PlotBoxWhisker(Title,best,Y_ALL_compact)
FUNC_PlotViolin(Title, best, Y_ALL_compact, labels);
box

[~,p_value] = ttest2(Y_ALL_compact{1},Y_ALL_compact{2})
[~,p_value] = ttest2(Y_ALL_compact{1},Y_ALL_compact{3})
[~,p_value] = ttest2(Y_ALL_compact{2},Y_ALL_compact{3})
p = ranksum(Y_ALL_compact{1},Y_ALL_compact{2})
p = ranksum(Y_ALL_compact{1},Y_ALL_compact{3})
p = ranksum(Y_ALL_compact{2},Y_ALL_compact{3})

%% Percent Area of the Cell that is the largest condensate
% From area
Title = 'Percent Area Condensate';

X = categorical({'Wild Type1','Wild Type2','Wild Type3','key1 mutant1','key1 mutant2','KEY1 rescue1','KEY1 rescue2'});
X = reordercats(X,{'Wild Type1','Wild Type2','Wild Type3','key1 mutant1','key1 mutant2','KEY1 rescue1','KEY1 rescue2'});

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
        if PROPERTIES(j).number>0
            y1 = [y1, max(PROPERTIES(j).area)];
        else
            y1 = [y1, 0];
        end
        a = [a, (CELL(j).area).*umperpix_x*umperpix_y];
    end
    

    y = y1./a*100;
%     y1(y1>21) = [];
    Y_ALL{i} = y;
    Y = [Y,mean(y)];
    Yerr = [Yerr,std(y)];
end

labels = {'Wild Type1', 'WildType2','Wild Type3','key1 mutant1', 'key1 mutant2', 'KEY1 rescue1', 'KEY1 rescue2'};
best = [1 2 3 4 5 6 7];

FUNC_PlotBoxWhisker(Title,best,Y_ALL)
FUNC_PlotViolin(Title, best, Y_ALL, labels);
box


Y_ALL_compact = Y_ALL;

Y_ALL_compact{1} = [Y_ALL{2}, Y_ALL{3}]
Y_ALL_compact{2} = [Y_ALL{4}, Y_ALL{5}];
Y_ALL_compact{3} = [Y_ALL{6}, Y_ALL{7}];

% Y_ALL_compact{2} = [Y_ALL{2}; Y_ALL{3}];
% Y_ALL_compact{3} = [Y_ALL{4}; Y_ALL{5}];

labels = {'Wild Type', 'key1 mutant','KEY1 rescue'};
best = [1 2 3]

FUNC_PlotBoxWhisker(Title,best,Y_ALL_compact)
FUNC_PlotViolin(Title, best, Y_ALL_compact, labels);
box

[~,p_value] = ttest2(Y_ALL_compact{1},Y_ALL_compact{2})
[~,p_value] = ttest2(Y_ALL_compact{1},Y_ALL_compact{3})
[~,p_value] = ttest2(Y_ALL_compact{2},Y_ALL_compact{3})
p = ranksum(Y_ALL_compact{1},Y_ALL_compact{2})
p = ranksum(Y_ALL_compact{1},Y_ALL_compact{3})
p = ranksum(Y_ALL_compact{2},Y_ALL_compact{3})