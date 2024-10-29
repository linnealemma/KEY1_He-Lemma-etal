% Script to Plot average condensate properties for key1 mutant, KEY1 rescue
% and wild-type (from 2D analysis)

%% RBCS1-Venus

DataPath = 'E:\PrincetonData\20231103';

DataSubPath_mutant1 = [DataPath, '\key1\R0C0'];

DataSubPath_mutant2 = [DataPath, '\key1\R0C1'];


DataSubPath_rescue1 = [DataPath, '\KEY1-rescue\R0C0'];


DataSubPath_rescue2 = [DataPath, '\KEY1-rescue\R1C1'];

DataSubPath_wt1 = [DataPath, '\wild-type']; %bin 1, all the rest of the data is bin 2

DataSubPath_wt2 = [DataPath, '\wild-type\bin2\R0C0'];

DataSubPath_wt3 = [DataPath, '\wild-type\bin2\R0C1'];


%%

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
%% Radius of All Condensates
Title = 'Average Radius of All Condensates';

X = categorical({'Wild Type1','Wild Type2','Wild Type3','key1 mutant1','key1 mutant2','KEY1 rescue1','KEY1 rescue2'});
X = reordercats(X,{'Wild Type1','Wild Type2','Wild Type3','key1 mutant1','key1 mutant2','KEY1 rescue1','KEY1 rescue2'});

Y_ALL = cell(size(DataPaths));
Y = [];
Yerr = [];
for i = 1:length(DataPaths)
    load([DataPaths{i},'\PROPERTIES.mat'])
    % Create a logical index for structures you want to keep
    indexToKeep = [PROPERTIES.number] ~= 0;
    
    % Use this index to create a new array
    PROPERTIES = PROPERTIES(indexToKeep);

    y = cat(1,PROPERTIES.radius);
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

Y_ALL_compact{1} = [Y_ALL{2}; Y_ALL{3}]
Y_ALL_compact{2} = [Y_ALL{4}; Y_ALL{5}];
Y_ALL_compact{3} = [Y_ALL{6}; Y_ALL{7}];

% Y_ALL_compact{2} = [Y_ALL{2}; Y_ALL{3}];
% Y_ALL_compact{3} = [Y_ALL{4}; Y_ALL{5}];

labels = {'Wild Type', 'key1 mutant','KEY1 rescue'};
best = [1 2 3]

FUNC_PlotBoxWhisker(Title,best,Y_ALL_compact)
FUNC_PlotViolin(Title, best, Y_ALL_compact, labels);
box

%% Radius of Largest

Title = 'Radius of Largest';

X = categorical({'Wild Type1','Wild Type2','Wild Type3','key1 mutant1','key1 mutant2','KEY1 rescue1','KEY1 rescue2'});
X = reordercats(X,{'Wild Type1','Wild Type2','Wild Type3','key1 mutant1','key1 mutant2','KEY1 rescue1','KEY1 rescue2'});
Y_ALL = cell(size(DataPaths));
Y = [];
Yerr = [];
for i = 1:length(DataPaths)
    load([DataPaths{i},'\PROPERTIES.mat']);
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
        y1 = [y1, max(PROPERTIES(j).radius)];
    end
    y = y1;
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

Y_ALL_compact{1} = [Y_ALL{1}, Y_ALL{2}, Y_ALL{3}]
Y_ALL_compact{2} = [Y_ALL{4}, Y_ALL{5}];
Y_ALL_compact{3} = [Y_ALL{6}, Y_ALL{7}];

% Y_ALL_compact{2} = [Y_ALL{2}; Y_ALL{3}];
% Y_ALL_compact{3} = [Y_ALL{4}; Y_ALL{5}];

labels = {'Wild Type', 'key1 mutant','KEY1 rescue'};
best = [1 2 3]

FUNC_PlotBoxWhisker(Title,best,Y_ALL_compact)
FUNC_PlotViolin(Title, best, Y_ALL_compact, labels);
box

%% Volume of Largest
% From the equivalent diameter
Title = 'Volume of Largest';

X = categorical({'Wild Type','key1 mutant1','key1 mutant2','KEY1 rescue1','KEY1 rescue2'});
X = reordercats(X,{'Wild Type','key1 mutant1','key1 mutant2','KEY1 rescue1','KEY1 rescue2'});

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
        y1 = [y1, max(PROPERTIES(j).radius)];
    end
    y = 4/3*pi*y1.^3;
    Y_ALL{i} = y;
    Y = [Y,mean(y)];
    Yerr = [Yerr,std(y)];
end

FUNC_PlotBoxWhisker(Title,best,Y_ALL)
FUNC_PlotViolin(Title, best, Y_ALL, labels);
box


%% Area of Largest
% From area
Title = 'Area of Largest';

X = categorical({'Wild Type','key1 mutant1','key1 mutant2','KEY1 rescue1','KEY1 rescue2'});
X = reordercats(X,{'Wild Type','key1 mutant1','key1 mutant2','KEY1 rescue1','KEY1 rescue2'});

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
        y1 = [y1, max(PROPERTIES(j).area)];
    end
    y = y1;
    y1(y1>21) = [];
    Y_ALL{i} = y1;
    Y = [Y,mean(y)];
    Yerr = [Yerr,std(y)];
end

FUNC_PlotBoxWhisker(Title,best,Y_ALL)
FUNC_PlotViolin(Title, best, Y_ALL, labels);
box

%% Percent Area of the Cell that is condensate
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

%% Percent Area of the Cell that is condensate
% From radius
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
            y1 = [y1, max(PROPERTIES(j).radius)];
        else
            y1 = [y1, 0];
        end
        a = [a, (CELL(j).area).*umperpix_x*umperpix_y];
    end
    

    y = (pi*y1.^2)./a*100;
    y(y>20)=[];
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

%% Volume of Largest Normalized by Cell Volume

% Would like to normalize this by cell size, but will have to store cell
% size in PROPERTIES.
% Oh, maybe not. It should be in CELL.mat.

Title = 'Volume of Largest Normalized by Cell Volume';


X = categorical({'Wild Type','key1 mutant1','key1 mutant2','KEY1 rescue1','KEY1 rescue2'});
X = reordercats(X,{'Wild Type','key1 mutant1','key1 mutant2','KEY1 rescue1','KEY1 rescue2'});

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
        a = [a, (CELL(j).volume*voxel)];
    end
   
    y = y1./a;
    y(y>(21/mean(a))) = [];
    Y_ALL{i} = y;
end

FUNC_PlotBoxWhisker(Title,best,Y_ALL)
FUNC_PlotViolin(Title, best, Y_ALL, labels);
box

FUNC_PlotBoxWhisker(Title,[1 2 3 4 5 6 7],Y_ALL)
FUNC_PlotViolin(Title, [1 2 3 4 5 6 7], Y_ALL, labels);
box

%% Radius of Largest Normalized by Cell Volume

% Would like to normalize this by cell size, but will have to store cell
% size in PROPERTIES.
% Oh, maybe not. It should be in CELL.mat.

Title = 'Radius of Largest Normalized by Cell Volume';

X = categorical({'Wild Type','key1 mutant1','key1 mutant2','KEY1 rescue1','KEY1 rescue2'});
X = reordercats(X,{'Wild Type','key1 mutant1','key1 mutant2','KEY1 rescue1','KEY1 rescue2'});

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
        y1 = [y1, max(PROPERTIES(j).radius)];
        a = [a, (CELL(j).volume)*voxel];
    end
    y = y1/mean(a);
    Y_ALL{i} = y;
    Y = [Y,mean(y)];
    Yerr = [Yerr,std(y)];
end

FUNC_PlotBoxWhisker(Title,best,Y_ALL)
FUNC_PlotViolin(Title, best, Y_ALL, labels);
box

FUNC_PlotBoxWhisker(Title,[1 2 3 4 5 6 7],Y_ALL)
FUNC_PlotViolin(Title, [1 2 3 4 5 6 7], Y_ALL, labels);
box

%% Volume of All Condensates 

% Would like to normalize this by cell size, but will have to store cell
% size in PROPERTIES.
% Oh, maybe not. It should be in CELL.mat.
Title = 'Volume of All Condensates';

X = categorical({'Wild Type','key1 mutant1','key1 mutant2','KEY1 rescue1','KEY1 rescue2'});
X = reordercats(X,{'Wild Type','key1 mutant1','key1 mutant2','KEY1 rescue1','KEY1 rescue2'});

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

FUNC_PlotBoxWhisker(Title,best,Y_ALL)
FUNC_PlotViolin(Title, best, Y_ALL, labels);
box

FUNC_PlotBoxWhisker(Title,[1 2 3 4 5 6 7],Y_ALL)
FUNC_PlotViolin(Title, [1 2 3 4 5 6 7], Y_ALL, labels);
box

%% Volume of Largest Normalized by Cell Volume

% Would like to normalize this by cell size, but will have to store cell
% size in PROPERTIES.
% Oh, maybe not. It should be in CELL.mat.
Title = 'Volume of Largest'

X = categorical({'Wild Type','key1 mutant1','key1 mutant2','KEY1 rescue1','KEY1 rescue2'});
X = reordercats(X,{'Wild Type','key1 mutant1','key1 mutant2','KEY1 rescue1','KEY1 rescue2'});

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
        y1 = [y1, max(PROPERTIES(j).volume)];
        a = [a, (CELL(j).bb(4)*CELL(j).bb(5)*CELL(j).bb(6))*voxel];
    end
    y = y1;
    Y_ALL{i} = y;
    Y = [Y,mean(y)];%./mean(a);
    Yerr = [Yerr,std(y)];%./mean(a);
end

FUNC_PlotBoxWhisker(Title,best,Y_ALL)
FUNC_PlotViolin(Title, best, Y_ALL, labels);
box

% ylim([-0.5 6])

%% Cell Area by cell type
% This is using the pixels that are positive so maybe should use bb instead
% because pixels seems sensitive to thresholding
Title = 'Cell Area';

X = categorical({'Wild Type','key1 mutant1','key1 mutant2','KEY1 rescue1','KEY1 rescue2'});
X = reordercats(X,{'Wild Type','key1 mutant1','key1 mutant2','KEY1 rescue1','KEY1 rescue2'});

Y_ALL = cell(size(DataPaths));
Y = [];
Yerr = [];
for i = 1:length(DataPaths)
    load([DataPaths{i},'\CELL.mat']);
    load([DataPaths{i},'\PhysicalParameters.mat']);
    y = cat(1,CELL.area);
    Y_ALL{i} = y.*umperpix_x*umperpix_y;
    Y = [Y,mean(y)];
    Yerr = [Yerr,std(y)];
end

labels = {'Wild Type', 'key1 mutant1', 'key1 mutant2', 'KEY1 rescue1', 'KEY1 rescue2'};
best = [1 2 3 4 5]

FUNC_PlotBoxWhisker(Title,best,Y_ALL)
FUNC_PlotViolin(Title, best, Y_ALL, labels);
box

Y_ALL_compact = Y_ALL;

Y_ALL_compact{2} = [Y_ALL{2}; Y_ALL{3}];
Y_ALL_compact{3} = [Y_ALL{4}; Y_ALL{5}];

labels = {'Wild Type', 'key1 mutant','KEY1 rescue'};
best = [1 2 3]

FUNC_PlotBoxWhisker(Title,best,Y_ALL_compact)
FUNC_PlotViolin(Title, best, Y_ALL_compact, labels);
box



%% Cell Volume by cell type bounding box
% Here we will use bb instead 
Title = 'Cell Volume (bb)';

DataPaths = {DataSubPath_wt1,DataSubPath_mutant1,DataSubPath_mutant2,DataSubPath_mutant3,DataSubPath_rescue1,DataSubPath_rescue2, DataSubPath_rescue3};

X = categorical({'Wild Type','key1 mutant1','key1 mutant2','key1 mutant3','KEY1 rescue1','KEY1 rescue2','KEY1 rescue3'});
X = reordercats(X,{'Wild Type','key1 mutant1','key1 mutant2','key1 mutant3','KEY1 rescue1','KEY1 rescue2','KEY1 rescue3'});

Y_ALL = cell(size(DataPaths));
Y = [];
Yerr = [];
for i = 1:length(DataPaths)
    load([DataPaths{i},'\CELL.mat'])
    y = [];
    for j = 1:size(CELL,2)
        y = [y, CELL(j).bb(4)*CELL(j).bb(5)*CELL(j).bb(6)*voxel];
    end
    %y = cat(1,CELL.bb);
    Y_ALL{i} = y;
    Y = [Y,mean(y)];
    Yerr = [Yerr,std(y)];
end


FUNC_PlotBoxWhisker(Title,best,Y_ALL)
FUNC_PlotViolin(Title, best, Y_ALL, labels);
box


FUNC_PlotBoxWhisker(Title,[1 2 3 4 5 6 7],Y_ALL)
FUNC_PlotViolin(Title, [1 2 3 4 5 6 7], Y_ALL, labels);
box



%% Total Volume Normalized by Cell Volume

Title = 'Total Volume Normalized by Cell Volume'

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
        a = [a, (CELL(j).volume)*voxel];
    end
    y = y1;
    Y_ALL{i} = y./mean(a);
    Y = [Y,mean(y)]./mean(a);
    Yerr = [Yerr,std(y)]./mean(a);
end

FUNC_PlotBoxWhisker(Title,best,Y_ALL)
FUNC_PlotViolin(Title, best, Y_ALL, labels);
box


FUNC_PlotBoxWhisker(Title,[1 2 3 4 5 6 7],Y_ALL)
FUNC_PlotViolin(Title, [1 2 3 4 5 6 7], Y_ALL, labels);
box

%% Total Intensity in Condensate

Title = 'Total Intensity Inside Condensate';

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
    a = [];
    y1 = [];
    for j = 1:size(PROPERTIES,2)
        y1 = [y1, sum(PROPERTIES(j).SumIntIn)];
        a = [a, (CELL(j).volume)*voxel];
    end
    y = y1;
    Y_ALL{i} = y;
    Y = [Y,mean(y)];
    Yerr = [Yerr,std(y)];
end

FUNC_PlotBoxWhisker(Title,best,Y_ALL)
FUNC_PlotViolin(Title, best, Y_ALL, labels);
box


FUNC_PlotBoxWhisker(Title,[1 2 3 4 5 6 7],Y_ALL)
FUNC_PlotViolin(Title, [1 2 3 4 5 6 7], Y_ALL, labels);
box
%%


%% Total Intensity in Condensate Normalized by Cell Volume

Title = 'Total Intensity Inside Condensate Normalized by Cell Volume';

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
    a = [];
    y1 = [];
    for j = 1:size(PROPERTIES,2)
        y1 = [y1, sum(PROPERTIES(j).SumIntIn)];
        a = [a, (CELL(j).volume)*voxel];
    end
    y = y1;
    Y_ALL{i} = y./mean(a);
    Y = [Y,mean(y)]./mean(a);
    Yerr = [Yerr,std(y)]./mean(a);
end

FUNC_PlotBoxWhisker(Title,best,Y_ALL)
FUNC_PlotViolin(Title, best, Y_ALL, labels);
box


FUNC_PlotBoxWhisker(Title,[1 2 3 4 5 6 7],Y_ALL)
FUNC_PlotViolin(Title, [1 2 3 4 5 6 7], Y_ALL, labels);
box


%% Total Intensity in Condensate/Total Intensity out of condensate

Title = 'Total Intensity in Condensate/Total Intensity out of Condensate';

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
    a = [];
    y1 = [];
    for j = 1:size(PROPERTIES,2)
        y1 = [y1, sum(PROPERTIES(j).SumIntIn)./sum(PROPERTIES(j).SumIntOut)];
    end
    y = y1;
    Y_ALL{i} = y;
    Y = [Y,mean(y)];
    Yerr = [Yerr,std(y)];
end

FUNC_PlotBoxWhisker(Title,best,Y_ALL)
FUNC_PlotViolin(Title, best, Y_ALL, labels);
box


FUNC_PlotBoxWhisker(Title,[1 2 3 4 5 6 7],Y_ALL)
FUNC_PlotViolin(Title, [1 2 3 4 5 6 7], Y_ALL, labels);
box



%% Total Intensity in Condensate/Total Intensity in cell

Title = 'Total Intensity in Condensate/Total Intensity';

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
    a = [];
    y1 = [];
    for j = 1:size(PROPERTIES,2)
        y1 = [y1, sum(PROPERTIES(j).SumIntIn)./(sum(PROPERTIES(j).SumIntOut)+sum(PROPERTIES(j).SumIntIn))];
    end
    y = y1;
    Y_ALL{i} = y;
    Y = [Y,mean(y)];
    Yerr = [Yerr,std(y)];
end

FUNC_PlotBoxWhisker(Title,best,Y_ALL)
FUNC_PlotViolin(Title, best, Y_ALL, labels);
box


FUNC_PlotBoxWhisker(Title,[1 2 3 4 5 6 7],Y_ALL)
FUNC_PlotViolin(Title, [1 2 3 4 5 6 7], Y_ALL, labels);
box


%% Partition Coefficient

Title = 'Partition Coefficient (Mean Intensity In Condensate/Mean Intensity Out Condensate)';

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

    a = [];
    y1 = [];
    for j = 1:size(PROPERTIES,2)
        y1 = [y1, mean(PROPERTIES(j).MeanIntIn)/mean(PROPERTIES(j).MeanIntOut)];
    end
    y = y1;
    Y_ALL{i} = y;
    Y = [Y,mean(y)];
    Yerr = [Yerr,std(y)];
end

FUNC_PlotBoxWhisker(Title,best,Y_ALL)
FUNC_PlotViolin(Title, best, Y_ALL, labels);
box


FUNC_PlotBoxWhisker(Title,[1 2 3 4 5 6 7],Y_ALL)
FUNC_PlotViolin(Title, [1 2 3 4 5 6 7], Y_ALL, labels);
box

%% Average Intensity Inside the Condensate

Title = 'Average Intensity Inside the Condensates';

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
    a = [];
    y1 = [];
    for j = 1:size(PROPERTIES,2)
        y1 = [y1, mean(PROPERTIES(j).MeanIntIn)];
    end
    y = y1;
    Y_ALL{i} = y;
    Y = [Y,mean(y)];
    Yerr = [Yerr,std(y)];
end
FUNC_PlotBoxWhisker(Title,best,Y_ALL)
FUNC_PlotViolin(Title, best, Y_ALL, labels);
box


FUNC_PlotBoxWhisker(Title,[1 2 3 4 5 6 7],Y_ALL)
FUNC_PlotViolin(Title, [1 2 3 4 5 6 7], Y_ALL, labels);
box

% ylim([-0.5 6])

%% Plotting a bar graph old code

X = categorical({'Wild Type','key1 mutant1','key1 mutant2','key1 mutant3','KEY1 rescue1','KEY1 rescue2','KEY1 rescue3'});
X = reordercats(X,{'Wild Type','key1 mutant1','key1 mutant2','key1 mutant3','KEY1 rescue1','KEY1 rescue2','KEY1 rescue3'});
bar(X,Y);
hold on
er = errorbar(Y,Yerr);
er.Color = [0 0 0 ];
er.LineStyle = 'none';
title('Number of Pyrenoids Per Cell')

figure2 = figure('color',[1 1 1]);

best = [1,4,5];
% bad = [2,3,6,7];
% X2 = X;
% X2(bad) = [];
X2 = categorical({'Wild Type','key1 mutant','KEY1 rescue',});
X2 = reordercats(X2,{'Wild Type','key1 mutant','KEY1 rescue'});
%X2 = X(best);
%X2 = reordercats(X2,{X2(1),X2(2),X2(3)});
h = bar(X2,Y(best));
hold on
er = errorbar(Y(best),Yerr(best));
er.Color = [0 0 0 ];
er.LineStyle = 'none';

figure1 = figure('color',[1 1 1]);

X = categorical({'Wild Type','key1 mutant1','key1 mutant2','key1 mutant3','KEY1 rescue1','KEY1 rescue2','KEY1 rescue3'});
X = reordercats(X,{'Wild Type','key1 mutant1','key1 mutant2','key1 mutant3','KEY1 rescue1','KEY1 rescue2','KEY1 rescue3'});
bar(X,Y);
hold on
er = errorbar(Y,Yerr);
er.Color = [0 0 0 ];
er.LineStyle = 'none';
title(plotTitle)

figure2 = figure('color',[1 1 1]);

best = [1,4,5];
% bad = [2,3,6,7];
% X2 = X;
% X2(bad) = [];
X2 = categorical({'Wild Type','key1 mutant','KEY1 rescue',});
X2 = reordercats(X2,{'Wild Type','key1 mutant','KEY1 rescue'});
%X2 = X(best);
%X2 = reordercats(X2,{X2(1),X2(2),X2(3)});
h = bar(X2,Y(best));
hold on
er = errorbar(Y(best),Yerr(best));
er.Color = [0 0 0 ];
er.LineStyle = 'none';

title(plotTitle)

for i = 1:length(best)
    j = best(i);
    y = Y_ALL{j};
    x = ones(size(y)).*h.XEndPoints(i);
    swarmchart(x,y,'MarkerFaceColor',[0 0 0],'MarkerFaceAlpha',0.4,'MarkerEdgeColor',[0 0 0], ...
        'MarkerEdgeAlpha',0.4)
    hold on
end

