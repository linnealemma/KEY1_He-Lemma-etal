function FUNC_PlotBoxWhisker(Title,best,Y_ALL)

data = []; % This will hold the corresponding data

%g1 = repmat(labels(best(1)),size(Y_ALL{1},2),1);
X=[];

for i = 1:length(best)
    j = best(i);
    y = reshape(Y_ALL{j},1,[]);
    data = [data, y]; % Append the data
    X = [X,  repmat(i, 1, length(Y_ALL{j}))];
end

figure3 = figure('color', [1 1 1]);
boxchart(X,data,'Notch','off');
hold on
title(Title);
box
end