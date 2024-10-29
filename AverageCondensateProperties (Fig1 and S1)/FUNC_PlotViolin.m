function [V] = FUNC_PlotViolin(Title,best,Y_ALL,labels)

data = []; % This will hold the corresponding data


group=[];

for i = 1:length(best)
    j = best(i);
    y = reshape(Y_ALL{j},1,[]);
    data = [data, y]; % Append the data
    group = [group, repmat(labels(j),1, length(Y_ALL{j}))]; % Append the labels
%     if i<length(best)
%         temp = Y_ALL{best(i):best(i+1)}
%         group = [group, repmat(labels(j),1,length(temp))]
%     else
%         temp = Y_ALL{best(end):length(labels)};
%         group = [group, repmat(labels(j),1,length(temp))];
%     end
end

figure3 = figure('color', [1 1 1]);
V = violinplot(data,group,'GroupOrder',labels(best),'ViolinColor',[0 0 1]);
hold on
title(Title);
box
end


