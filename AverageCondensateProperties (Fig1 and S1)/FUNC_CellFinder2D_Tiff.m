function [CELL] = FUNC_CellFinder2D_Tiff(ChlorFile, z_size, chlor_background)

%% Load data (chlorophyll channel), threshold appropriately, get bounding boxes
% script to plot the bounding boxes of each cell 

[DATA] = tiffreadVolume(ChlorFile);
%%
DATA_2c = DATA-chlor_background;

DATA_mip = max(DATA_2c,[],3);

imcontrastscale(1)=min(DATA_mip,[],'all');
imcontrastscale(2)=mean(DATA_mip,'all')*10;

figure0 = figure('color',[1 1 1]);
imshow(DATA_mip,imcontrastscale);

sigma = [1 1]; % width of gaussianlc

DATA_g = imgaussfilt(DATA_mip,sigma);

DATA_3c = uint16(DATA_g);
%Tc = graythresh(DATA_3c(:,:,:));
Tc = 0.05;
%T = graythresh(DATA_3c);
DATA_BW_Temp = imbinarize(DATA_3c,Tc);

figure2 = figure('color',[1 1 1]);
imshow(DATA_BW_Temp);

 %%
% se = strel('disk',5);
% erodedI = imerode(DATA_BW_Temp,se);
% imshow(erodedI)
% % figure2 = figure('color',[1 1 1]);
% % imshow(DATA_BW_Temp(:,:,12))
% CC = bwconncomp(erodedI);
% title('eroded')
%%
CC = bwconncomp(DATA_BW_Temp);
    S = regionprops(CC,'centroid','area','boundingbox');
    
    cent=cat(1,S.Centroid);
    area=cat(1,S.Area);
    bb=cat(1,S.BoundingBox);
    
    data = [cent(:,1),cent(:,2),area,bb];
    data(data(:,3)<1000,:)=[];
    
    bb = data(:,4:7);

    %plot all the cells and bounding boxes as a sanity check
    % for j = 1:z_size
    colors = parula(length(data));
    
    figure2 = figure('color',[1 1 1]);
    
    imcontrastscale(1)=min(DATA_mip,[],'all');
    imcontrastscale(2)=mean(DATA_mip,'all')*5;
    flipud(imshow(DATA_mip,imcontrastscale));
    %pad = 50;
    pad = 0;
    %flipud(imshow(DATA(:,:,j,1),imcontrastscale));
    hold on
    for i=1:length(bb)
        %scatter(DATA_PLOT(i,1),DATA_PLOT(i,2),50,[0 0 1],'filled');
        rectangle('Position',[bb(i,1)-pad, bb(i,2)-pad, bb(i,3)+1.5*pad, bb(i,4)+1.5*pad],'edgecolor',colors(i,:));
    end
    %end

%% save the data for each cell in a structure
%     if j>1
%         c = length(CELL);
%     else
%         c=0;
%     end
c=0;
    for i = 1:length(data)
        CELL(i+c).position = data(i,1:2);
        CELL(i+c).area = data(i,3);
        CELL(i+c).frame = j;
        CELL(i+c).bb = data(i,4:7);
    end

% disp(size(CELL,2))
%end

%%

% Get rid of cells that are way too big or way too small (based on bounding
% box area).
% maximumarea = 5*10^5;
% minarea = 3*10^2;
% areas = [CELL(1:end).area];
% badarea = find(area>maximumarea);
% badarea = [badarea, find(areas<minarea)];
% CELL(badarea)=[];

% disp(size(CELL,2))

% %To crop image
for i =1:length(CELL)
    x1=ceil(CELL(i).bb(1));
    y1=ceil(CELL(i).bb(2));
   % z1=ceil(CELL(i).bb(3));
    x2=round(x1+CELL(i).bb(3));
    y2=round(y1+CELL(i).bb(4));
   % z2=round(y1+CELL(i).bb(6));
    if x2>size(DATA,2)
        x2=size(DATA,2);
    end
    if y2 > size(DATA,1)
        y2 = size(DATA,1);
    end
    test=DATA(y1:y2,x1:x2);
    figure1 = figure('color',[1 1 1]);
    imshow(test,imcontrastscale);
    title(['i=',num2str(i)])
end
disp('Delete badly identified cells from CELL')

pause
prompt = "Make a list of the indices of badly identified cells";
badcells = input(prompt);
%badcells = [ 1, 4, 5, 8, 7, 9, 19, 23];
CELL(badcells)=[];

% disp(size(CELL,2))

% Replot all "good" cells over original image -- if there are still badly
% identified ones, re-run and prune more aggressively
imcontrastscale(1)=min(DATA,[],'all');
imcontrastscale(2)=mean(DATA,'all')*5;
flipud(imshow(DATA,imcontrastscale));
%flipud(imshow(DATA(:,:,j,1),imcontrastscale));
hold on
for i=1:length(CELL)
    %scatter(DATA_PLOT(i,1),DATA_PLOT(i,2),50,[0 0 1],'filled');
    rectangle('Position',[CELL(i).bb(1), CELL(i).bb(2), CELL(i).bb(3), CELL(i).bb(4)],'edgecolor',colors(i,:));
end

%%
SavePath = fileparts(ChlorFile);
save([SavePath,'\CELL.mat'],'CELL')
