%% Find cells and measure the intensity of KEY1-SNAP signal (no time series)
% Data is in volumetric Tiff files

%% Data that has been analyzed by this script
DataPath = 'E:\PrincetonData\20231103';
DataSubPath = [DataPath, '\wild-type\bin2\R0C0'];

DataPath = 'E:\PrincetonData\20231103';
DataSubPath = [DataPath, '\wild-type\bin2\R0C1'];

DataPath = 'E:\PrincetonData\20231103';
DataSubPath = [DataPath, '\KEY1-rescue\R0C0'];

DataPath = 'E:\PrincetonData\20231103';
DataSubPath = [DataPath, '\KEY1-rescue\R1C1'];

DataPath = 'E:\PrincetonData\20231103';
DataSubPath = [DataPath, '\woRBM\a5_bin2'];

DataPath = 'E:\PrincetonData\20231103';
DataSubPath = [DataPath, '\woRBM\a5_2_bin2'];
%% folder where data to analyze is located

DataPath = 'E:\PrincetonData\20231103';
DataSubPath = [DataPath, '\KEY1-rescue\R0C1'];

%%
ChlorPath = dir([DataSubPath,'\2D*_chlorophyll_CLAHE*']);
ChlorFile = [DataSubPath,'\',ChlorPath.name];

VenusPath = dir([DataSubPath,'\*_Venus_CLAHE*']);
VenusFile = [DataSubPath,'\',VenusPath.name];

SnapPath = dir([DataSubPath,'\*_SNAP_raw*']);
SnapFile = [DataSubPath,'\',SnapPath.name];

snap_background = 950;

%save([DataSubPath,'\PhysicalParameters.mat'],'snap_background','-append');
%% input physical parameters and other metadata
% ********Run the first time you are analyzing data************
z_size = 24;

umperpix_x = 1/6.5217;
umperpix_y = umperpix_x;
umperpix_z = 0.3;

voxel = umperpix_x*umperpix_y*umperpix_z;

%tiles = 4;

CellType = 'key1;RBCS1-Venus;KEY1-SNAP TP Constant Light';

save([DataSubPath, '\PhysicalParameters.mat']);

chlor_background = 835;
venus_background = 813;

snap_background = 950;

save([DataSubPath,'\PhysicalParameters.mat'],'venus_background','chlor_background','snap_background','-append');

%% Crop cells using FUNC_CellFinder

[CELL] = FUNC_CellFinder2D_Tiff(ChlorFile,z_size,chlor_background);

%% 
Tc = 0.03; % determined in FUNC_CellFinder. Need to find a way to store this
%pad = 50; % if imerode was large, pad>0
pad = 3;

save([DataSubPath,'\PhysicalParameters.mat'],'pad','Tc','-append');

%% Load cells 
% Do this even when returning to data
load([DataSubPath,'\PhysicalParameters.mat']);
load([DataSubPath,'\CELL.mat']);

% CELL.mat is a structure where each {} element is a different Chlamy.
% Includes .position [x, y, z], .volume [pixels], .frame [list of frames
% the cell existed in], and .bb [bounding box for the region]

% CELL.mat has been manually pruned for badly identified cells

%% Start analysis - smooth data and find threshold

%% Load data
% raw data for Venus channel to measure things about the condensates
VenusPath_raw = dir([DataSubPath,'\*_Venus_raw*']);
VenusFile_raw = [DataSubPath,'\',VenusPath_raw.name];

DATA_raw = tiffreadVolume(VenusFile_raw);

% CLAHE chlorophyll data for masking and plotting

DATA_c = tiffreadVolume(ChlorFile); % 2D

DATA_s = tiffreadVolume(SnapFile);


%% subtract off background intensity
% DATA_2_snap will be the 3D full field of view data for the SNAP channel.
% We can then crop from that data. 

DATA_c = DATA_c;

DATA_v = DATA_raw-uint16(venus_background);

DATA_2_snap = DATA_s-uint16(snap_background);

DATA_2_snap = sum(DATA_2_snap,3);

%% Import mask for condensates

% Used Threshold function in FIJI and manually threshold Venus CLAHE'd stack
% so that condensates are identified (and out-of-plane light is rejected).
% Then despeckled in FIJI. Inverted LUT and saved to
% FIJI_threshold_condensate.tif

% Load that file here

CondensatePath = dir([DataSubPath,'\*_FIJI*']);
CondensateFile = [DataSubPath,'\',CondensatePath.name];
Condensate_mask = tiffreadVolume(CondensateFile);

%% then, for each cell, c, crop DATA_2 and measure 
% (1) average SNAP intensity 
% (2) total SNAP intensity 
% (3) average intensity inside condensate 
% (4) total SNAP intensity inside condensate 

sigma = [1 1];

DATA_2 = imgaussfilt(DATA_2_snap,sigma);
f = 1;

clearvars 'PROPERTIES'

for c = 1:size(CELL,2)
    x1=ceil(CELL(c).bb(1))-pad; %BX
    y1=ceil(CELL(c).bb(2))-pad; %BY
    x2=round(x1 + CELL(c).bb(3)+1.5*pad); %BX+W
    y2=round(y1 + CELL(c).bb(4)+1.5*pad); %BY+H
    
    if x1 < 1
        x1 = 1;
    end

    if y1 <1
        y1 = 1;
    end

    if y2 > size(DATA_2,1)
        y2 = size(DATA_2,1);
    end

    if x2 > size(DATA_2,2)
        x2 = size(DATA_2,2);
    end


   % CellTemp_raw = DATA_2(y1:y2,x1:x2);
    CellTemp = Condensate_mask(y1:y2,x1:x2); % Thresholded in FIJI 
    CellTemp_c = DATA_c(y1:y2,x1:x2);

    CellTemp_snap = DATA_2_snap(y1:y2,x1:x2);

    CellTemp_Venus = DATA_v(y1:y2, x1:x2);
    
    % make mask of chlorophyll channel CellTemp_c
    CellMask = imbinarize(imgaussfilt(CellTemp_c,sigma),Tc);
    CellTemp_plot = CellTemp_snap;
     
    imcontrastscale(1) = min(CellTemp_plot,[],'all');
    imcontrastscale(2) = mean(CellTemp_plot,'all')*2;


     mip = max(CellTemp_plot,[],3); % create maximum intensity projection
     mip_c = max(CellTemp_c,[],3);
     mip_v = max(CellTemp_Venus, [], 3);

%      mis = sum(CellTemp_plot, [], 3);% create sum intensity projection
%      mis = sum(CellTemp_c, [], 3);
%      mip_v = sum(CellTemp_Venus, [], 3);


    figure1 = figure('color',[1 1 1]);
    test = imfuse(mip_c,mip_v);
    imshow(test)
    hold on



    DATA_BW_Temp = CellTemp;

    CC = bwconncomp(DATA_BW_Temp);
    S = regionprops(CC,'boundingbox','centroid','Area','EquivDiameter');
    S2 = regionprops(CC,CellTemp_snap,'PixelValues');
    S3 = regionprops(CC,CellTemp_Venus,'PixelValues');

    number = size(S,1);
    area = cat(1,S.Area);
    r = cat(1,S.EquivDiameter)./2;
    intensity_snap = cat(1,S2.PixelValues);
    intensity_Venus = cat(1,S3.PixelValues);
    bb = cat(1,S.BoundingBox);
    cent = cat(1,S.Centroid);

    index = CC.PixelIdxList; % index for pixels that are in the droplet

    for i=1:size(bb,1)
        %scatter(DATA_PLOT(i,1),DATA_PLOT(i,2),50,[0 0 1],'filled');
        rectangle('Position',[bb(i,1), bb(i,2), bb(i,3), bb(i,4)],'edgecolor',[1 1 1],'LineWidth',0.5);
    end

    title(['c=',num2str(c)])

    hold off


    % SNAP intensity inside each droplet
    MeanIntensityIn = [];
    SumIntensityIn = [];
    if number == 1
        for i = 1:size(intensity_snap,2)
            MeanIntensityIn = [MeanIntensityIn, mean(intensity_snap(:,i))];
            SumIntensityIn = [SumIntensityIn, sum(intensity_snap(:,i))];
        end
    else
        for i = 1:size(S2,1)
            MeanIntensityIn = [MeanIntensityIn, mean(S2(i).PixelValues)];
            SumIntensityIn = [SumIntensityIn, sum(S2(i).PixelValues)];
        end
    end

    % Venus intensity inside each droplet

    MeanIntensityIn_V = [];
    SumIntensityIn_V = [];
    if number == 1
        for i = 1:size(intensity_snap,2)
            MeanIntensityIn_V = [MeanIntensityIn_V, mean(intensity_snap(:,i))];
            SumIntensityIn_V = [SumIntensityIn_V, sum(intensity_snap(:,i))];
        end
    else
        for i = 1:size(S2,1)
            MeanIntensityIn_V = [MeanIntensityIn_V, mean(S3(i).PixelValues)];
            SumIntensityIn_V = [SumIntensityIn_V, sum(S3(i).PixelValues)];
        end
    end
    
    % Total SNAP intensity in the cell

    MeanIntensityCell = [];
    SumIntensityCell = [];
    for i = 1:size(intensity_snap,2)
        MeanIntensityCell = mean(CellTemp_snap,'all','omitnan');
        SumIntensityCell = sum(CellTemp_snap,'all','omitnan');
    end

    % Total Venus intensity in the cell

    MeanIntensityCell_V = [];
    SumIntensityCell_V = [];
    for i = 1:size(intensity_snap,2)
        MeanIntensityCell_V = mean(CellTemp_Venus,'all','omitnan');
        SumIntensityCell_V = sum(CellTemp_Venus,'all','omitnan');
    end

%     TooSmall = find(r<1);
%     r(TooSmall)=[];
%     number = number - length(TooSmall);
%     area(TooSmall)=[];
%     MeanIntensityIn(TooSmall) = [];
%     SumIntensityIn(TooSmall) = [];

    %if number < 20
        PROPERTIES(f).number = number;
        PROPERTIES(f).area = area.*umperpix_x*umperpix_y;
        PROPERTIES(f).radius = r.*umperpix_x;
        %PROPERTIES(c).intensity = mean(intensity{1},'all');
        PROPERTIES(f).MeanIntSNAPIn = MeanIntensityIn;
        PROPERTIES(f).MeanIntVenusIn = MeanIntensityIn_V;
        PROPERTIES(f).MeanIntSNAPCell = MeanIntensityCell;
        PROPERTIES(f).MeanIntVenusCell = MeanIntensityCell_V;
        PROPERTIES(f).SumIntSNAPIn = SumIntensityIn;
        PROPERTIES(f).SumIntSNAPCell = SumIntensityCell;
        PROPERTIES(f).SumIntVenusIn = SumIntensityIn_V;
        PROPERTIES(f).SumIntVenusCell = SumIntensityCell_V;
        f = f+1;
    %end

    %disp (f)
end

disp('paused')
pause
close all

SavePath = fileparts(ChlorFile);
save([SavePath,'\PROPERTIES_SNAP.mat'],'PROPERTIES')
% PROPERTIES will be in physical units of um or um^3
