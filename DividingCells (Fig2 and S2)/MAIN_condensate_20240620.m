%% Exponential Fit

%% folder where data to analyze is located
% ********Run the first time you are analyzing data************

%% 20240620
DataPath = 'E:\PrincetonData\20240418\WT;RBCS1-Venus';
DataSubPath = [DataPath, '\Cell_34_register'];

% DataPath = 'E:\PrincetonData\20240410\KEY1_rescue';
% DataSubPath = [DataPath, '\Cell_1_register'];

%fileseparate_one(DataSubPath);

VenusPath = [DataSubPath,'\Venus'];

%startframe = 1;
%endframe = 54;

%% 20240708 

% DataPath = 'E:\PrincetonData\fromColton\20231126';
% DataSubPath = [DataPath, '\Cell13_not-dividing'];
% 
% fileseparate_one(DataSubPath);

% DataPath = 'E:\PrincetonData\20220103';
% DataSubPath = [DataPath, '\Cell1'];

% DataPath = 'E:\PrincetonData\20220117';
% DataSubPath = [DataPath, '\k4;EPYC1_not-dividing_cell1'];

% DataPath = 'E:\PrincetonData\20220209\WT;RBCS1-Venus';
% DataSubPath = [DataPath, '\Cell7_not-dividing'];

% DataPath = 'E:\PrincetonData\20221005';
% DataSubPath = [DataPath, '\key1_no-rescue_cell_1'];

% DataPath = 'E:\PrincetonData\20240418\WT;RBCS1-Venus';
% DataSubPath = [DataPath, '\Cell_43_register'];

DataPath = 'E:\PrincetonData\20240410\KEY1_rescue';
DataSubPath = [DataPath, '\Cell_1_register'];

%fileseparate_one(DataSubPath);

VenusPath = [DataSubPath,'\Venus'];
%% input physical parameters and other metadata
% Run the first time you are analyzing data
z_size = 26;


FrameInterval = [1.5]; %in minutes

umperpix_x = 1/6.5217;
umperpix_y = umperpix_x;
umperpix_z = 0.5;

DivisionFrames = [ 1 28];

initial_cell_number = 1;

endframe = 28;

%save physical parameters in a mat structure to be loaded later

T=0.00015;

voxel = umperpix_x*umperpix_y*umperpix_z;

background = 800;
%T= 0.009; %input this based on the above tests
CellType = 'key1-1;EPYC1-Venus';

save([DataSubPath, '\PhysicalParameters.mat']);


%% import data
% START HERE IF YOU ARE COMING BACK TO DATA***************

%%
load([DataSubPath,'\PhysicalParameters.mat']);
VenusPath = [DataSubPath,'\Venus'];
[DATA] = import3D(VenusPath,z_size);

DATA = DATA(:,:,:,1:endframe);

%% Try to correct for bleaching by fitting to an exponential

% (1) Only include points that are inside the cell

% create a mask from the chlorophyll channel (exclude regions outside the
% cell)

ChlorPath = [DataSubPath,'\chlorophyll'];

[DATA_c] = import3D(ChlorPath,z_size);


DATA_c = DATA_c(:,:,:,1:endframe);
% normalize each frame by mean intensity

INTENSITY_c = [];
startframe = 1;

startframe=1;
for i = startframe:endframe
    intensity = mean(DATA_c(:,:,:,i),'all');
    max_intensity = max(DATA_c(:,:,:,i),[],'all');
    INTENSITY_c = [INTENSITY_c intensity];
end

%DATA_2=DATA;
DATA_2c = DATA_c;

for i = 1:endframe
    DATA_2c(:,:,:,i) = DATA_c(:,:,:,i)./(INTENSITY_c(i));
end

sigma = [3 3 1]; % width of gaussian
DATA_g = imgaussfilt3(DATA_c(:,:,:,10),sigma);

% If you don't want to blur 
%DATA_g = DATA;
%close all

% We want to threshold the image so that bright things are 1s and
% everything else is 0s. We start by trying to use the MATLAB function
% graythresh. 

%T=graythresh(DATA_g(:,:,:,1)); 
%Tc= 0.008;
Tc = 0.004;
DATA_BW_Temp = imbinarize(DATA_g(:,:,:,1),Tc);

% plotting to check if the threshold is good
figure0 = figure('color', [1 1 1]);
imshow(DATA_BW_Temp(:,:,10,1))

figure1=figure('color',[1 1 1]);
imcube=DATA_2c(:,:,:,1);

cmp = redpeachblue(256);
% Generate your volume object
V = volshow(DATA_g(:,:,:,1),...
    'Renderer', 'MaximumIntensityProjection',...
    'Colormap', cmp,...
    'BackgroundColor', [1, 1, 1]);

% define Data_mask and store Cell size
DATA_mask = DATA_c;
CELL_VOLUME = [];
for i = 1:size(DATA_c,4)
    DATA_BW_Temp = imbinarize(imgaussfilt3(DATA_c(:,:,:,i),sigma)); 
    CC = bwconncomp(DATA_BW_Temp);
    S = regionprops3(CC,'centroid','volume','EquivDiameter');
    %cell_volume = sum(cat(1,4/3*pi*S.EquivDiameter.^3));
    cell_volume = sum(cat(1, S.Volume));
    CELL_VOLUME = [CELL_VOLUME, cell_volume];
    DATA_mask(:,:,:,i) = DATA_BW_Temp;
end

% (2) Mask data so that any padding from registration is NaN

DATA_cell = double(DATA.*DATA_mask);
DATA_cell(DATA_cell==0) = NaN;


% (3) Find the average intensity vs time
totalframes = endframe;

INTENSITY = [];

for i = 1:totalframes
    intensity = mean(DATA_cell(:,:,:,i)-background,'all','omitnan');
    INTENSITY = [INTENSITY intensity];
end


% (4) Fit an exponential decay model to the average intensity 
% Be sure to exclude frames or values with total black due to registration

time = (1:totalframes);
time = reshape(time, [],1);
INTENSITY = reshape(INTENSITY,[],1);

[xData, yData] = prepareCurveData( time, INTENSITY);

% exponential fit with offset
ft = fittype('a*exp(-b*t) + c', 'independent', 't');
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf 0 -Inf];
%opts.StartPoint = [0.27692298496089 0.0461713906311539 0.0971317812358475];

% Fit model to data.
decayModel = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( decayModel, xData, yData );
legend( h, 'INTENSITY vs. time', 'single exponential', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'time', 'Interpreter', 'none' );
ylabel( 'INTENSITY', 'Interpreter', 'none' );
grid on

% Extract the decay parameters
a = decayModel.a;
b = decayModel.b;
c = decayModel.c;

% Calculate the correction factor for each frame
correctionFactor = a * exp(-b * (time)) + c;
correctionFactor = correctionFactor ;%/ correctionFactor(1);

% Apply the correction to each frame
DATA_2 = DATA_cell;

%comment out if no bleach correction is needed
for i = 1:totalframes
    DATA_2(:,:,:,i) = (DATA_cell(:,:,:,i)-background)./ correctionFactor(i);
end

correctedINTENSITY = [];
for i = 1:totalframes
    intensity = mean(DATA_2(:,:,:,i),'all','omitnan');
    correctedINTENSITY = [correctedINTENSITY intensity];
end

correctedINTENSITY = reshape(correctedINTENSITY,[],1);


correctedINTENSITY = INTENSITY./correctionFactor;
% Save the corrected image stack
% save('correctedImageStack.mat', 'correctedImageStack');

% Display the original and corrected intensity profiles
figure;
plot(time, INTENSITY, 'b-', 'LineWidth', 2);
hold on;
plot(time, correctedINTENSITY, 'r--', 'LineWidth', 2);



% Plot the exponential fit
fitIntensity = a * exp(-b * time) + c;
%plot(time, fitIntensity, 'k-.', 'LineWidth', 2);
plot(decayModel);
legend('Original', 'Corrected', 'Exponential Fit');
xlabel('Frame');
ylabel('Average Intensity');
title('Photobleaching Correction');
grid on;



DATA_2 = uint16(DATA_2);


%%
% Depending on image collection, it can be useful to gaussian blur the
% image (usually resonant scan no post processing)
sigma = [3 3 1]; % width of gaussian
DATA_g = imgaussfilt3(DATA_2(:,:,:,1),sigma);

% If you don't want to blur 
%DATA_g = DATA;
close all
% We want to threshold the image so that bright things are 1s and
% everything else is 0s. We start by trying to use the MATLAB function
% graythresh. 
%T=graythresh(DATA_g(:,:,:,1)); 
%  T= 0.00015; % default
% T=0.00015;
%T = 0.00003;
%T = 8.0*10^(-5);
%T=0.00004;
T= 0.00006;
DATA_BW_Temp = imbinarize(DATA_g(:,:,:,1),T);

figure0 = figure('color',[1 1 1]);
volshow(DATA_BW_Temp);

figure1=figure('color',[0 0 0]);
imcube=DATA(:,:,:,10);

cmp = redpeachblue(256);
% Generate your volume object
V = volshow(DATA_g(:,:,:,1),...
    'Renderer', 'MaximumIntensityProjection',...
    'Colormap', cmp,...
    'BackgroundColor', [1, 1, 1]);

%% When satisfied, save T to Physical Parameters file
% Only on first run through

%background = 17;

save([DataSubPath,'\PhysicalParameters.mat'], 'T','Tc','sigma','-append');

close all


%% Measure size of droplets

DATA_3 = DATA;

%DATA = DATA_2;

%VOLUME will be an nx2 array with the [frame,volume] for each found
%condensate
VOLUME = [];
RADIUS = [];
for i = 1:size(DATA_2,4)
    %DATA_BW_Temp = imbinarize(DATA(:,:,:,i),T);
     DATA_BW_Temp = imbinarize(imgaussfilt3(DATA_2(:,:,:,i),sigma),T);  % use this if taken in resonant scanning mode to blur out shot noise
   % DATA_BW_Temp = imbinarize(imgaussfilt3(DATA_2(:,:,:,i),sigma));     
    CC = bwconncomp(DATA_BW_Temp);
    S = regionprops3(CC,'centroid','volume','EquivDiameter'); 
    cent = cat(1,S.Centroid);
    vol = cat(1,S.Volume);
    r = cat(1,S.EquivDiameter./2);
    time = i*ones(size(vol));
    VOLUME = vertcat(VOLUME, [time,vol]);
    RADIUS = vertcat(RADIUS, [time,r]);
end


%get rid of very small volumes?
threshold = 0;
VOLUME_CLEAN = VOLUME;
VOLUME_CLEAN(VOLUME_CLEAN(:,2)<threshold,:)=[]; 


% plot average values vs time
MEAN_RADIUS = [];
MEAN_VOLUME = [];
TOTAL_VOLUME = [];
STD_VOLUME = [];
MAX_VOLUME = [];
NUMBER = [];
RADII = {};
startframe = 1;
%endframe = max(VOLUME_CLEAN(:,1));
endframe = max(VOLUME_CLEAN(:,1));
%endframe = 100;
for i = startframe:endframe 
    frame_volumes = VOLUME_CLEAN(VOLUME_CLEAN(:,1)==i,:);
    frame_radius = RADIUS(RADIUS(:,1)==i,:);
    mean_radius = mean(frame_radius(:,2));
    mean_volume = mean(frame_volumes(:,2));
    total_volume = sum(frame_volumes(:,2));
    std_volume = std(frame_volumes(:,2));
    max_volume = max(frame_volumes(:,2));
    if isempty(max_volume) == 1
        max_volume = NaN;
    end
    number = size(frame_volumes,1);
    MEAN_VOLUME = [MEAN_VOLUME, mean_volume];
    STD_VOLUME = [STD_VOLUME, std_volume];
    MAX_VOLUME = [MAX_VOLUME, max_volume];
    TOTAL_VOLUME = [TOTAL_VOLUME, total_volume];
    NUMBER = [NUMBER, number];
    MEAN_RADIUS = [MEAN_RADIUS, mean_radius];
    RADII{i} = frame_radius;
end

% set the x axis as time in units of minutes
if size(FrameInterval,2) == 1; 
    % For uniform time interval: 
    DivisionFramesMin = DivisionFrames.*FrameInterval;
    time = linspace(0,endframe-startframe,endframe-startframe+1).*FrameInterval;
    time = time-DivisionFramesMin(2);
else
    % For non-uniform time interval

    % set up linearly spaced time index
    time = linspace(0,endframe-startframe,endframe-startframe+1);
    
    % offset by the first division
    % time = time-DivisionFrames(2)+1; 
    
    % convert to minutes
    i=1;
    t1 = FrameIntervalPhaseLength(i);
    time(1:t1) = time(1:t1).*FrameInterval(i);
    time(t1+1:end) = time(t1).*ones(size(time(t1+1:end)))+FrameInterval(i+1).*(time(t1+1:end)-t1);

    % off set by the first division

    time = time - time(DivisionFrames(2));
    DivisionFramesMin = time(DivisionFrames);
    ZeroShiftMin = find(time==0);
end
%%

startindex = startframe;
endindex = endframe;

voxel = umperpix_x*umperpix_y*umperpix_z;

x = time(startindex:endindex);
y = MEAN_RADIUS(startindex:endindex).*voxel^(1/3);

%DivisionFramesMin = DivisionFrames.*FrameInterval;

figure1 = figure('color',[1 1 1]);
%plot(time,MEAN_RADIUS.*voxel^(1/3));
scatter(x,y,'filled');
title('<R>')
set(gca,'fontsize',14);
xlabel('Time (min)');
ylabel('Radius (\mum)')
hold on
for i =2:length(DivisionFramesMin)-1
    plot([DivisionFramesMin(i)-DivisionFramesMin(2),DivisionFramesMin(i)-DivisionFramesMin(2)],[0 6000],'color','red');
end
ylim([0 max(y)+0.1*max(y)])
box


x = time;
y = NUMBER;

figure2 = figure('color',[1 1 1]);
plot(x,y,'linewidth',1.3);
xlabel('Time (min)');
ylabel('Number')
set(gca,'fontsize',14);
hold on
for i =2:length(DivisionFrames)-1
    plot([DivisionFramesMin(i)-DivisionFramesMin(2),DivisionFramesMin(i)-DivisionFramesMin(2)],[0 6000],'color','red');
end
ylim([0 max(y)+0.1*max(y)])
box


y = MEAN_VOLUME.*voxel;

figure3 = figure('color',[1 1 1]);
scatter(x,y,'filled');
xlabel('Time (min)')
ylabel('Mean Volume (\mum^3)')
set(gca,'fontsize',14);
hold on
for i =2:length(DivisionFrames)-1
    plot([DivisionFramesMin(i)-DivisionFramesMin(2),DivisionFramesMin(i)-DivisionFramesMin(2)],[0 6000],'color','red');
end
ylim([0 max(y)+0.1*max(y)])
box

y = STD_VOLUME.*voxel;

figure4 = figure('color',[1 1 1]);
scatter(x,y,'filled');
xlabel('Time (min)')
ylabel('Standard Deviation Volume (\mum^3)')
set(gca,'fontsize',14);
hold on
for i =2:length(DivisionFrames)-1
    plot([DivisionFramesMin(i)-DivisionFramesMin(2),DivisionFramesMin(i)-DivisionFramesMin(2)],[0 6000],'color','red');
end
ylim([0 max(y)+0.1*max(y)])
box

y=MAX_VOLUME.*voxel;

figure5 = figure('color',[1 1 1]);
scatter(x,y,'filled');
xlabel('Time (min)')
ylabel('Max Volume (\mum^3)')
set(gca,'fontsize',14);
hold on
for i =2:length(DivisionFrames)-1
    plot([DivisionFramesMin(i)-DivisionFramesMin(2),DivisionFramesMin(i)-DivisionFramesMin(2)],[0 6000],'color','red');
end
ylim([0 max(y)+0.1*max(y)])
box


CELL_NUMBER = ones(size(DATA_2,4));
%DivisionFramesFrame=DivisionFrames./1.5;
DivisionFramesFrame = DivisionFrames;
for i =1:length(DivisionFrames)-1
    CELL_NUMBER(DivisionFramesFrame(i):DivisionFramesFrame(i+1)) = 2^(i-1).*CELL_NUMBER(DivisionFramesFrame(i):DivisionFramesFrame(i+1));
end

y = NUMBER./CELL_NUMBER(startframe:endframe);

figure6 = figure('color',[1 1 1]);
plot(x,y,'linewidth',1.3);
%scatter(time,NUMBER./CELL_NUMBER,'filled');
xlabel('Time (min)');
ylabel('Number/cell')
set(gca,'fontsize',14);
hold on
for i =2:length(DivisionFrames)-1
    plot([DivisionFramesMin(i)-DivisionFramesMin(2),DivisionFramesMin(i)-DivisionFramesMin(2)],[0 6000],'color','red');
end
ylim([0 max(y)+0.1*max(y)])
box

%% Save Plots

FigName = { 'MeanRadius','Number','MeanVolume','StdevVolume','MaxVolume','NumberPerCell'};

SaveDirectory = [DataSubPath, '\plots_20240708_exp2'];

mkdir(SaveDirectory);

% SaveName = [SaveDirectory,'\Radii'];
% saveas(figure0,SaveName);

i=1;
 SaveName = [SaveDirectory,'\',FigName{i}];
 saveas(figure1,SaveName);
i=2;
 SaveName = [SaveDirectory,'\',FigName{i}];
 saveas(figure2,SaveName);
 i=3;
 SaveName = [SaveDirectory,'\',FigName{i}];
 saveas(figure3,SaveName);
 i=4;
 SaveName = [SaveDirectory,'\',FigName{i}];
 saveas(figure4,SaveName);
 i=5;
 SaveName = [SaveDirectory,'\',FigName{i}];
 saveas(figure5,SaveName);
 i=6;
 SaveName = [SaveDirectory,'\',FigName{i}];
 saveas(figure6,SaveName);
 


save([SaveDirectory,'\plot_data.mat'],'MEAN_RADIUS','NUMBER','MEAN_VOLUME','STD_VOLUME','MAX_VOLUME','NUMBER','time','CELL_NUMBER');

close all



%%

% create INTENSITY:a list of the intensities of each region identified by
%regionprops3.
INTENSITY = cell(1,endframe);
RADIUS1 = cell(1,endframe);
for i = 2:endframe
    %DATA_BW_Temp = imbinarize(DATA(:,:,:,i),T);
    DATA_BW_Temp = imbinarize(imgaussfilt3(DATA_2(:,:,:,i),sigma),T);  % use this if taken in resonant scanning mode to blur out shot noise
    % DATA_BW_Temp = imbinarize(imgaussfilt3(DATA_2(:,:,:,i),sigma)); 
    CC = bwconncomp(DATA_BW_Temp);
    DATA_Temp = DATA_2(:,:,:,i); % Using the raw, un-normalized data for this
    S = regionprops3(CC,'centroid','volume','EquivDiameter'); 
    S2 = regionprops3(CC,DATA_Temp,'VoxelValues');
    intensity = cat(1,S2.VoxelValues);
    r = cat(1,S.EquivDiameter./2);
    RADIUS1{i} = r;
    INTENSITY{i} = intensity; 
end


% calculate the mean intensity of all droplets/regions in a particular frame

MEAN_INTENSITY=[];
SUM_INTENSITY = [];
INTEGRATED_INTENSITY = [];
for i = 1:length(INTENSITY)
    intensity = INTENSITY{i};
    radii = RADIUS1{i};
    condensed_volume = sum(4/3*pi*radii.^3);
    intensity2 = [];
    if iscell(intensity)==1
        for j = 1:length(intensity)
            intensity2 = vertcat(intensity2, intensity{j});
        end
        mean_intensity = mean(intensity2,'omitnan');
        sum_intensity = sum(intensity2);
        integrated_intensity = condensed_volume.*mean_intensity;
        MEAN_INTENSITY = [MEAN_INTENSITY, mean_intensity];
        SUM_INTENSITY = [SUM_INTENSITY, sum_intensity];
        if isempty(RADIUS1{i}) == 1;
            INTEGRATED_INTENSITY = [INTEGRATED_INTENSITY, 0];
        else
            INTEGRATED_INTENSITY = [INTEGRATED_INTENSITY, integrated_intensity];
        end
    else
        for j = 1:length(intensity)
            intensity2 = vertcat(intensity2, intensity(j));
        end
        mean_intensity = mean(intensity2,'omitnan');
        sum_intensity = sum(intensity2);
        MEAN_INTENSITY = [MEAN_INTENSITY, mean_intensity];
        SUM_INTENSITY = [SUM_INTENSITY, sum_intensity];
    end
end

% plot average intensity inside regions vs time

% time = linspace(0,endframe-startframe,endframe-startframe+1).*FrameInterval;
% 
% time = time-DivisionFramesMin(2);

x = time(startindex:endindex);

y = MEAN_INTENSITY(startindex:endindex);


figure7 = figure('color',[1 1 1]);
plot(x,y);
title('mean intensity inside')
xlabel('Time (min)')
ylabel('Mean Intensity Inside Droplet (AU)')
set(gca,'fontsize',14);
box
hold on
for i =2:length(DivisionFrames)-1
    plot([DivisionFramesMin(i)-DivisionFramesMin(2),DivisionFramesMin(i)-DivisionFramesMin(2)],[0 6000],'color','red');
end
ylim([0 max(y)+0.1*max(y)])
box

% calculate average intensity outside the regionprops
% should exclude values that are not at all fluorescent

MEAN_INTENSITY_OUT = [];
SUM_INTENSITY_OUT = [];
for i = 1:endframe
    %DATA_BW_Temp = imbinarize(DATA(:,:,:,i),T);
    DATA_BW_Temp = imbinarize(imgaussfilt3(DATA_2(:,:,:,i),sigma),T);  % use this if taken in resonant scanning mode to blur out shot noise
    CC = bwconncomp(DATA_BW_Temp);
    DATA_Temp = DATA_2(:,:,:,i);
    DATA_Temp = double(DATA_Temp);
    index = CC.PixelIdxList;

    for j = 1:length(index)
        idx = index{j};
        DATA_Temp(idx) = 0;
    end
    DATA_Temp = DATA_Temp.*double(DATA_mask(:,:,:,i));
    DATA_Temp(DATA_Temp<=0.5) = NaN;
    mean_intensity = mean(mean(mean(DATA_Temp,'omitnan'),'omitnan'),'omitnan');
    sum_intensity = sum(DATA_Temp ,'all','omitnan');
    MEAN_INTENSITY_OUT = [MEAN_INTENSITY_OUT, mean_intensity];
    SUM_INTENSITY_OUT = [SUM_INTENSITY_OUT, sum_intensity];
end

% plot average intensity outside regions vs time

y = MEAN_INTENSITY_OUT(startindex:endindex);


figure8 = figure('color',[1 1 1]);
plot(x, y);
title('mean intensity outside')
set(gca,'fontsize',14);
xlabel('Time (min)')
ylabel('Mean Intensity Outside Droplet (AU)')
box
hold on
for i =2:length(DivisionFrames)-1
    plot([DivisionFramesMin(i)-DivisionFramesMin(2),DivisionFramesMin(i)-DivisionFramesMin(2)],[0 6000],'color','red');
end
ylim([0 max(y)+0.1*max(y)])
box


% plot average ratio of inside to outside over time

y = (MEAN_INTENSITY(startindex:endindex))./(MEAN_INTENSITY_OUT(startindex:endindex));

figure9 = figure('color',[1 1 1]);
plot(x, y);
title('inside:outside intensity')
xlabel('Time (min)')
ylabel('Ratio of Inside to Outisde Intensity')
set(gca,'fontsize',14);
hold on
for i =2:length(DivisionFrames)-1
    plot([DivisionFramesMin(i)-DivisionFramesMin(2),DivisionFramesMin(i)-DivisionFramesMin(2)],[0 6000],'color','red');
end
ylim([min(y)-0.1*min(y) max(y)+0.1*max(y)])
box


% Plot average total intensity of the cell 

%background = 20;
%dark = 15;
MEAN_INTENSITY_TOTAL = [];
MAX_INTENSITY = [];
for i = 1:endframe
    DATA_Temp = DATA_2(:,:,:,i);
    DATA_Temp = double(DATA_Temp).*double(DATA_mask(:,:,:,i));
    DATA_Temp(DATA_Temp==0) = NaN;
    mean_intensity = mean(mean(mean(DATA_Temp,'omitnan'),'omitnan'),'omitnan');
    max_intensity = max(DATA_Temp,[],'all');
    MAX_INTENSITY = [MAX_INTENSITY, max_intensity];
    MEAN_INTENSITY_TOTAL = [MEAN_INTENSITY_TOTAL, mean_intensity];
end

y = (MEAN_INTENSITY_TOTAL(startindex:endindex));

figure10 = figure('color',[1 1 1]);
plot(x(2:end), y(2:end),'linewidth',1.3);
title('total intensity')
xlabel('Time (min)')
ylabel('Total Intensity (AU)')
set(gca,'fontsize',14);
hold on
for i =2:length(DivisionFrames)-1
    plot([DivisionFramesMin(i)-DivisionFramesMin(2),DivisionFramesMin(i)-DivisionFramesMin(2)],[0 6000],'color','red');
end
ylim([0 max(y)+0.1*max(y)])
box


y = (MAX_INTENSITY(startindex:endindex));

figure11 = figure('color',[1 1 1]);
plot(x(2:end), y(2:end),'linewidth',1.3);
title('max intensity')
xlabel('Time (min)')
ylabel('Max Intensity (AU)')
set(gca,'fontsize',14);
hold on
for i =2:length(DivisionFrames)-1
    plot([DivisionFramesMin(i)-DivisionFramesMin(2),DivisionFramesMin(i)-DivisionFramesMin(2)],[0 6000],'color','red');
end
ylim([0 max(y)+0.1*max(y)])
box


y = SUM_INTENSITY(startindex:endindex);


figure12 = figure('color',[1 1 1]);
plot(x,y);
title('sum intensity inside')
xlabel('Time (min)')
ylabel('Mean Intensity Inside Droplet (AU)')
set(gca,'fontsize',14);
box
hold on
for i =2:length(DivisionFrames)-1
    plot([DivisionFramesMin(i)-DivisionFramesMin(2),DivisionFramesMin(i)-DivisionFramesMin(2)],[0 6000],'color','red');
end
ylim([0 max(y)+0.1*max(y)])
box


y = SUM_INTENSITY_OUT(startindex:endindex);


figure13 = figure('color',[1 1 1]);
plot(x, y);
title('sum intensity outside')
set(gca,'fontsize',14);
xlabel('Time (min)')
ylabel('Mean Intensity Outside Droplet (AU)')
box
hold on
for i =2:length(DivisionFrames)-1
    plot([DivisionFramesMin(i)-DivisionFramesMin(2),DivisionFramesMin(i)-DivisionFramesMin(2)],[0 6000],'color','red');
end
ylim([0 max(y)+0.1*max(y)])
box

y = (SUM_INTENSITY(startindex:endindex))+(SUM_INTENSITY_OUT(startindex:endindex));

figure14 = figure('color',[1 1 1]);
plot(x, y);
title('sum intensity')
xlabel('Time (min)')
ylabel('Ratio of Inside to Outisde Intensity')
set(gca,'fontsize',14);
hold on
for i =2:length(DivisionFrames)-1
    plot([DivisionFramesMin(i)-DivisionFramesMin(2),DivisionFramesMin(i)-DivisionFramesMin(2)],[0 6000],'color','red');
end
ylim([min(y)-0.1*min(y) max(y)+0.1*max(y)])
box


y = SUM_INTENSITY(startindex:endindex)./(SUM_INTENSITY(startindex:endindex)+SUM_INTENSITY_OUT(startindex:endindex)).*100;
figure15 = figure('color',[1 1 1]);
plot(x(2:end),y(2:end));
set(gca,'fontsize',14);
box
hold on
for i =2:length(DivisionFrames)-1
plot([DivisionFramesMin(i)-DivisionFramesMin(2),DivisionFramesMin(i)-DivisionFramesMin(2)],[0 60000],'color','black');
end
ylim([0 max(y)+0.1*max(y)])
box
y = SUM_INTENSITY_OUT(startindex:endindex)./(SUM_INTENSITY(startindex:endindex)+SUM_INTENSITY_OUT(startindex:endindex)).*100;
hold on
plot(x(2:end), y(2:end));
title('Percent Fluorescence')
set(gca,'fontsize',14);
xlabel('Time (min)')
ylabel('Percent inside (outside) drops')
box
for i =2:length(DivisionFrames)-1
plot([DivisionFramesMin(i)-DivisionFramesMin(2),DivisionFramesMin(i)-DivisionFramesMin(2)],[0 60000],'color','black');
end
ylim([0 max(y)+0.1*max(y)])
box

% Volume fraction inside dense phase
y = TOTAL_VOLUME(startindex:endindex)./CELL_VOLUME(startindex:endindex);
figure16 = figure('color',[1 1 1]);
plot(x, y);
title('Volume Fraction')
xlabel('Time (min)')
ylabel('Volume dense phase/Cell volume')
set(gca,'fontsize',14);
hold on
for i =2:length(DivisionFrames)-1
    plot([DivisionFramesMin(i)-DivisionFramesMin(2),DivisionFramesMin(i)-DivisionFramesMin(2)],[0 6000],'color','red');
end
ylim([min(y)-0.1*min(y) max(y)+0.1*max(y)])
box

%%
FigName = { FigName{:},'InsideIntensity','OutsideIntensity','Ratio', 'MeanIntensity','MaxIntensity'};


%mkdir(SaveDirectory);


% 
 i=7;
 SaveName = [SaveDirectory,'\',FigName{i}];
 saveas(figure7,SaveName);

i=8;
 SaveName = [SaveDirectory,'\',FigName{i}];
 saveas(figure8,SaveName);
i=9;
 SaveName = [SaveDirectory,'\',FigName{i}];
 saveas(figure9,SaveName);
 i=10;
 SaveName = [SaveDirectory,'\',FigName{i}];
 saveas(figure10,SaveName);
 i=11;
 SaveName = [SaveDirectory,'\',FigName{i}];
 saveas(figure11,SaveName);

 %save variables as matlab
 
save([SaveDirectory,'\plot_data.mat'],'-append','INTENSITY','MEAN_INTENSITY_TOTAL','MAX_INTENSITY','MEAN_INTENSITY_OUT','time','TOTAL_VOLUME','CELL_VOLUME');


FigName = { FigName{:},'SumInsideIntensity','SumOutsideIntensity','SumTotalIntensity','PercentFluorescence','VolumeFraction'};

i=12;
 SaveName = [SaveDirectory,'\',FigName{i}];
 saveas(figure12,SaveName);
 i=13;
 SaveName = [SaveDirectory,'\',FigName{i}];
 saveas(figure13,SaveName);
 i=14;
 SaveName = [SaveDirectory,'\',FigName{i}];
 saveas(figure14,SaveName);
  i=15;
 SaveName = [SaveDirectory,'\',FigName{i}];
 saveas(figure15,SaveName);

  i=16;
 SaveName = [SaveDirectory,'\',FigName{i}];
 saveas(figure16,SaveName);

 save([SaveDirectory,'\plot_data.mat'],'-append','SUM_INTENSITY','SUM_INTENSITY_OUT');



%% Plot standard deviation of intensity inside the whole cell vs time


DATA_cell = DATA_mask.*DATA_2;

MEAN_INTENSITY_ALL = [];
SUM_INTENSITY_ALL = [];
STD_INTENSITY_ALL = [];
for i = 1:endframe
    %DATA_BW_Temp = imbinarize(DATA(:,:,:,i),T);
  %  DATA_BW_Temp = imbinarize(imgaussfilt3(DATA_2(:,:,:,i),sigma),T);  % use this if taken in resonant scanning mode to blur out shot noise
%     CC = bwconncomp(DATA_BW_Temp);
      DATA_Temp = DATA_2(:,:,:,i);
      DATA_Temp = double(DATA_Temp);
%     index = CC.PixelIdxList;
% 
%     for j = 1:length(index)
%         idx = index{j};
%         DATA_Temp(idx) = 0;
%     end
    DATA_Temp = DATA_Temp.*double(DATA_mask(:,:,:,i));
%     DATA_Temp(DATA_Temp<=0.5) = NaN;
    mean_intensity = mean(mean(mean(DATA_Temp,'omitnan'),'omitnan'),'omitnan');
    sum_intensity = sum(DATA_Temp ,'all','omitnan');
    std_intensity = std(std(std(DATA_Temp,'omitnan'),'omitnan'),'omitnan');
    MEAN_INTENSITY_ALL = [MEAN_INTENSITY_ALL, mean_intensity];
    SUM_INTENSITY_ALL = [SUM_INTENSITY_ALL, sum_intensity];
    STD_INTENSITY_ALL = [STD_INTENSITY_ALL, std_intensity];
end


y = STD_INTENSITY_ALL/mean(MEAN_INTENSITY_ALL);

figure17 = figure('color',[1 1 1]);
plot(time, y);
title('Standard Deviation of Intensity')
xlabel('Time (min)')
ylabel('Standard Deviation of Intensity/ Mean of Intensity')
set(gca,'fontsize',14);
hold on
for i =2:length(DivisionFrames)-1
    plot([DivisionFramesMin(i)-DivisionFramesMin(2),DivisionFramesMin(i)-DivisionFramesMin(2)],[0 6000],'color','red');
end
ylim([min(y)-0.1*min(y) max(y)+0.1*max(y)])
box

%%

FigName = { FigName{:},'StandardDeviationIntensity'};

SaveDirectory = [DataSubPath, '\plots_20240419'];

  i=17;
 SaveName = [SaveDirectory,'\',FigName{i}];
 saveas(figure17,SaveName);

 save([SaveDirectory,'\plot_data.mat'],'-append','MEAN_INTENSITY_ALL','SUM_INTENSITY_ALL','STD_INTENSITY_ALL');


%% Max Intensity Projection Analysis

mip = max(DATA_2,[],3);

DATA_BW_Temp = imbinarize(mip);


%% Plot histogram of intensity inside cell for each time point

HistIntensityPath = [SaveDirectory,'\IntensityHistograms2'];

mkdir(HistIntensityPath);

BinEdges = linspace(0,20,41);

for i=1:10:endframe
    DATA_Temp = DATA_2(:,:,:,i);
    DATA_Temp = double(DATA_Temp).*double(DATA_mask(:,:,:,i));
    DATA_Temp(DATA_Temp==0) = NaN;
    DATA_Temp = reshape(DATA_Temp,[],1);
    figure12 = figure('color', [1 1 1]);
    histogram(DATA_Temp, BinEdges,'Normalization','Probability');
%     xlim([0 20]);
%     ylim([0 0.1])
    %ylim([0 22000]);
    SaveName = [HistIntensityPath,'\hist_',sprintf('%03d',i),'.png'];
    saveas(figure12,SaveName);
    close all
end














%% bleach measured from non dividing cells
ND = dir([DataPath,'\*not-dividing*']);

A = [];
B = [];
C = [];

for i = 1:length(ND)
    load([DataPath,'\',ND(i).name,'\bleaching_fit.mat']);
    A = [A, fitresult.a];
    B = [B, fitresult.b];
    C = [C, fitresult.c];
end

mA = mean(A);
mB = mean(B);
mC = mean(C);



%% correct intensity for bleaching over time
DATA2=DATA;
for t = 1:length(DATA)
    mIntensity = max(DATA(:,:,:,t));
    DATA2(:,:,:,t) = DATA(:,:,:,t)./mIntensity.*(mA*exp(mB*t)+mC);
end
%%
%check that this worked
INTENSITY = [];
MAX_INTENSITY = [];
startframe=1;

for i = startframe:endframe
    intensity = mean(DATA(:,:,:,i),'all');
    max_intensity = max(DATA(:,:,:,i),[],'all');
    INTENSITY = [INTENSITY intensity];
    MAX_INTENSITY = [MAX_INTENSITY, max_intensity];
end

norm_mean = max(INTENSITY);

time = linspace(0,endframe-startframe,endframe-startframe+1).*FrameInterval;
figure1 = figure ('color', [1 1 1]);
scatter(time, INTENSITY/norm_mean);
title('Mean Intensity')
figure2 = figure('color',[1 1 1]);
scatter(time, MAX_INTENSITY);
title('Max Intensity')
