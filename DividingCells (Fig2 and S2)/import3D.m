function [DATA] = import3D(folder,z_size);
% Basic import 3D data code from Bezia 2022-01-17
% Get Image Data, 
% need to specify z_size, the number of z-slices in the stack, beforehand
% Need to specify folder, the path holding tiff files beforehand
files = dir([folder '\*.tif']); 
for i = 1:length(files)
    z = mod(i,z_size);
    if z == 0
       z = z_size;
    end
    t = ceil(i/z_size);
    DATA(:,:,z,t) = imread([folder ,'\', files(i).name]);
end