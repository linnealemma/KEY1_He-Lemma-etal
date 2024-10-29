function fileseparate_one(directory) 
%This program separates different channels of fluorescence imaging into
%separate folders 
%%




    directory1=directory;
  
    %make a folder for Chlorophyll images 
    Chlordirectory=strcat(directory1,'\chlorophyll');

    %make a folder for Rubisco images
    Rubiscodirectory=strcat(directory1,'\Venus');

    mkdir(Chlordirectory);
    mkdir(Rubiscodirectory);

    %find all the images in each channel
    Drub=dir(strcat(directory1,'\*_c*1*.tif'));
    Dchlor=dir(strcat(directory1,'\*_c*2*.tif'));



    for i=1:size(Drub,1);
        file=strcat(directory1,'\',Drub(i).name);
        newfile=strcat(Rubiscodirectory);

        movefile(file,newfile);
    end

    clearvars file newfile

        %KEEP BLUE: if you don't want to delete blue pictures
    
    for j=1:size(Dchlor,1);
        fileb=strcat(directory1,'\',Dchlor(j).name);
        newfileb=Chlordirectory;

        movefile(fileb,newfileb);
    end
    














