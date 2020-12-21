


%% Set the directory
path = uigetdir;

filefolder = fullfile (path);
cd(path)

d = dir;
names = {d.name};
names(:,1:3) = [];

info = dicominfo (fullfile (filefolder, names {5})); % get the relevant DICOM info from a random image slice
numimages = length (5:length(names));

Original = zeros (info.Rows, info.Columns, numimages, class(info.BitDepth)); % creating an empty matrix with the relevant infomation


%% rearranging the 3D matrix according to the instance number

h = waitbar(0, 'Loading CT data'); % create waitbar
order = [];

% for all images in the folder
for i = 1: length (names)
    
% read the instance number
fname = fullfile(filefolder, names{i});
info = dicominfo(fname);
instance_no = info.InstanceNumber;

% insert the image in the proper location of the matrix
Original (:,:,instance_no) = uint16 (dicomread (fname));

order (instance_no) = i; % make sure that the DICOM images load in sequence

waitbar (i/length(names),h);
end

close (h) % close waitbar

%% converting from gray values to HU
Original = int16 (Original); % change from unsigned to signed integer to accomodate negative values
Original = (Original*info.RescaleSlope) + info.RescaleIntercept; % changing to HU
