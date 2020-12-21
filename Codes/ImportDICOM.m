
%  This code loads DICOM files into the workspace
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     This file is part of Core-CT. Copyright (C) 2020  Yu Ting Yan
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.
%      
%     Please report any bug, error or suggestion to yuting004@e.ntu.edu.sg
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%% making sure that the DICOM is loaded the right way up
% Original = flip(Original,3); % flip if required
