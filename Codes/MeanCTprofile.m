% This code plots a mean CT number profile of the cropped DICOM volume

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

%% create empty matrices to store results
CTmean =[];
SD =[];
Ulim =[];
Llim =[];
            
for i = 1:dim(3) % going through each image slice
              
img = input (:,:,i);
bi_img = imbinarize (img); % create a binary image
ind = find (bi_img == 1); % find all the 1s
                
ct_no = double (img(ind));
ct_mean = mean(ct_no); % calculate mean
ct_sd = std (ct_no); % calculate STD
                
CTmean = [CTmean; ct_mean]; % store the results
SD = [SD; ct_sd]; % store the results
                
ul = ct_mean + ct_sd; % calculate upper limit
ll = ct_mean - ct_sd; % calculate lower limit
                
Ulim = [Ulim; ul]; % store the results
Llim = [Llim; ll]; % store the results
end

%% identify the threshold for detecting changes in data
change_thresh = Threshold;
            
mean_cal = all_mean (~isnan(all_mean)); % removing Nan values
            
% identifying points of changes
pts = findchangepts (mean_cal, 'Statistic', 'rms',...
                'MinThreshold',change_thresh);
            
% convert the points to depth
pts_depth = pts*info.SliceThickness/10;
