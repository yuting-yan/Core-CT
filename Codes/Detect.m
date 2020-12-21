            
%% identify the threshold for detecting changes in data
change_thresh = Threshold;
            
mean_cal = all_mean (~isnan(all_mean)); % removing Nan values
            
% identifying points of changes
pts = findchangepts (mean_cal, 'Statistic', 'rms',...
                'MinThreshold',change_thresh);
            
% convert the points to depth
pts_depth = pts*info.SliceThickness/10;
