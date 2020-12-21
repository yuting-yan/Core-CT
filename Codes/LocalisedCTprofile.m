% This code plots a CT number profile according to a user-defined line on a CT image slice

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

L = drawline (Image_plot, 'InteractionsAllowed','none'); % draw line on selected image
hold (Image_plot, 'on');
            
% get the positions of the line
% L.Position: [x1 y1; x2 y2];
x_line = [L.Position(1,1), L.Position(2,1)];
y_line = [L.Position(1,2), L.Position(2,2)];
            
% finding the distance of the line in mm
x_dist = (x_line(2) - x_line(1))*info.PixelSpacing(1)/10;
y_dist = (y_line(2) - y_line(1))*info.SliceThickness/10;
dist = sqrt((x_dist)^2 + (y_dist)^2);
            
CT_plot = subplot(1,2,2,'Parent',app.Figure3Panel);
set (CT_plot, 'position', [0.6    0.06    0.35    0.90]); % set position of plot
disableDefaultInteractivity(CT_plot);
            
% get the profile of the CT number from the line drawn
[cx, cy, c] = improfile(disp_image,x_line, y_line);
app.scale = dist/length(c);
app.c_depth = (1:length(c))*scale;
 
% plot the localised CT number profile
plot (CT_plot,c, c_depth);
set (CT_plot, 'YDir', 'reverse')
grid (CT_plot, 'minor');

%% Identifying the peaks
% threshold value is the minimum peak distance in pixels
% converting the threshold in cm to pixel

p_thresh = MinDistance/scale; % MinDistance is the average distance between each peak
            
[pks, locs] = findpeaks (c, 'MinPeakDistance', p_thresh);
c_loc = locs*scale;
            
% number of peaks
peaks_no = length(pks);
            
%% finding distance between each peak

            
peaks_no = length (app.pks);
for i = 1:(peaks_no-1)
                
% finding distance between each peak in cm
dist_diff = (locs(i+1) - locs (i))*scale;
             
end
