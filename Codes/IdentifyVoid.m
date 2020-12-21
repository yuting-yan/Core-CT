% convert the cropped image to binary
level = graythresh(img_cropped);
BW_cropped = imbinarize (img_cropped, level);
            
% calculating the area of the cropped image
totalarea = sum(sum(BW_cropped)); % Total number of white pixels
area_mm = totalarea * (info.PixelSpacing(1) *info.PixelSpacing(2)); % convert to mm^2
area_cm = area_mm/100; % convert to cm^2
            
img_filled = imfill(BW_cropped,'holes'); % filling in the hole of the image (will show a white circle)
            
% calculating the area if the void is filled up
totalarea_filled = sum(sum(img_filled)); % Total number of white pixels
area_filled_mm = totalarea_filled * (info.PixelSpacing (1)*info.PixelSpacing (2)); % convert to mm^2
area_filled_cm = area_filled_mm/100; % convert to cm^2

% Showing the void area in false color
imshowpair (BW_cropped,img_filled, 'falsecolor',...
                'Parent',app.area_plot);
            title (app.area_plot, 'Identified void');
            
%% Calculating volume of void
% cropping the volume (only use the image slices that has the void)
new_vol = subvolume (Original, [nan nan nan nan Top_slice Bottom_slice]);

vol = [];
filled_vol = [];
            
% calculate the area of void per image slice
for i = 1:((Bottom_slice - Top_slice)+1)

image = new_vol (:,:,i);
                
BW_cropped = imbinarize (image); % convert the cropped image to binary
totalarea = sum(sum (BW_cropped)); % Total number of white pixels
                
img_filled = imfill(BW_cropped,'holes'); % filling in the hole of the image (will show a white circle)
               
% imshowpair (BW_cropped, img_filled) % comparing the images
% calculating the area of the cropped image
                
totalarea_filled = sum(sum(img_filled)); % Total number of white pixels
                
area_diff = totalarea_filled - totalarea;
                
filled_vol (:,:,i) = img_filled;
vol(i) = area_diff;
end
            
%%
vol_sum = sum (vol);
            
% computing the volume
pixelvol = info.PixelSpacing (1) * info.PixelSpacing (2) *info.SliceThickness; % actual volume of each pixel
vol = pixelvol*vol_sum;
vol_cm = vol/1000;
