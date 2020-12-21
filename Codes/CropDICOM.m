%% For a circle mask
% size of the circle that I want
Rthres = 170; 

% make an x-vector the size of one dimension of your dataset
x=1:size(dicom_test,1); 
x=x-mean(x); % subtract the mean to ensure x=0 is in the middle
y=1:size(dicom_test,2); % make an y-vector the size of the other dimension of your dataset
y=y-mean(y); % subtract the mean to ensure y=0 is in the middle

% mask is 1 if the radius is smaller than the threshold, 0 otherwise
mask = uint16(sqrt(x'.^2+y.^2)<=Rthres);

%% For a rectangular mask
x1 = r_x; % r_x and r_y are bottom left coordinate of rectangle
x2 = r_x + r_width; % r_width is width of rectangle
y1 = r_y;
y2 = r_y + r_height; % r_height is height of rectangle
                
dim = size (Original); % size of DICOM volume
rect_mask = zeros (round (dim(1)),round(dim(2))); % generate grid of ones
rect_mask (y1:y2,x1:x2) = 1; % rect Mask(Y values, X values)
rect_mask = logical(rect_mask);
rect_mask = double (rect_mask); % convert from logical to numerical
rect_mask (rect_mask == 0 ) = nan; % assign all 0 to nan (all values out of ROI = nan)

%% Cropping the DICOM volume using mask created
cropDICOM = bsxfun(@times, Original, cast(mask, class(Original)));
