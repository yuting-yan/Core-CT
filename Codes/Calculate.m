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
