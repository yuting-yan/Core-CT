% choosing the folder to save the cropped volume
dest = uigetdir;
dest_folder = NEWname;
mkdir (dest, dest_folder);
cd([dest filesep dest_folder]);
                
% the Dicom data is not saved in order so we need to load the dicom files again according to its names so that we can export it with the corresponding metadata
                
order = order'; % the order of the name according to instance number
[~,idx] = sort(order(:,1)); % sort according to the name
                
s = waitbar (0, 'Saving'); % creating waitbar
                
%% saving cropped volume as dicom files
% sort the whole matrix using the sort indices
                    
sortedmat = sortedDICOM(:,:,idx); %sortedDICOM is the new cropped DICOM volume
                    
for i = 1:length (idx)
BaseName= NEWname;
filename =[BaseName,num2str(i)];
info2 = dicominfo (fullfile (filefolder,names {i})); % assign DICOM metadata to the correct image slice
dicomwrite (sortedmat(:,:,i),filename, info2);
                        
waitbar (i/length(names),s);
end

close (s); % close waitbar once all data are saved
msgbox ('Cropped CT images saved','success');
