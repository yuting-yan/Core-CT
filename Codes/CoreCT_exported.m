classdef CoreCT_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        ImportyourdataPanel             matlab.ui.container.Panel
        ImportButton                    matlab.ui.control.Button
        ResetButton                     matlab.ui.control.Button
        FilenameEditFieldLabel          matlab.ui.control.Label
        FilenameEditField               matlab.ui.control.EditField
        FlipuprightCheckBox             matlab.ui.control.CheckBox
        SelecttypeofmaskPanel           matlab.ui.container.Panel
        ButtonGroup                     matlab.ui.container.ButtonGroup
        CircleButton                    matlab.ui.control.RadioButton
        RectangleButton                 matlab.ui.control.RadioButton
        DrawmaskButton                  matlab.ui.control.Button
        Panel                           matlab.ui.container.Panel
        RadiuscmEditFieldLabel          matlab.ui.control.Label
        RadiuscmEditField               matlab.ui.control.NumericEditField
        Panel_2                         matlab.ui.container.Panel
        HeightcmEditFieldLabel          matlab.ui.control.Label
        HeightcmEditField               matlab.ui.control.NumericEditField
        WidthcmEditFieldLabel           matlab.ui.control.Label
        WidthcmEditField                matlab.ui.control.NumericEditField
        CropthedataPanel                matlab.ui.container.Panel
        CropButton                      matlab.ui.control.Button
        ViewallslicesButton             matlab.ui.control.Button
        VolumeViewerButton              matlab.ui.control.Button
        AdjustifneededPanel             matlab.ui.container.Panel
        XangleSliderLabel               matlab.ui.control.Label
        XangleSlider                    matlab.ui.control.Slider
        YangleSliderLabel               matlab.ui.control.Label
        YangleSlider                    matlab.ui.control.Slider
        PlotthedataPanel                matlab.ui.container.Panel
        PlotButton                      matlab.ui.control.Button
        TopmmEditFieldLabel             matlab.ui.control.Label
        TopmmEditField                  matlab.ui.control.NumericEditField
        EndmmEditFieldLabel             matlab.ui.control.Label
        EndmmEditField                  matlab.ui.control.NumericEditField
        ExportPanel                     matlab.ui.container.Panel
        SaveDICOMfileasEditFieldLabel   matlab.ui.control.Label
        SaveDICOMfileasEditField        matlab.ui.control.EditField
        ExportDICOMButton               matlab.ui.control.Button
        SavedataButton                  matlab.ui.control.Button
        SaveFigure1Button               matlab.ui.control.Button
        SaveFigure3Button               matlab.ui.control.Button
        GetpositionsButton              matlab.ui.control.Button
        Figure1Panel                    matlab.ui.container.Panel
        Figure3Panel                    matlab.ui.container.Panel
        AnalysisPanel                   matlab.ui.container.Panel
        DetectButton                    matlab.ui.control.Button
        NumberofcontactsEditFieldLabel  matlab.ui.control.Label
        NumberofcontactsEditField       matlab.ui.control.NumericEditField
        ThresholdEditFieldLabel         matlab.ui.control.Label
        ThresholdEditField              matlab.ui.control.NumericEditField
        IdentifyingpeaksPanel           matlab.ui.container.Panel
        NoofPeaksEditFieldLabel         matlab.ui.control.Label
        NoofPeaksEditField              matlab.ui.control.NumericEditField
        IdentifypeaksButton             matlab.ui.control.Button
        MinDistancecmEditFieldLabel     matlab.ui.control.Label
        MinDistancecmEditField          matlab.ui.control.NumericEditField
        CalculatedistanceButton         matlab.ui.control.Button
        ShowpeaksButton                 matlab.ui.control.Button
        SelectplanePanel                matlab.ui.container.Panel
        ButtonGroup_2                   matlab.ui.container.ButtonGroup
        CoronalButton                   matlab.ui.control.RadioButton
        SagitalButton                   matlab.ui.control.RadioButton
        DrawlineButton                  matlab.ui.control.Button
        SlicenumberSpinnerLabel         matlab.ui.control.Label
        SlicenumberSpinner              matlab.ui.control.Spinner
        ColorCheckBox                   matlab.ui.control.CheckBox
        Figure2Panel_2                  matlab.ui.container.Panel
        VoidareaEditFieldLabel          matlab.ui.control.Label
        VoidareaEditField               matlab.ui.control.NumericEditField
        AxialsliceSpinnerLabel          matlab.ui.control.Label
        AxialsliceSpinner               matlab.ui.control.Spinner
        IdentifyareaButton              matlab.ui.control.Button
        CalculateVolumeButton           matlab.ui.control.Button
        VoidvolumeEditFieldLabel        matlab.ui.control.Label
        VoidvolumeEditField             matlab.ui.control.NumericEditField
        TopsliceEditFieldLabel          matlab.ui.control.Label
        TopsliceEditField               matlab.ui.control.NumericEditField
        BottomsliceEditFieldLabel       matlab.ui.control.Label
        BottomsliceEditField            matlab.ui.control.NumericEditField
        ViewVolumeButton                matlab.ui.control.Button
    end

%     
%     MATLAB application for the quantitative analysis of CT scan images for geosciences
%     Copyright (C) 2020  Yu Ting Yan
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
    properties (Access = private)
        Property % Description
        
        info = []; % metadata of the dicom images
        names = []; % the names of all the dicom images
        filefolder =[]; % file of the dicom folder
        imported = []; % 3D matrix arranged by instance number
        sorted = []; % correct order of 3D matrix
        order = []; % order of the instance number
        
        cbar = []; % colorbar for the Axial image
        pages = []; % number of image slices
        depth =[]; % for plotting the graphs
        dim = []; % size of the matrix
        
        % different views of the core
        Axial = []; % cross sectional image
        Axial_plot =[]; % plot of the axial image
        coronal = []; % Coronal image
        Cor_plot =[]; % plot of the coronal image
        sagital = []; % Sagittal image
        Sag_plot = []; % plot of the sagital image
        
        circle_mask = []; % mask of the circle created (2D)
        Rthres = []; % threshold of radius of the circle mask
        c_center =[]; % center of circle
        
        rect_mask =[]; % mask of the rectangle created (2D)
        r_x = []; % x-coordinate of the left bottom point of rectangle
        r_y = []; % y-coordinate of the left bottom point of rectangle
        r_width = []; % width of rectangle
        r_height  = []; % height of rectangle
        r_center = []; % center of rectangle
        
        circle_cor = []; % coronal view of the circle mask
        circle_sag = []; % sagital view of the circle mask
        rect_cor = []; % coronal view of the rectangle mask
        rect_sag =[]; % sagital view of the rectangle mask
        
        cor_pos = []; % position of first and last slice of cropped volume in coronal
        sag_pos = []; % position of first and last slice of cropped volume in sagittal
        
        cylinder = []; % 3D circle mask
        cylinder_sorted =[]; % 3D matrix of the cropped cylinder
        cube = []; % 3D rectangle mask
        cube_sorted =[]; % 3D matrix of the cropped cuboid
        
        all_mean =[]; % mean of the CT numbers
        all_sd =[]; % standard deviation of the CT numbers
        ulim =[]; % upper limit (ulim) = mean + std
        llim =[]; % lower limit (llim) = mean - std
        hist_data = []; % data to plot the distribution plot
        
        Image_plot = []; 
        CT_plot =[]; % plot of mean CT number profile
        Dis_plot = []; % cropped coronal/sagital image
        disp_axial = []; % cropped axial image
        
        area_plot = [];
        void_vol = [];
        
        pts_depth =[]; % depth of the change position
        
        disp_image = []; % display image of cropped DICOM
        c = []; % CT profile of line
        cx =[]; % spatial coordinates of the CT number profile
        cy =[];
        c_depth = []; % depth of CT profile
        c_smooth = []; % smoothed CT profile
        scale = []; % scale to convert image pixel to mm
        pks = []; % peaks identified
        locs = []; % position of peaks
        
    end
    
    methods (Access = public)
        
        function info = get_info(app)
            path = uigetdir; % choosing the dicom folder
            app.filefolder = fullfile(path);
            cd (path);
            
            d = dir;
            app.names = {d.name};
            app.names(:,1:5) = []; % removing the first three rows which are not essential
            
            info = dicominfo (fullfile (app.filefolder, app.names {6}));
        end
        
        function new_mask = rotate_mask (app,input_mask, rot_x, rot_y)
            
            
            emp_mask = uint16(ones(size(app.sorted)));
            
            % making the 3d mask
            mask3d = bsxfun(@times, emp_mask, cast(input_mask, class(emp_mask)));
            
            %% rotating around X-axis
            % creating transformation matrix
            x_t = [cos(rot_x)  0      -sin(rot_x)   0
                0             1              0     0
                sin(rot_x)    0       cos(rot_x)   0
                0             0              0     1];
            
            tform_x = affine3d(x_t);
            X_rotated = imwarp(mask3d,tform_x);
            
            X_permute = permute (X_rotated, [2 1 3]);
            
            %% Rotating around Y-axis
            % creating transformation matrix
            y_t = [cos(rot_y)  0      -sin(rot_y)   0
                0             1              0     0
                sin(rot_y)    0       cos(rot_y)   0
                0             0              0     1];
            
            tform_y = affine3d (y_t);
            
            
            %%
            Straight = imwarp (X_permute, tform_y);
            
            %% trimming the adjusted mask to the original size
            mask_straight = permute (Straight, [2 1 3]);
            
            new_size = size (mask_straight);
            
            x_center = (new_size (1)/2);
            y_center = (new_size(2)/2);
            z_center = (new_size(3)/2);
            
            newx1 = x_center - (app.dim(1)/2);
            newx2 = x_center + (app.dim(1)/2);
            newy1 = y_center - (app.dim(2)/2);
            newy2 = y_center + (app.dim(2)/2);
            newz1 = z_center - (app.dim(3)/2);
            newz2 = z_center + (app.dim(3)/2);
            
            new_mask = mask_straight(newx1+1:newx2, newy1+1:newy2, newz1+1:newz2);
        end
        
        function [CTmean, SD, Ulim, Llim] = Calculate(app,input)
            
            CTmean =[];
            SD =[];
            Ulim =[];
            Llim =[];
            
            for i = 1:app.pages
%                 ct_no = reshape (input(:,:,i),1,[]); % rearranging all elements
%                 ct_no = double (ct_no);
%                 ct_no (ct_no == 0) = nan; % removing the 0 values
                
                img = input (:,:,i);
                bi_img = imbinarize (img); % create a binary image
                ind = find (bi_img == 1); % find all the 1s
                
                ct_no = double (img(ind));
              
                ct_mean = mean(ct_no);
                ct_sd = std (ct_no);
                
                CTmean = [CTmean; ct_mean];
                SD = [SD; ct_sd];
                
                ul = ct_mean + ct_sd;
                ll = ct_mean - ct_sd;
                
                Ulim = [Ulim; ul];
                Llim = [Llim; ll];
            end
        end
        
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: ImportButton
        function ImportButtonPushed(app, event)
          tic
            try
                % loading the metadata of the dicom file
                app.info = app.get_info();
                
                % printing the name of the dicom file chosen
                [~,name,~]=fileparts(pwd);
                app.FilenameEditField.Value = name;
                
                % creating a matrix of zero to store all the dicom files
                app.imported =  zeros (app.info.Rows, app.info.Columns, ...
                    length (app.names), class(app.info.BitDepth));
                
                % creating a waitbar
                h = waitbar (0, 'Loading CT data');
                
                for i = 1: length (app.names)
                    
                    % read the instance number
                    fname = fullfile(app.filefolder, app.names{i});
                    app.info = dicominfo(fname);
                    instance_no = app.info.InstanceNumber;
                    
                    % load the dicom files into the matrix created according to
                    % the instance number so that the CT images are in order
                    app.imported (:,:,instance_no) = int16 (dicomread (fname));
                    
                    % arrangement of the name in order of instance number
                    app.order (instance_no) = i;
                    waitbar (i/length(app.names),h);
                    
                end
                
                % close the waitbar once the data loads finish
                close (h);
                app.imported = int16 (app.imported);
                %% converting from gray values to HU
                % change from unsigned to signed integer to accomodate negative values
                app.imported = (app.imported * app.info.RescaleSlope) + app.info.RescaleIntercept; % changing to HU
                
                  value = app.FlipuprightCheckBox.Value;
            
            switch (value)
                case 0
                 app.sorted = app.imported;
                    
                case 1
                 app.sorted = flip(app.imported,3);
            end
            
                % Getting the cross-section image from the middle of the core
                app.dim = size (app.sorted);
                half_image = round (length(app.names)/2); % finding the middle slice
                app.Axial = app.imported (:,:, half_image);
                
                % plotting the image slice
                % permute the matrix to get cross-sectional view down core
                % coronal image
                side = permute (app.sorted, [3 2 1]);
                middle_cor = round (size (side,2)/2);
                app.coronal = side (:,:, middle_cor);
                
                % sagital image
                side2 = permute (app.sorted, [3 1 2]);
                middle_sag = round (size (side2, 2)/2);
                app.sagital = side2 (:,:,middle_sag);
                
                % creating figures for the images and plot them
                app.Figure1Panel.AutoResizeChildren = 'off';
                app.Axial_plot = subplot(3,2,[1,2],'Parent',app.Figure1Panel);
                app.Cor_plot = subplot(3,2,[3,5],'Parent',app.Figure1Panel);
                app.Sag_plot = subplot (3,2,[4,6], 'Parent', app.Figure1Panel);
                
                % plotting the axial image
                imshow (app.Axial, 'Parent',app.Axial_plot, ...
                    'Display',[], 'YData',[0,app.dim(1)],'XData',[0,(app.dim(2))]);
                disableDefaultInteractivity(app.Axial_plot);
                Ori_size = get(app.Axial_plot, 'Position');
                app.cbar = colorbar (app.Axial_plot,'Location', 'eastoutside'); % show the colorbar
                app.cbar.Label.String = 'Grayscale value'; % label the colorbar
                % pos = get(c,'Position');
                % c.Label.Position = [pos(1) pos(2)+1.8]; % changing the colorbar label position
                set (app.Axial_plot, 'Position', Ori_size)
                
                ylabel (app.Axial_plot, 'y');
                xlabel (app.Axial_plot, 'x');
                title (app.Axial_plot, 'Axial View');
                
                % plotting the coronal image
                imagesc(app.coronal, 'Parent', app.Cor_plot, ...
                    'YData',[0 length(app.sorted)]);
                ylim (app.Cor_plot, [0, length(app.sorted)]);
                colormap (app.Cor_plot, gray);
                disableDefaultInteractivity(app.Cor_plot);
                set (app.Cor_plot, 'xtick',[], 'xticklabel',[],...
                    'ytick',[], 'yticklabel',[]);
                ylabel (app.Cor_plot, 'z');
                xlabel (app.Cor_plot, 'x');
                title (app.Cor_plot, 'Coronal View');
                
                % plotting the sagittal image
                imagesc(app.sagital, 'Parent', app.Sag_plot, ...
                    'YData',[0 length(app.sorted)]);
                ylim (app.Sag_plot, [0, length(app.sorted)]);
                colormap (app.Sag_plot, gray);
                disableDefaultInteractivity(app.Sag_plot);
                set (app.Sag_plot, 'xtick', [], 'xticklabel',[],...
                    'ytick',[], 'yticklabel',[]);
                ylabel (app.Sag_plot, 'z');
                xlabel (app.Sag_plot, 'y');
                title (app.Sag_plot, 'Sagital View');
                
            catch
                % load a error message if there are no dicom files found
                errordlg('No DICOM files found','File Error');
            end
            
            % fixing the positions of the axes
            set (app.Axial_plot, 'position', [0.1300    0.7093    0.5497    0.257]);
            set (app.Cor_plot, 'position',[0.1300    0.0300    0.35   0.6]);
            set (app.Sag_plot, 'position', [0.5703    0.0300    0.35    0.6]);
            
            %% printing the top and end values on the editfield based on the DICOM images loaded
            app.TopmmEditField.Value = 0;
            app.EndmmEditField.Value = app.dim(3) * app.info.SliceThickness; % values are in mm
            app.SlicenumberSpinner.Value = middle_cor;
            
            % allow user to draw mask only after CT images are imported
            app.DrawmaskButton.Enable = 'on';
            app.CircleButton.Enable = 'on';
            app.RectangleButton.Enable = 'on';
            app.RadiuscmEditField.Enable = 'on';
            app.RadiuscmEditFieldLabel.Enable = 'on';
            app.FlipuprightCheckBox.Enable = 'on';
      
            toc
        end

        % Value changed function: FlipuprightCheckBox
        function FlipuprightCheckBoxValueChanged(app, event)
            tic
            value = app.FlipuprightCheckBox.Value;
            
            switch (value)
                case 0
                 app.sorted = app.imported;
                    
                case 1
                 app.sorted = flip(app.imported,3);
             end
            
            % Getting the cross-section image from the middle of the core
                app.dim = size (app.sorted);
                half_image = round (length(app.names)/2); % finding the middle slice
                app.Axial = app.sorted (:,:, half_image);
                
                % plotting the image slice
                % permute the matrix to get cross-sectional view down core
                % coronal image
                side = permute (app.sorted, [3 2 1]);
                middle_cor = round (size (side,2)/2);
                app.coronal = side (:,:, middle_cor);
                
                % sagital image
                side2 = permute (app.sorted, [3 1 2]);
                middle_sag = round (size (side2, 2)/2);
                app.sagital = side2 (:,:,middle_sag);
                
                % creating figures for the images and plot them
                app.Figure1Panel.AutoResizeChildren = 'off';
                app.Axial_plot = subplot(3,2,[1,2],'Parent',app.Figure1Panel);
                app.Cor_plot = subplot(3,2,[3,5],'Parent',app.Figure1Panel);
                app.Sag_plot = subplot (3,2,[4,6], 'Parent', app.Figure1Panel);
                
                % plotting the axial image
                imshow (app.Axial, 'Parent',app.Axial_plot, ...
                    'Display',[], 'YData',[0,app.dim(1)],'XData',[0,(app.dim(2))]);
                disableDefaultInteractivity(app.Axial_plot);
                Ori_size = get(app.Axial_plot, 'Position');
                app.cbar = colorbar (app.Axial_plot,'Location', 'eastoutside'); % show the colorbar
                app.cbar.Label.String = 'Grayscale value'; % label the colorbar
                % pos = get(c,'Position');
                % c.Label.Position = [pos(1) pos(2)+1.8]; % changing the colorbar label position
                set (app.Axial_plot, 'Position', Ori_size)
                
                ylabel (app.Axial_plot, 'y');
                xlabel (app.Axial_plot, 'x');
                title (app.Axial_plot, 'Axial View');
                
                % plotting the coronal image
                imagesc(app.coronal, 'Parent', app.Cor_plot, ...
                    'YData',[0 length(app.sorted)]);
                ylim (app.Cor_plot, [0, length(app.sorted)]);
                colormap (app.Cor_plot, gray);
                disableDefaultInteractivity(app.Cor_plot);
                set (app.Cor_plot, 'xtick',[], 'xticklabel',[],...
                    'ytick',[], 'yticklabel',[]);
                ylabel (app.Cor_plot, 'z');
                xlabel (app.Cor_plot, 'x');
                title (app.Cor_plot, 'Coronal View');
                
                % plotting the sagittal image
                imagesc(app.sagital, 'Parent', app.Sag_plot, ...
                    'YData',[0 length(app.sorted)]);
                ylim (app.Sag_plot, [0, length(app.sorted)]);
                colormap (app.Sag_plot, gray);
                disableDefaultInteractivity(app.Sag_plot);
                set (app.Sag_plot, 'xtick', [], 'xticklabel',[],...
                    'ytick',[], 'yticklabel',[]);
                ylabel (app.Sag_plot, 'z');
                xlabel (app.Sag_plot, 'y');
                title (app.Sag_plot, 'Sagital View');
                
                  % fixing the positions of the axes
            set (app.Axial_plot, 'position', [0.1300    0.7093    0.5497    0.257]);
            set (app.Cor_plot, 'position',[0.1300    0.0300    0.35   0.6]);
            set (app.Sag_plot, 'position', [0.5703    0.0300    0.35    0.6]);
           
            toc
        end

        % Selection changed function: ButtonGroup
        function ButtonGroupSelectionChanged(app, event)
           
            if (app.CircleButton.Value)
                app.Panel.Visible = 'on';
                app.Panel_2.Visible = 'off';
            end
            
            if (app.RectangleButton.Value)
                app.Panel.Visible = 'off';
                app.Panel_2.Visible = 'on';
            end
            
        end

        % Button pushed function: DrawmaskButton
        function DrawmaskButtonPushed(app, event)
            tic
            pixel_size = app.info.PixelSpacing (1);
            
            if (app.CircleButton.Value)
                
                imshow (app.Axial, 'Parent',app.Axial_plot, ...
                    'Display',[],'YData',[0,app.dim(1)],'XData',[0,(app.dim(2))]);
                app.cbar = colorbar (app.Axial_plot,'Location', 'eastoutside');
                app.cbar.Label.String = 'Grayscale value';
                ylabel (app.Axial_plot, 'y');
                xlabel (app.Axial_plot, 'x');
                title (app.Axial_plot, 'Axial View');
                app.c = drawcircle (app.Axial_plot, 'InteractionsAllowed','none');
                
                app.Rthres = app.c.Radius;
                app.c_center = app.c.Center;
                
                % printing the radius value on the edit field
                app.RadiuscmEditField.Value = (app.Rthres*pixel_size)/10;
                app.RadiuscmEditField.Enable = 'on';
                app.RadiuscmEditField.Editable = 'off';
                app.RadiuscmEditFieldLabel.Enable = 'on';
                
                % make an x-vector the size of one dimension of your dataset
                x = 1:size(app.sorted,1);
                x = x-app.c_center (2); % subtract the mean to ensure x=0 is in the middle
                
                % make an y-vector the size of the other dimension of your dataset
                y = 1:size(app.sorted,2);
                y = y-app.c_center (1); % subtract the mean to ensure y=0 is in the middle
                
                % mask is 1 if the radius is smaller than the threshold, 0 otherwise
                app.circle_mask = uint16(sqrt(x'.^2+y.^2)<= app.Rthres);
                app.circle_mask = double (app.circle_mask); % convert from logical to numerical
                app.circle_mask (app.circle_mask == 0 ) = nan; % assign all 0 to nan (all values out of ROI = nan) 
                
                cx1 = app.c_center(1) - app.Rthres;
                cx2 = app.c_center(1) + app.Rthres;
                
                sx1 = app.c_center(2) - app.Rthres;
                sx2 = app.c_center (2) + app.Rthres;
                
                cy1 = 0;
                cy2 = length (app.sorted);
                
                % plotting the cylinder for the coronal view
                imagesc(app.coronal, 'Parent', app.Cor_plot, ...
                    'YData',[0 length(app.sorted)]);
                ylim (app.Cor_plot, [0, length(app.sorted)]);
                colormap (app.Cor_plot, gray);
                disableDefaultInteractivity(app.Cor_plot);
                set (app.Cor_plot, 'xtick',[], 'xticklabel',[],...
                    'ytick',[], 'yticklabel',[]);
                ylabel (app.Cor_plot, 'z');
                xlabel (app.Cor_plot, 'x');
                title (app.Cor_plot, 'Coronal View');
                
                hold(app.Cor_plot, 'on');
                app.circle_cor = plot (app.Cor_plot,...
                    [cx1 cx1 cx2 cx2 cx1], [cy1 cy2 cy2 cy1 cy1],...
                    'LineWidth',3, 'Color', [0 0.4470 0.7410]);
                
                % plotting the cylinder for the sagital view
                imagesc(app.sagital, 'Parent', app.Sag_plot, ...
                    'YData',[0 length(app.sorted)]);
                ylim (app.Sag_plot, [0, length(app.sorted)]);
                colormap (app.Sag_plot, gray);
                disableDefaultInteractivity(app.Sag_plot);
                set (app.Sag_plot, 'xtick', [], 'xticklabel',[],...
                    'ytick',[], 'yticklabel',[]);
                ylabel (app.Sag_plot, 'z');
                xlabel (app.Sag_plot, 'y');
                title (app.Sag_plot, 'Sagital View');
                
                hold (app.Sag_plot, 'on');
                app.circle_sag = plot(app.Sag_plot, [sx1 sx1 sx2 sx2 sx1],...
                    [cy1 cy2 cy2 cy1 cy1],'LineWidth',3,...
                    'Color', [0 0.4470 0.7410]);
            end
            
            
            if (app.RectangleButton.Value)
                
                imshow (app.Axial, 'Parent',app.Axial_plot, ...
                    'Display',[],'YData',[0,app.dim(1)],'XData',[0,(app.dim(2))]);
                app.cbar = colorbar (app.Axial_plot,'Location', 'eastoutside');
                app.cbar.Label.String = 'Grayscale value';
                ylabel (app.Axial_plot, 'y');
                xlabel (app.Axial_plot, 'x');
                title (app.Axial_plot, 'Axial View');
                r = drawrectangle (app.Axial_plot, 'InteractionsAllowed','none');
                
                app.WidthcmEditField.Value = (r.Position (3)*pixel_size)/10;
                app.HeightcmEditField.Value = (r.Position (4)*pixel_size)/10;
                
                % position of 4 points of the cuboid
                app.r_x = r.Position (1);
                app.r_y = r.Position (2);
                app.r_width = r.Position (3);
                app.r_height = r.Position (4);
                
                app.r_center = [(r.Position(1)+(r.Position(3)/2)), ...
                    (r.Position(2)+(r.Position (4)/2))];
                
                x1 = app.r_x;
                x2 = app.r_x + app.r_width;
                y1 = app.r_y;
                y2 = app.r_y + app.r_height;
                
                % create a rectangular mask
                app.rect_mask = zeros (round (app.dim(1)),round(app.dim(2))); % generate grid of ones
                app.rect_mask (y1:y2,x1:x2) = 1; % rect Mask(Y values, X values)
                app.rect_mask = logical(app.rect_mask);
                app.rect_mask = double (app.rect_mask); % convert from logical to numerical
                app.rect_mask (app.rect_mask == 0 ) = nan; % assign all 0 to nan (all values out of ROI = nan)
                
                % plotting the mask downcore on the sagital and coronal
                
                ry1 = 0;
                ry2 = length (app.sorted);
                
                % plotting the cylinder for the coronal view
                imagesc(app.coronal, 'Parent', app.Cor_plot, ...
                    'YData',[0 length(app.sorted)]);
                ylim (app.Cor_plot, [0, length(app.sorted)]);
                colormap (app.Cor_plot, gray);
                disableDefaultInteractivity(app.Cor_plot);
                set (app.Cor_plot, 'xtick',[], 'xticklabel',[],...
                    'ytick',[], 'yticklabel',[]);
                ylabel (app.Cor_plot, 'z');
                xlabel (app.Cor_plot, 'x');
                title (app.Cor_plot, 'Coronal View');
                
                hold(app.Cor_plot, 'on');
                app.rect_cor = plot(app.Cor_plot,...
                    [x1 x1 x2 x2 x1], [ry1 ry2 ry2 ry1 ry1],...
                    'Linewidth', 3, 'Color', [0 0.4470 0.7410]) ;
                plot (app.Cor_plot, app.rect_cor);
                
                % plotting the cylinder for the sagital view
                imagesc(app.sagital, 'Parent', app.Sag_plot, ...
                    'YData',[0 length(app.sorted)]);
                ylim (app.Sag_plot, [0, length(app.sorted)]);
                colormap (app.Sag_plot, gray);
                disableDefaultInteractivity(app.Sag_plot);
                set (app.Sag_plot, 'xtick', [], 'xticklabel',[],...
                    'ytick',[], 'yticklabel',[]);
                ylabel (app.Sag_plot, 'z');
                xlabel (app.Sag_plot, 'y');
                title (app.Sag_plot, 'Sagital View');
                
                hold (app.Sag_plot, 'on');
                app.rect_sag = plot(app.Sag_plot, ...
                    [y1 y1 y2 y2],[ry1 ry2 ry2 ry1], ...
                    'LineWidth',3, 'Color', [0 0.4470 0.7410]);
                plot (app.Sag_plot, app.rect_sag);
                
            end
            
            app.CropButton.Enable = 'on';
            app.XangleSlider.Enable = 'on';
            app.XangleSliderLabel.Enable = 'on';
            app.YangleSlider.Enable = 'on';
            app.YangleSliderLabel.Enable = 'on';
            
            toc
        end

        % Value changing function: XangleSlider
        function XangleSliderValueChanging(app, event)
          tic
            x_ang = event.Value;
            
            if (app.CircleButton.Value)
                
                cx1 = app.c_center(1) - app.Rthres;
                cx2 = app.c_center(1) + app.Rthres;
                
                cy1 = 0;
                cy2 = length (app.sorted);
                
                % plotting the cylinder for the coronal view
                imagesc(app.coronal, 'Parent', app.Cor_plot, ...
                    'YData',[0 length(app.sorted)]);
                ylim (app.Cor_plot, [0, length(app.sorted)]);
                colormap (app.Cor_plot, gray);
                disableDefaultInteractivity(app.Cor_plot);
                set (app.Cor_plot, 'xtick',[], 'xticklabel',[],...
                    'ytick',[], 'yticklabel',[]);
                ylabel (app.Cor_plot, 'z');
                xlabel (app.Cor_plot, 'x');
                title (app.Cor_plot, 'Coronal View');
                
                hold(app.Cor_plot, 'on');
                app.circle_cor = plot (app.Cor_plot,...
                    [cx1 cx1 cx2 cx2 cx1], [cy1 cy2 cy2 cy1 cy1],...
                    'LineWidth',3, 'Color', [0 0.4470 0.7410]);
                
                rotate(app.circle_cor,[0 0 1],x_ang);
            end
            
            if (app.RectangleButton.Value)
                x1 = app.r_x;
                x2 = app.r_x + app.r_width;
                y1 = app.r_y;
                y2 = app.r_y + app.r_height;
                
                % create a rectangular mask
                app.rect_mask = zeros (round (app.dim(1)),round(app.dim(2))); % generate grid of ones
                app.rect_mask (y1:y2,x1:x2) = 1; % rect Mask(Y values, X values)
                app.rect_mask = logical(app.rect_mask);
                
                % plotting the mask downcore on the sagital and coronal
                ry1 = 0;
                ry2 = length (app.sorted);
                
                % plotting the cylinder for the coronal view
                imagesc(app.coronal, 'Parent', app.Cor_plot, ...
                    'YData',[0 length(app.sorted)]);
                ylim (app.Cor_plot, [0, length(app.sorted)]);
                colormap (app.Cor_plot, gray);
                disableDefaultInteractivity(app.Cor_plot);
                set (app.Cor_plot, 'xtick',[], 'xticklabel',[],...
                    'ytick',[], 'yticklabel',[]);
                ylabel (app.Cor_plot, 'z');
                xlabel (app.Cor_plot, 'x');
                title (app.Cor_plot, 'Coronal View');
                
                hold(app.Cor_plot, 'on');
                app.rect_cor = plot(app.Cor_plot,...
                    [x1 x1 x2 x2 x1], [ry1 ry2 ry2 ry1 ry1],...
                    'Linewidth', 3, 'Color', [0 0.4470 0.7410]) ;
                plot (app.Cor_plot, app.rect_cor);
                
                rotate(app.rect_cor,[0 0 1],x_ang);
            end
            toc
        end

        % Value changing function: YangleSlider
        function YangleSliderValueChanging(app, event)
            y_ang = event.Value;
            
            if (app.CircleButton.Value)
                
                sx1 = app.c_center(2) - app.Rthres;
                sx2 = app.c_center (2) + app.Rthres;
                
                cy1 = 0;
                cy2 = length (app.sorted);
                
                % plotting the cylinder for the sagital view
                imagesc(app.sagital, 'Parent', app.Sag_plot, ...
                    'YData',[0 length(app.sorted)]);
                ylim (app.Sag_plot, [0, length(app.sorted)]);
                colormap (app.Sag_plot, gray);
                disableDefaultInteractivity(app.Sag_plot);
                set (app.Sag_plot, 'xtick', [], 'xticklabel',[],...
                    'ytick',[], 'yticklabel',[]);
                ylabel (app.Sag_plot, 'z');
                xlabel (app.Sag_plot, 'y');
                title (app.Sag_plot, 'Sagital View');
                
                hold (app.Sag_plot, 'on');
                app.circle_sag = plot(app.Sag_plot, [sx1 sx1 sx2 sx2 sx1],...
                    [cy1 cy2 cy2 cy1 cy1],'LineWidth',3,...
                    'Color', [0 0.4470 0.7410]);
                
                rotate(app.circle_sag,[0 0 1],y_ang);
            end
            
            if (app.RectangleButton.Value)
                x1 = app.r_x;
                x2 = app.r_x + app.r_width;
                y1 = app.r_y;
                y2 = app.r_y + app.r_height;
                
                % create a rectangular mask
                app.rect_mask = zeros (round (app.dim(1)),round(app.dim(2))); % generate grid of ones
                app.rect_mask (y1:y2,x1:x2) = 1; % rect Mask(Y values, X values)
                app.rect_mask = logical(app.rect_mask);
                
                % plotting the mask downcore on the sagital and coronal
                ry1 = 0;
                ry2 = length (app.sorted);
                
                % plotting the cylinder for the sagital view
                imagesc(app.sagital, 'Parent', app.Sag_plot, ...
                    'YData',[0 length(app.sorted)]);
                ylim (app.Sag_plot, [0, length(app.sorted)]);
                colormap (app.Sag_plot, gray);
                disableDefaultInteractivity(app.Sag_plot);
                set (app.Sag_plot, 'xtick', [], 'xticklabel',[],...
                    'ytick',[], 'yticklabel',[]);
                ylabel (app.Sag_plot, 'z');
                xlabel (app.Sag_plot, 'y');
                title (app.Sag_plot, 'Sagital View');
                
                hold (app.Sag_plot, 'on');
                app.rect_sag = plot(app.Sag_plot, ...
                    [y1 y1 y2 y2],[ry1 ry2 ry2 ry1], ...
                    'LineWidth',3, 'Color', [0 0.4470 0.7410]);
                plot (app.Sag_plot, app.rect_sag);
                
                
                rotate(app.rect_sag,[0 0 1],y_ang);
            end
        end

        % Button pushed function: CropButton
        function CropButtonPushed(app, event)
            tic
            app.pages = app.dim(3);
            app.depth = (1:app.pages);
            pixel_size = app.info.SliceThickness;
            app.AxialsliceSpinner.Limits = [0 inf];
            
            % getting the depth information in cm
            app.depth = app.depth * (pixel_size/10);
            
            % getting the angles (in radian) from the slider
            x_ang = -deg2rad(app.XangleSlider.Value);
            y_ang = -deg2rad(app.YangleSlider.Value);
            
            f = waitbar (0, 'Loading the data'); % create waitbar
            
            % if a cylinder mask is chosen
            if (app.CircleButton.Value)
                new_mask = rotate_mask (app,app.circle_mask, x_ang, y_ang);
                
                waitbar (0.5, f,'Processing');
                
                %% apply the mask to the original volume
                app.cylinder_sorted = bsxfun(@times, app.sorted, ...
                    cast(new_mask, class(app.sorted)));
                
                app.cylinder_sorted (app.cylinder_sorted ==0) = -1000;
                
                waitbar(0.9 ,f,'Finishing');
                
                % plotting the image slice
                % permute the matrix to get cross-sectional view down core
                side = permute (app.cylinder_sorted, [3 2 1]);
                
                % finding the middle slice
                mid_slice = round (app.c_center (2));
                
                % cropping to get just the core image
                image_slice = side (:,:, mid_slice);
                I = mat2gray (image_slice);
                
                
                % creating the figures and plots
                app.Figure3Panel.AutoResizeChildren = 'off';
                app.Image_plot = subplot(1,3,1,'Parent',app.Figure3Panel);
                set (app.Image_plot, 'position',[0.1  0.06  0.35  0.90]);
                
                imagesc (I, 'Parent', app.Image_plot, ...
                    'YData',[0 (max(app.pages))]);
                
                disableDefaultInteractivity(app.Image_plot);
                
                % creating the Axial image
                app.Figure2Panel_2.AutoResizeChildren = 'off';
                app.Dis_plot = subplot(1,2,1,'Parent',app.Figure2Panel_2);
                set (app.Dis_plot, 'position',[0.0  0.12  0.45  0.75]);
                
                app.disp_axial = app.cylinder_sorted (:,:,round(app.dim(3)/2));
                imshow (app.disp_axial, 'Parent',app.Dis_plot, ...
                    'Display',[], 'YData',[0,app.dim(1)],'XData',[0,(app.dim(2))]);
                title (app.Dis_plot, 'Cropped Axial image');
                
                
                % limit values of the spinner based on cropped volume
                app.cor_pos = [round(app.c_center(2) - app.Rthres) round(app.c_center(2) + app.Rthres)];
                app.sag_pos = [round(app.c_center(1) - app.Rthres) round(app.c_center(1) + app.Rthres)];
                
                
                app.SlicenumberSpinner.Enable = 'on';
                app.SlicenumberSpinnerLabel.Enable = 'on';
                app.SlicenumberSpinner.Limits = [app.cor_pos(1) app.cor_pos(2)];
                app.SlicenumberSpinner.Value = round (mean (app.cor_pos));
                
                app.AxialsliceSpinner.Enable = 'on';
                app.AxialsliceSpinnerLabel.Enable = 'on';
                app.AxialsliceSpinner.Value = round (app.dim(3)/2);
                app.AxialsliceSpinner.Limits = [1 app.dim(3)];
                
            end
            
            
            if (app.RectangleButton.Value)
                
                new_mask = rotate_mask (app,app.rect_mask, x_ang, y_ang);
                
                waitbar (0.5, f,'Processing');
                
                % cropping the matrix using the rectangle mask
                app.cube_sorted = bsxfun(@times, app.sorted, ...
                    cast(new_mask, class(app.sorted)));
                
                app.cube_sorted (app.cube_sorted ==0) = -1000;
                
                waitbar(0.9 ,f,'Finishing');
                
                % plotting the image slice
                side = permute (app.cube_sorted, [3 2 1]);
                
                mid_slice = round (app.r_center (2));
                
                % plotting the image side view
                image_slice = side (:,:, mid_slice);
                
                app.Figure3Panel.AutoResizeChildren = 'off';
                app.Image_plot = subplot(1,2,1,'Parent',app.Figure3Panel);
                set (app.Image_plot, 'position',[0.1  0.06  0.35  0.90]);
                
                imagesc (image_slice, 'Parent', app.Image_plot, ...
                    'YData',[0 (max(app.pages))]);
                
                 % creating the Axial image
                app.Figure2Panel_2.AutoResizeChildren = 'off';
                app.Dis_plot = subplot(1,2,1,'Parent',app.Figure2Panel_2);
                set (app.Dis_plot, 'position',[0.0  0.12  0.45  0.75]);
                
                app.disp_axial = app.cube_sorted (:,:,round(app.dim(3)/2));
                imshow (app.disp_axial, 'Parent',app.Dis_plot, ...
                    'Display',[], 'YData',[0,app.dim(1)],'XData',[0,(app.dim(2))]);
                title (app.Dis_plot, 'Cropped Axial image');
                
                disableDefaultInteractivity(app.Dis_plot);
                disableDefaultInteractivity(app.Image_plot);
                
                % limit values of the spinner based on cropped volume
                app.cor_pos = [round(app.r_center(2) - app.r_width) round(app.r_center(2) + app.r_width)];
                app.sag_pos = [round(app.r_center(1) - app.r_height) round(app.r_center(1) + app.r_height)];
                
                app.SlicenumberSpinner.Enable = 'on';
                app.SlicenumberSpinnerLabel.Enable = 'on';
                app.SlicenumberSpinner.Limits = [app.cor_pos(1) app.cor_pos(2)];
                app.SlicenumberSpinner.Value = round(mean (app.cor_pos));
                
                app.AxialsliceSpinner.Enable = 'on';
                app.AxialsliceSpinnerLabel.Enable = 'on';
                app.AxialsliceSpinner.Value = round (app.dim(3)/2);
                app.AxialsliceSpinner.Limits = [1 app.dim(3)];
                
            end
            
            close(f); % close the waitbar
            
            % formatting the axes for the cropped image
            ylim (app.Image_plot, [0, app.pages]);
            colormap (app.Image_plot, gray);
            title (app.Image_plot, 'Cropped CT image');
            set (app.Image_plot, 'xtick', [], 'xticklabel',[],...
                'ytick',[], 'yticklabel',[]);
            
            % enable the following buttons
            app.ViewallslicesButton.Enable = 'on';
            app.VolumeViewerButton.Enable = 'on';
            app.ExportDICOMButton.Enable = 'on';
            app.SaveDICOMfileasEditField.Enable = 'on';
            app.SaveDICOMfileasEditFieldLabel.Enable = 'on';
            app.PlotButton.Enable = 'on';
            app.TopmmEditField.Enable = 'on';
            app.TopmmEditFieldLabel.Enable = 'on';
            app.EndmmEditField.Enable = 'on';
            app.EndmmEditFieldLabel.Enable = 'on';
            app.ColorCheckBox.Enable = 'on';
            app.CoronalButton.Enable = 'on';
            app.SagitalButton.Enable = 'on';
            app.DrawlineButton.Enable = 'on';
            app.SaveFigure1Button.Enable = 'on';
            app.IdentifyareaButton.Enable = 'on';
            
            
            disableDefaultInteractivity(app.Dis_plot);
            disableDefaultInteractivity(app.Image_plot);
            toc
        end

        % Value changed function: AxialsliceSpinner
        function AxialsliceSpinnerValueChanged(app, event)
            Axial_slice = app.AxialsliceSpinner.Value;
            
                if (app.CircleButton.Value)
                    
                    app.disp_axial = app.cylinder_sorted (:,:, Axial_slice);
                end
                
                if (app.RectangleButton.Value)
                    
                    app.disp_axial = app.cube_sorted(:,:, Axial_slice);
                end
            
                app.Figure2Panel_2.AutoResizeChildren = 'off';
                app.Dis_plot = subplot(1,2,1,'Parent',app.Figure2Panel_2);
                set (app.Dis_plot, 'position',[0.01  0.12  0.45  0.75]);
                
                imshow (app.disp_axial, 'Parent',app.Dis_plot, ...
                    'Display',[], 'YData',[0,app.dim(1)],'XData',[0,(app.dim(2))]);
                
                disableDefaultInteractivity(app.Dis_plot);
        end

        % Button pushed function: IdentifyareaButton
        function IdentifyareaButtonPushed(app, event)
          tic
            if (app.CircleButton.Value)
                img_cropped = app.cylinder_sorted (:,:,app.AxialsliceSpinner.Value);
            end
            
            if (app.RectangleButton.Value)
                img_cropped = app.cube_sorted (:,:,app.AxialsliceSpinner.Value);
            end
            
            % convert the cropped image to binary
            level = graythresh(img_cropped);
            BW_cropped = imbinarize (img_cropped, level);
            
            % calculating the area of the cropped image
            totalarea = sum(sum(BW_cropped)); % Total number of white pixels
            area_mm = totalarea * ...
                (app.info.PixelSpacing(1) * app.info.PixelSpacing(2)); % convert to mm^2
            area_cm = area_mm/100; % convert to cm^2
            
            img_filled = imfill(BW_cropped,'holes'); % filling in the hole of the image (will show a white circle)
            
            % calculating the area if the void is filled up
            totalarea_filled = sum(sum(img_filled)); % Total number of white pixels
            area_filled_mm = totalarea_filled * ...
                (app.info.PixelSpacing (1)*app.info.PixelSpacing (2)); % convert to mm^2
            area_filled_cm = area_filled_mm/100; % convert to cm^2
            
            app.Figure2Panel_2.AutoResizeChildren = 'off';
            app.area_plot = subplot(1,2,2,'Parent',app.Figure2Panel_2);
            set (app.area_plot, 'position',[0.35 0.12  0.45  0.75]);
            
            
            % Showing the void area in false color
            imshowpair (BW_cropped,img_filled, 'falsecolor',...
                'Parent',app.area_plot);
            title (app.area_plot, 'Identified void');
            
             % printing the values of the calculated void area on interface
            app.VoidareaEditField.Enable = 'on';
            app.VoidareaEditFieldLabel.Enable = 'on';
            app.VoidareaEditField.Editable ='off';
            app.VoidareaEditField.Value = area_filled_cm - area_cm;
            
            disableDefaultInteractivity(app.area_plot);
            
            app.CalculateVolumeButton.Enable = 'on';
            app.TopsliceEditField.Enable = 'on';
            app.TopsliceEditFieldLabel.Enable = 'on';
            app.BottomsliceEditField.Enable = 'on';
            app.BottomsliceEditFieldLabel.Enable = 'on';
            app.VoidvolumeEditField.Enable = 'on';
            app.VoidvolumeEditFieldLabel.Enable = 'on';
            
            toc
        end

        % Button pushed function: CalculateVolumeButton
        function CalculateVolumeButtonPushed(app, event)
          tic
            Top_slice = app.TopsliceEditField.Value;
            Bottom_slice = app.BottomsliceEditField.Value;
            
                if (app.CircleButton.Value)
                 Original = app.cylinder_sorted;
                end
                
                if (app.RectangleButton.Value)
                Original = app.cube_sorted;
                end
           
                try 
            % cropping the volume
            new_vol = subvolume (Original, [nan nan nan nan Top_slice Bottom_slice]);
            
            vol = [];
            filled_vol = [];
            
            % calculate the area of void per image slice
            for i = 1:((Bottom_slice - Top_slice)+1)
                image = new_vol (:,:,i);
                
                % convert the cropped image to binary
                BW_cropped = imbinarize (image);
                
                totalarea = sum(sum (BW_cropped)); % Total number of white pixels
                
                img_filled = imfill(BW_cropped,'holes'); % filling in the hole of the image (will show a white circle)
               %  imshowpair (BW_cropped, img_filled) % comparing the images
                % calculating the area of the cropped image
                
                totalarea_filled = sum(sum(img_filled)); % Total number of white pixels
                
                area_diff = totalarea_filled - totalarea;
                
                filled_vol (:,:,i) = img_filled;
                vol(i) = area_diff;
            end
            
            %%
            vol_sum = sum (vol);
            
            % computing the volume
            pixelvol = app.info.PixelSpacing (1) * app.info.PixelSpacing (2) *app.info.SliceThickness;
            vol = pixelvol*vol_sum;
            vol_cm = vol/1000;
            app.VoidvolumeEditField.Value = vol_cm;
            
            new = imbinarize (new_vol);
            app.void_vol = filled_vol - new;
            
                
            app.ViewVolumeButton.Enable = 'on';
            
                catch 
                    % load a error message if there are no dicom files found
                errordlg('Please specify image slices','File Error');
                end
                toc
        end

        % Button pushed function: ViewVolumeButton
        function ViewVolumeButtonPushed(app, event)
     tic  
            figure % creating a figure
            V = smooth3 (app.void_vol,'box',5); % smoothing the volume
            isosurface (V);
            set(gca,'BoxStyle','full','Box','on')
            
            grid minor
            
            xt = get(gca, 'XTick'); % 'XTick' Values
            set(gca, 'XTick', xt, 'XTickLabel', xt*app.info.PixelSpacing (1))
            yt = get(gca, 'YTick'); % 'XTick' Values
            set(gca, 'YTick', yt, 'YTickLabel', yt*app.info.PixelSpacing (2))
            zt = get(gca, 'ZTick'); % 'XTick' Values
            set(gca, 'ZTick', zt, 'ZTickLabel', zt*app.info.SliceThickness)

            xlabel ('mm');
            ylabel ('mm');
            zlabel ('mm');
            
            toc
        end

        % Button pushed function: ViewallslicesButton
        function ViewallslicesButtonPushed(app, event)
            
            if (app.CircleButton.Value)
                
                % this loads volume Viewer to visualise the cropped volume
                sliceViewer (app.cylinder_sorted);
                
            end
            
            if (app.RectangleButton.Value)
                
                % this loads volume Viewer to visualise the cropped volume
                sliceViewer (app.cube_sorted);
                
            end
            
        end

        % Button pushed function: VolumeViewerButton
        function VolumeViewerButtonPushed(app, event)
            
            if (app.CircleButton.Value)
                
                % this loads volume Viewer to visualise the cropped volume
                volshow (app.cylinder_sorted);
                
            end
            
            if (app.RectangleButton.Value)
                
                % this loads volume Viewer to visualise the cropped volume
                volshow(app.cube_sorted);
                
            end
        end

        % Button pushed function: PlotButton
        function PlotButtonPushed(app, event)
            tic
            app.dim = size (app.sorted);
            app.pages = app.dim(3); % total number of image slice in original dicom
            app.depth = (1:app.pages);
            pixel_size = app.info.SliceThickness;
            
            % getting the depth information in cm
            app.depth = (app.depth)' * (pixel_size/10);
            
            f = waitbar (0.1, 'Loading the data');
            
            new = [];
            
            if (app.CircleButton.Value)
                
                [CTmean,SD, Ulim, Llim] = ...
                    Calculate(app,app.cylinder_sorted);
                
                waitbar (0.3,  f,'Calculating');
                
                app.all_mean = CTmean;
                app.all_sd = SD;
                app.ulim = int16(Ulim);
                app.llim = int16(Llim);
                
                waitbar (0.7, f, 'Processing');
            end
            
            if (app.RectangleButton.Value)
                [CTmean, SD, Ulim, Llim] = ...
                    Calculate(app,app.cube_sorted);
                
                waitbar (0.3, f, 'Calculating');
                
                app.all_mean = CTmean;
                app.all_sd = SD;
                app.ulim = (Ulim);
                app.llim = (Llim);
                
                waitbar (0.7, f, 'Processing');
            end
            
            % plotting the mean and standard deviation
            
            app.CT_plot = subplot(1,2,2,'Parent',app.Figure3Panel);
            set (app.CT_plot, 'position', [0.6    0.06    0.35    0.90]);
            
            
            % subset the data according to what was selected by the user
            top_slice = round (app.TopmmEditField.Value/ app.info.SliceThickness);
            bottom_slice = round (app.EndmmEditField.Value/app.info.SliceThickness);
            
            app.depth ([1:top_slice, bottom_slice:length(app.depth)])= [];
            app.all_mean ([1:top_slice, bottom_slice:length(app.all_mean)])= [];
            app.all_sd ([1:top_slice, bottom_slice:length(app.all_sd)])= [];
            app.ulim ([1:top_slice, bottom_slice:length(app.ulim)])= [];
            app.llim ([1:top_slice, bottom_slice:length(app.llim)])= [];
            
            % plotting the data
            plot (app.CT_plot, app.all_mean, app.depth,...
                'Color',[0, 0.4470, 0.7410], 'LineWidth', 2);
            hold (app.CT_plot, 'on');
            patch (app.CT_plot,[app.llim; app.ulim(end:-1:1); app.llim(1)],...
                [app.depth; app.depth(end:-1:1); app.depth(1)],[0, 0.4470, 0.7410],...
                'EdgeColor','none', 'FaceAlpha', 0.2);
            hold (app.CT_plot, 'on');
            plot (app.CT_plot, app.ulim, app.depth,...
                'Color',[0.3010, 0.7450, 0.9330], 'LineWidth',0.5);
            hold (app.CT_plot, 'on');
            plot (app.CT_plot, app.llim, app.depth,...
                'Color',[0.3010, 0.7450, 0.9330], 'LineWidth', 0.5);
            title (app.CT_plot, 'CT numbers')
            ylim (app.CT_plot, [0 max(app.pages*pixel_size/10)]);
            ylabel (app.CT_plot, 'Depth (cm)');
            xlabel (app.CT_plot, 'CT numbers');
            set (app.CT_plot, 'YDir', 'reverse');
            grid (app.CT_plot, 'minor');
            hold (app.CT_plot, 'off');
            
            close (f); % close the waitbar
            
            % User can only save figure 2 after all the plots are done
            app.SaveFigure3Button.Enable = 'on';
            app.SavedataButton.Enable = 'on';
            
            app.ThresholdEditField.Enable = 'on';
            app.ThresholdEditFieldLabel.Enable = 'on';
            app.DetectButton.Enable = 'on';
            app.NumberofcontactsEditField.Enable = 'on';
            app.NumberofcontactsEditFieldLabel.Enable = 'on';
            
            disableDefaultInteractivity(app.CT_plot);
            
            toc
        end

        % Button pushed function: DetectButton
        function DetectButtonPushed(app, event)
            tic
            
            % identify the threshold for detecting changes in data
            change_thresh = app.ThresholdEditField.Value;
            
            % removing Nan values
            mean_cal = app.all_mean (~isnan(app.all_mean));
            
            % identifying points of changes
            pts = findchangepts (mean_cal, 'Statistic', 'rms',...
                'MinThreshold',change_thresh);
            
            % convert the points to depth
            app.pts_depth = pts*app.info.SliceThickness/10;
            
            % plot the CT mean
            delete (app.CT_plot)
            
            app.CT_plot = subplot(1,2,2,'Parent',app.Figure3Panel);
            set (app.CT_plot, 'position', [0.6    0.06    0.35    0.90]);
            
            plot (app.CT_plot, app.all_mean, app.depth,...
                'Color',[0, 0.4470, 0.7410], 'LineWidth', 2);
            hold (app.CT_plot, 'on');
            patch (app.CT_plot,[app.llim; app.ulim(end:-1:1); app.llim(1)],...
                [app.depth; app.depth(end:-1:1); app.depth(1)],[0, 0.4470, 0.7410],...
                'EdgeColor','none', 'FaceAlpha', 0.2);
            hold (app.CT_plot, 'on');
            plot (app.CT_plot, app.ulim, app.depth,...
                'Color',[0.3010, 0.7450, 0.9330], 'LineWidth',0.5);
            hold (app.CT_plot, 'on');
            plot (app.CT_plot, app.llim, app.depth,...
                'Color',[0.3010, 0.7450, 0.9330], 'LineWidth', 0.5);
            title (app.CT_plot, 'CT numbers')
            ylim (app.CT_plot, [0 max(app.pages*app.info.SliceThickness/10)]);
            ylabel (app.CT_plot, 'Depth (cm)');
            xlabel (app.CT_plot, 'CT numbers');
            set (app.CT_plot, 'YDir', 'reverse');
            grid (app.CT_plot, 'minor');
            hold (app.CT_plot, 'off');
            
            % plotting the lines on the CT plot
            for i = 1:length(app.pts_depth)
                
                yline(app.CT_plot, app.pts_depth (i),...
                    'Color',[0.8500, 0.3250, 0.0980]);
                hold (app.CT_plot, 'on');
                
            end
            
            % display the possible number of contacts on the app
            app.NumberofcontactsEditField.Value = size(pts,1);
            
            app.GetpositionsButton.Enable = 'on';
            
            toc
        end

        % Value changed function: ColorCheckBox
        function ColorCheckBoxValueChanged(app, event)
            value = app.ColorCheckBox.Value;
            
            % change the colormap of the image plot using the checkbox
            switch (value)
                case 0
                    colormap (app.Image_plot, 'gray');
                    
                case 1
                    colormap (app.Image_plot, 'jet');
            end
            
        end

        % Selection changed function: ButtonGroup_2
        function ButtonGroup_2SelectionChanged(app, event)
            
            if (app.CoronalButton.Value)
                
                 % changing spinner limits according to cropped volume
                app.SlicenumberSpinner.Limits = [app.cor_pos(1) app.cor_pos(2)];
                app.SlicenumberSpinner.Value = round (mean (app.cor_pos));
                
                if (app.CircleButton.Value)
                    
                    side = permute (app.cylinder_sorted, [3 2 1]);
                    app.disp_image = side (:,:, app.SlicenumberSpinner.Value);
                    
                end
                
                if (app.RectangleButton.Value)
                    
                    side = permute (app.cube_sorted, [3 2 1]);
                    app.disp_image = side (:,:, app.SlicenumberSpinner.Value);
                end
                
            end
            
            if (app.SagitalButton.Value)
                
                % changing spinner limits according to cropped volume
                app.SlicenumberSpinner.Limits = [app.sag_pos(1) app.sag_pos(2)];
                app.SlicenumberSpinner.Value = round (mean (app.sag_pos));
                
                if (app.CircleButton.Value)
                    
                    side = permute (app.cylinder_sorted, [3 1 2]);
                    app.disp_image = side (:,:, app.SlicenumberSpinner.Value);
                end
                if (app.RectangleButton.Value)
                    
                    side = permute (app.cube_sorted, [3 1 2]);
                    app.disp_image = side (:,:, app.SlicenumberSpinner.Value);
                end
                
            end
            
            app.Figure3Panel.AutoResizeChildren = 'off';
            app.Image_plot = subplot(1,2,1,'Parent',app.Figure3Panel);
            set (app.Image_plot, 'position', [0.1    0.06    0.35    0.90]);
            imagesc(app.disp_image, 'Parent', app.Image_plot);
            set (app.Image_plot, 'xtick', [], 'xticklabel',[],...
                'ytick',[], 'yticklabel',[]);
            colormap (app.Image_plot, 'gray');
            
            if (app.ColorCheckBox.Value)
                colormap (app.Image_plot, 'jet');
            end
            
            app.SlicenumberSpinner.Enable = 'on';
        end

        % Value changed function: SlicenumberSpinner
        function SlicenumberSpinnerValueChanged(app, event)
            Slice_number= app.SlicenumberSpinner.Value;
            
            if (app.CoronalButton.Value)
                if (app.CircleButton.Value)
                    
                    side = permute (app.cylinder_sorted, [3 2 1]);
                    app.disp_image = side (:,:, Slice_number);
                end
                
                if (app.RectangleButton.Value)
                    
                    side = permute (app.cube_sorted, [3 2 1]);
                    app.disp_image = side (:,:, Slice_number);
                end
                
            end
            
            if (app.SagitalButton.Value)
                if (app.CircleButton.Value)
                    
                    side = permute (app.cylinder_sorted, [3 1 2]);
                    app.disp_image = side (:,:, Slice_number);
                end
                if (app.RectangleButton.Value)
                    
                    side = permute (app.cube_sorted, [3 1 2]);
                    app.disp_image = side (:,:, Slice_number);
                end
            end
            
            app.Figure3Panel.AutoResizeChildren = 'off';
            app.Image_plot = subplot(1,2,1,'Parent',app.Figure3Panel);
            set (app.Image_plot, 'position', [0.1    0.06    0.35    0.90]);
            imagesc(app.disp_image, 'Parent', app.Image_plot);
            set (app.Image_plot, 'xtick', [], 'xticklabel',[],...
                'ytick',[], 'yticklabel',[]);
            
            % set default color of the plot to be grayscale
            colormap (app.Image_plot, 'gray');
            
            % only change to jet colormap when checked
            if (app.ColorCheckBox.Value)
                colormap (app.Image_plot, 'jet');
            end
        end

        % Button pushed function: DrawlineButton
        function DrawlineButtonPushed(app, event)
            tic
            Slice_number = app.SlicenumberSpinner.Value;
            
            if (app.CoronalButton.Value)
                if (app.CircleButton.Value)
                    
                    side = permute (app.cylinder_sorted, [3 2 1]);
                    app.disp_image = side (:,:, Slice_number);
                    
                end
                
                if (app.RectangleButton.Value)
                    
                    side = permute (app.cube_sorted, [3 2 1]);
                    app.disp_image = side (:,:, Slice_number);
                end
                
            end
            
            if (app.SagitalButton.Value)
                if (app.CircleButton.Value)
                    
                    side = permute (app.cylinder_sorted, [3 1 2]);
                    app.disp_image = side (:,:, Slice_number);
                end
                if (app.RectangleButton.Value)
                    
                    side = permute (app.cube_sorted, [3 1 2]);
                    app.disp_image = side (:,:, Slice_number);
                end
            end
            
            app.Image_plot = subplot(1,2,1,'Parent',app.Figure3Panel);
            set (app.Image_plot, 'position', [0.1    0.06    0.35    0.90]);
            
            imagesc(app.disp_image, 'Parent', app.Image_plot);
            set (app.Image_plot, 'xtick', [], 'xticklabel',[],...
                'ytick',[], 'yticklabel',[]);
            
            colormap (app.Image_plot, 'gray');
            
            if (app.ColorCheckBox.Value)
                colormap (app.Image_plot, 'jet');
            end
            
            L = drawline (app.Image_plot, 'InteractionsAllowed','none');
            hold (app.Image_plot, 'on');
            
            % get the positions of the line
            % L.Position: [x1 y1; x2 y2];
            x_line = [L.Position(1,1), L.Position(2,1)];
            y_line = [L.Position(1,2), L.Position(2,2)];
            
            % finding the distance of the line in mm
            x_dist = (x_line(2) - x_line(1))*app.info.PixelSpacing(1)/10;
            y_dist = (y_line(2) - y_line(1))*app.info.SliceThickness/10;
            dist = sqrt((x_dist)^2 + (y_dist)^2);
            
            app.CT_plot = subplot(1,2,2,'Parent',app.Figure3Panel);
            set (app.CT_plot, 'position', [0.6    0.06    0.35    0.90]);
            disableDefaultInteractivity(app.CT_plot);
            
            % get the profile of the CT number from the line drawn
            [app.cx, app.cy, app.c] = improfile(app.disp_image,x_line, y_line);
            app.scale = dist/length(app.c);
            app.c_depth = (1:length(app.c))*app.scale;
            
            plot (app.CT_plot,app.c, app.c_depth);
            set (app.CT_plot, 'YDir', 'reverse')
            grid (app.CT_plot, 'minor');
            
            hold (app.CT_plot, 'on');
            
            % smoothing the graph
            app.c_smooth = smooth (app.c, 'rlowess');
            plot (app.CT_plot, app.c_smooth, app.c_depth);
            ylabel (app.CT_plot, 'length (cm)');
            xlabel (app.CT_plot, 'CT number');
            
            hold (app.CT_plot, 'on');
            
            app.IdentifypeaksButton.Enable = 'on';
            app.MinDistancecmEditField.Enable = 'on';
            app.MinDistancecmEditFieldLabel.Enable = 'on';
            app.NoofPeaksEditField.Enable = 'on';
            app.NoofPeaksEditFieldLabel.Enable = 'on';
            
            app.SaveFigure3Button.Enable = 'on';
            
            toc
        end

        % Button pushed function: IdentifypeaksButton
        function IdentifypeaksButtonPushed(app, event)
          tic
          
            app.CT_plot = subplot(1,2,2,'Parent',app.Figure3Panel);
            set (app.CT_plot, 'position', [0.6    0.06    0.35    0.90]);
            
            plot (app.CT_plot, app.c, app.c_depth);
            set (app.CT_plot, 'YDir', 'reverse')
            grid (app.CT_plot, 'minor');
            hold (app.CT_plot, 'on');
            plot (app.CT_plot, app.c_smooth, app.c_depth);
            hold (app.CT_plot, 'on');
            
            ylabel (app.CT_plot, 'length (cm)');
            xlabel (app.CT_plot, 'CT number');
            
            % threshold value is the minimum peak distance in pixels
            % converting the threshold in cm to pixel
            p_thresh = app.MinDistancecmEditField.Value/app.scale;
            
            [app.pks, app.locs] = findpeaks (app.c, 'MinPeakDistance', p_thresh);
            c_loc = app.locs*app.scale;
            plot (app.CT_plot,app.pks+30,c_loc, '<b', 'MarkerFaceColor', 'b'); % plotting the peaks
            text(app.CT_plot,app.pks+50,c_loc,num2str((1:numel(app.pks))')); % numbering the peaks
            
            % number of peaks
            peaks_no = length(app.pks);
            app.NoofPeaksEditField.Value = peaks_no;
            
            app.ShowpeaksButton.Enable = 'on';
            
            toc
        end

        % Button pushed function: ShowpeaksButton
        function ShowpeaksButtonPushed(app, event)

           tic 
            %% labelling the peaks on the CT image
            scatter (app.Image_plot, (app.cx(app.locs) + 10) , app.cy (app.locs), '<', 'b', 'filled');
            text(app.Image_plot, app.cx(app.locs)+15, app.cy(app.locs),num2str((1:numel(app.pks))'));
            
            app.CalculatedistanceButton.Enable = 'on';
            
            toc
        end

        % Button pushed function: CalculatedistanceButton
        function CalculatedistanceButtonPushed(app, event)
            %% finding distance between each peak
            % print distance between each peaks
            tic
            peaks_no = length (app.pks);
            for i = 1:(peaks_no-1)
                
                % finding distance between each peak in cm
                dist_diff = (app.locs(i+1) - app.locs (i))*app.scale;
                
                % positions between each peak
                diff_pos = ((app.locs(i+1) + app.locs(i))/2)*app.scale;
                pos2 = ((app.pks(i+1) + app.pks(i))/2);
                text(app.CT_plot, pos2, diff_pos, num2str(dist_diff),...
                    'color', [0.4660, 0.6740, 0.1880]);
                
            end
toc
        end

        % Button pushed function: ExportDICOMButton
        function ExportDICOMButtonPushed(app, event)
          tic
            try
                % choosing the folder to save the cropped volume
                dest = uigetdir;
                dest_folder = app.SaveDICOMfileasEditField.Value;
                mkdir (dest, dest_folder);
                cd([dest filesep dest_folder]);
                
                % the Dicom data is not saved in order so we need to load the
                % dicom files again according to its names so that we can
                % export it with the corresponding metadata
                
                app.order = app.order'; % the order of the name according to instance number
                [~,idx] = sort(app.order(:,1)); % sort according to the name
                
                % creating waitbar
                s = waitbar (0, 'Saving');
                
                % saving cropped volume as dicom files
                if (app.CircleButton.Value)
                    % sort the whole matrix using the sort indices
                    
                    sortedmat = app.cylinder_sorted(:,:,idx);
                    
                    for i = 1:length (idx)
                        BaseName= app.SaveDICOMfileasEditField.Value;
                        filename =[BaseName,num2str(i)];
                        info2 = dicominfo (fullfile (app.filefolder,app.names {i}));
                        dicomwrite (sortedmat(:,:,i),filename, info2);
                        
                        waitbar (i/length(app.names),s);
                    end
                    
                    % close waitbar once all data are saved
                    close (s);
                    
                    msgbox ('Cropped CT images saved','success');
                end
                
                if (app.RectangleButton.Value)
                    
                    sortedmat = app.cube_sorted(:,:,idx);
                    
                    for i = 1:length (idx)
                        BaseName= app.SaveDICOMimagesasEditField.Value;
                        filename =[BaseName,num2str(i)];
                        info2 = dicominfo (fullfile (app.filefolder,app.names {i}));
                        dicomwrite (sortedmat(:,:,i),  filename, info2);
                        
                        waitbar (i/length(app.names),s);
                    end
                    close (s);
                    
                    msgbox ('Cropped CT images saved','success')
                    
                end
                
            catch
                % load a error message if the file is not saved
                errordlg('File not saved','File Error');
            end
            toc
        end

        % Button pushed function: SavedataButton
        function SavedataButtonPushed(app, event)
          tic
            try
                filter = {'*.csv'};
                [filename,filepath] = uiputfile(filter);
                fullfilename = fullfile(filepath, filename);
                
                if ischar(filename)
                    
                    data = horzcat(app.depth, app.all_mean, app.all_sd);
                    
                    col_name = {'Depth (cm)', 'Mean', 'Std'};
                    data = array2table (data,'VariableNames', col_name);
                    
                    writetable(data, fullfilename);
                    
                    msgbox ('File saved','success');
                end
                
            catch
                % load a error message if there are no dicom files found
                errordlg('File not saved','File Error');
            end
            toc
        end

        % Button pushed function: GetpositionsButton
        function GetpositionsButtonPushed(app, event)
        tic
            try
                filter = {'*.csv'};
                [filename,filepath] = uiputfile(filter);
                fullfilename = fullfile(filepath, filename);
                
                if ischar(filename)
                    
                    data = horzcat(app.pts_depth);
                    
                    col_name = {'Position (cm)'};
                    data = array2table (data,'VariableNames', col_name);
                    
                    writetable(data, fullfilename);
                    
                    msgbox ('File saved','success');
                end
                
            catch
                % load a error message if there are no dicom files found
                errordlg('File not saved','File Error');
            end
            toc
        end

        % Button pushed function: SaveFigure1Button
        function SaveFigure1ButtonPushed(app, event)
            try
                
                filter = {'*.jpg';'*.png';'*.tif';'*.pdf';'*.eps';'*.fig'};
                [filename,filepath] = uiputfile(filter);
                
                fig = figure;
                set (gcf, 'position', [440 300 400 600], 'resize', 'off');
                fig.Visible = 'off'; % hide this figure
                
                % copying the images from the app to new figure
                Axial_copy = copyobj ([app.Axial_plot], fig);
                subplot (3,2,[1,2], Axial_copy);
                ori_size = get (gca, 'position');
                c = colorbar; % show the colorbar
                c.Label.String = 'Grayscale value';
                set (gca, 'position', ori_size);
                
                Cor_copy = copyobj (app.Cor_plot, fig);
                subplot (3,2,[3,5], Cor_copy);
                
                Sag_copy = copyobj (app.Sag_plot, fig);
                subplot (3,2,[4,6], Sag_copy);
                
                if ischar(filename) % ensure that there is a filename input
                    
                    % save the figure
                    exportgraphics(fig,[filepath filename]);
                    
                    % inform user that figure has been saved
                    msgbox ('Figure 1 saved','success');
                end
                
                delete(fig) % delete the new figure
                
            catch
                % load a error message if there are no dicom files found
                errordlg('File not saved','File Error');
            end
        end

        % Button pushed function: SaveFigure3Button
        function SaveFigure3ButtonPushed(app, event)
            try
                filter = {'*.png';'*.jpg';'*.tif';'*.pdf';'*.eps';'*.fig'};
                [filename,filepath] = uiputfile(filter);
                
                fig2 = figure;
                set (gcf, 'position', [440 300 600 800], 'resize', 'off');
                fig2.Visible = 'off'; % hide this figure
                
                Img_copy = copyobj ([app.Image_plot], fig2);
                subplot (1,2,1, Img_copy);
                
                CT_copy = copyobj (app.CT_plot, fig2);
                subplot (1,2,2, CT_copy);
                
                if ischar(filename) % ensure that there is a filename input
                    
                    % save the figure
                    exportgraphics(fig2,[filepath filename]);
                    
                    % inform user that figure has been saved
                    msgbox ('Figure 2 saved','success');
                end
                
                delete(fig2) % delete the new figure
                
            catch
                % load a error message if there are no dicom files found
                errordlg('File not saved','File Error');
            end
        end

        % Button pushed function: ResetButton
        function ResetButtonPushed(app, event)
            % clear all the values in the edit field
            
            app.SlicenumberSpinner.Limits = [-inf inf];
            app.AxialsliceSpinner.Limits = [-inf inf];
            
            app.FilenameEditField.Value = '';
            app.XangleSlider.Value = 0;
            app.YangleSlider.Value = 0;
            app.RadiuscmEditField.Value = 0;
            app.WidthcmEditField.Value = 0;
            app.HeightcmEditField.Value = 0;
            app.ThresholdEditField.Value = 0;
            app.NumberofcontactsEditField.Value = 0;
            app.MinDistancecmEditField.Value = 0;
            app.NoofPeaksEditField.Value = 0;
            app.SaveDICOMfileasEditField.Value = '';
            app.TopmmEditField.Value = 0;
            app.EndmmEditField.Value = 0;
            app.SlicenumberSpinner.Value = 0;
            app.ColorCheckBox.Value = 0;
            app.CircleButton.Value = 1;
            app.RectangleButton.Value = 0;
            app.CoronalButton.Value = 1;
            app.SagitalButton.Value = 0;
            app.VoidareaEditField.Value = 0;
            app.AxialsliceSpinner.Value = 0;
            app.TopsliceEditField.Value = 0;
            app.BottomsliceEditField.Value = 0;
            app.VoidvolumeEditField.Value = 0;
            
            delete (app.CT_plot);
            delete (app.Image_plot);
            delete (app.Dis_plot);
            delete (app.area_plot);
            delete (app.Axial_plot);
            delete (app.Cor_plot);
            delete (app.Sag_plot);
            
            % disabling all the buttons
            app.CircleButton.Enable = 'off';
            app.RectangleButton.Enable = 'off';
            app.DrawmaskButton.Enable = 'off';
            app.RadiuscmEditField.Enable = 'off';
            app.RadiuscmEditFieldLabel.Enable = 'off';
            app.HeightcmEditField.Enable = 'off';
            app.HeightcmEditFieldLabel.Enable = 'off';
            app.WidthcmEditField.Enable = 'off';
            app.WidthcmEditFieldLabel.Enable = 'off';
            app.CropButton.Enable = 'off';
            app.ViewallslicesButton.Enable = 'off';
            app.XangleSlider.Enable = 'off';
            app.XangleSliderLabel.Enable = 'off';
            app.YangleSlider.Enable = 'off';
            app.YangleSliderLabel.Enable = 'off';
            app.PlotButton.Enable = 'off';
            app.TopmmEditField.Enable = 'off';
            app.TopmmEditFieldLabel.Enable = 'off';
            app.EndmmEditField.Enable = 'off';
            app.EndmmEditFieldLabel.Enable = 'off';
            app.SaveDICOMfileasEditField.Enable = 'off';
            app.SaveDICOMfileasEditFieldLabel.Enable = 'off';
            app.ExportDICOMButton.Enable = 'off';
            app.SavedataButton.Enable = 'off';
            app.SaveFigure1Button.Enable = 'off';
            app.SaveFigure3Button.Enable = 'off';
            app.GetpositionsButton.Enable = 'off';
            app.DetectButton.Enable = 'off';
            app.NumberofcontactsEditField.Enable = 'off';
            app.NumberofcontactsEditFieldLabel.Enable = 'off';
            app.ThresholdEditField.Enable = 'off';
            app.ThresholdEditFieldLabel.Enable = 'off';
            app.NoofPeaksEditField.Enable = 'off';
            app.NoofPeaksEditFieldLabel.Enable = 'off';
            app.IdentifypeaksButton.Enable = 'off';
            app.MinDistancecmEditField.Enable = 'off';
            app.MinDistancecmEditFieldLabel.Enable = 'off';
            app.CalculatedistanceButton.Enable = 'off';
            app.ShowpeaksButton.Enable = 'off';
            app.CoronalButton.Enable = 'off';
            app.SagitalButton.Enable = 'off';
            app.DrawlineButton.Enable = 'off';
            app.SlicenumberSpinner.Enable = 'off';
            app.SlicenumberSpinnerLabel.Enable = 'off';
            app.ColorCheckBox.Enable = 'off';
            app.VoidareaEditField.Enable = 'off';
            app.VoidareaEditFieldLabel.Enable = 'off';
            app.IdentifyareaButton.Enable = 'off';
            app.AxialsliceSpinner.Enable = 'off';
            app.AxialsliceSpinnerLabel.Enable = 'off';
            app.TopsliceEditField.Enable = 'off';
            app.TopsliceEditFieldLabel.Enable = 'off';
            app.BottomsliceEditField.Enable = 'off';
            app.BottomsliceEditFieldLabel.Enable = 'off';
            app.CalculateVolumeButton.Enable = 'off';
            app.VoidvolumeEditField.Enable = 'off';
            app.VoidvolumeEditFieldLabel.Enable = 'off';
            app.ViewVolumeButton.Enable = 'off';
            app.VolumeViewerButton.Enable = 'off';
            app.FlipuprightCheckBox.Enable = 'off';
            
            % clearing all variables
            app.info = []; % metadata of the dicom images
            app.names = []; % the names of all the dicom images
            app.filefolder =[]; % file of the dicom folder
            app.sorted = []; % 3D matrix arranged by instance number
            app.order = []; % order of the instance number
            
            app.cbar = []; % colorbar for the Axial image
            app.pages = [];
            app.depth =[]; % for plotting the graphs
            app.dim = []; % size of the matrix
            
            % different views of the core
            app.Axial = []; % cross sectional image
            app.Axial_plot =[]; % plot of the axial image
            app.coronal = [];
            app.Cor_plot =[]; % plot of the coronal image
            app.sagital = [];
            app.Sag_plot = []; % plot of the sagital image
            
            app.circle_mask = []; % mask of the circle created (2D)
            app.Rthres = []; % threshold of radius of the circle mask
            app.c_center =[]; % center of circle
            
            app.rect_mask =[]; % mask of the rectangle created (2D)
            app.r_x = []; % x-coordinate of the left bottom point of rectangle
            app.r_y = []; % y-coordinate of the left bottom point of rectang;e
            app.r_width = []; % width of rectangle
            app.r_height  = []; % height of rectangle
            app.r_center = []; % center of rectangle
            
            app.circle_cor = []; % coronal view of the circle mask
            app.circle_sag = []; % sagital view of the circle mask
            app.rect_cor = []; % coronal view of the rectangle mask
            app.rect_sag =[]; % sagital view of the rectangle mask
            
            app.cor_pos = []; % position of first and last slice of cropped volume in coronal
            app.sag_pos = [];
            
            app.cylinder = []; % 3D circle mask
            app.cylinder_sorted =[]; % 3D matrix of the cropped cylinder
            app.cube = []; % 3D rectangle mask
            app.cube_sorted =[]; % 3D matrix of the cropped cuboid
            
            app.all_mean =[]; % mean of the CT numbers
            app.all_sd =[]; % standard deviation of the CT numbers
            app.ulim =[]; % upper limit (ulim) = mean + std
            app.llim =[]; % lower limit (llim) = mean - std
            app.hist_data = []; % data to plot the distribution plot
            
            app.Image_plot = []; % plot of the display image
            app.CT_plot =[];
            app.Dis_plot = [];
            
            app.pts_depth =[]; % depth of the change position
            
            app.disp_image = []; % display image of the cropped DICOM
            app.c = []; % CT profile of line
            app.cx =[]; % spatial coordinates of the CT number profile
            app.cy =[];
            app.c_depth = []; % depth of CT profile
            app.c_smooth = []; % smoothed CT profile
            app.scale = []; % scale to convert image pixel to mm
            app.pks = []; % peaks identified
            app.locs = []; % position of peaks
            
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [200 200 1295 887];
            app.UIFigure.Name = 'MATLAB App';

            % Create ImportyourdataPanel
            app.ImportyourdataPanel = uipanel(app.UIFigure);
            app.ImportyourdataPanel.Title = '1. Import your data';
            app.ImportyourdataPanel.Position = [10 731 240 147];

            % Create ImportButton
            app.ImportButton = uibutton(app.ImportyourdataPanel, 'push');
            app.ImportButton.ButtonPushedFcn = createCallbackFcn(app, @ImportButtonPushed, true);
            app.ImportButton.Position = [11 94 100 22];
            app.ImportButton.Text = 'Import';

            % Create ResetButton
            app.ResetButton = uibutton(app.ImportyourdataPanel, 'push');
            app.ResetButton.ButtonPushedFcn = createCallbackFcn(app, @ResetButtonPushed, true);
            app.ResetButton.Position = [126 94 100 22];
            app.ResetButton.Text = 'Reset';

            % Create FilenameEditFieldLabel
            app.FilenameEditFieldLabel = uilabel(app.ImportyourdataPanel);
            app.FilenameEditFieldLabel.HorizontalAlignment = 'right';
            app.FilenameEditFieldLabel.Position = [10 64 57 22];
            app.FilenameEditFieldLabel.Text = 'File name';

            % Create FilenameEditField
            app.FilenameEditField = uieditfield(app.ImportyourdataPanel, 'text');
            app.FilenameEditField.Editable = 'off';
            app.FilenameEditField.Position = [12 43 213 22];

            % Create FlipuprightCheckBox
            app.FlipuprightCheckBox = uicheckbox(app.ImportyourdataPanel);
            app.FlipuprightCheckBox.ValueChangedFcn = createCallbackFcn(app, @FlipuprightCheckBoxValueChanged, true);
            app.FlipuprightCheckBox.Enable = 'off';
            app.FlipuprightCheckBox.Text = 'Flip upright';
            app.FlipuprightCheckBox.Position = [73 9 82 22];

            % Create SelecttypeofmaskPanel
            app.SelecttypeofmaskPanel = uipanel(app.UIFigure);
            app.SelecttypeofmaskPanel.Title = '2. Select type of mask';
            app.SelecttypeofmaskPanel.Position = [10 560 240 165];

            % Create ButtonGroup
            app.ButtonGroup = uibuttongroup(app.SelecttypeofmaskPanel);
            app.ButtonGroup.SelectionChangedFcn = createCallbackFcn(app, @ButtonGroupSelectionChanged, true);
            app.ButtonGroup.BorderType = 'none';
            app.ButtonGroup.Position = [10 61 100 81];

            % Create CircleButton
            app.CircleButton = uiradiobutton(app.ButtonGroup);
            app.CircleButton.Enable = 'off';
            app.CircleButton.Text = 'Circle';
            app.CircleButton.Position = [3 55 53 22];
            app.CircleButton.Value = true;

            % Create RectangleButton
            app.RectangleButton = uiradiobutton(app.ButtonGroup);
            app.RectangleButton.Enable = 'off';
            app.RectangleButton.Text = 'Rectangle';
            app.RectangleButton.Position = [3 18 76 22];

            % Create DrawmaskButton
            app.DrawmaskButton = uibutton(app.SelecttypeofmaskPanel, 'push');
            app.DrawmaskButton.ButtonPushedFcn = createCallbackFcn(app, @DrawmaskButtonPushed, true);
            app.DrawmaskButton.Enable = 'off';
            app.DrawmaskButton.Position = [64 15 100 22];
            app.DrawmaskButton.Text = 'Draw mask';

            % Create Panel
            app.Panel = uipanel(app.SelecttypeofmaskPanel);
            app.Panel.BorderType = 'none';
            app.Panel.Position = [87 107 142 30];

            % Create RadiuscmEditFieldLabel
            app.RadiuscmEditFieldLabel = uilabel(app.Panel);
            app.RadiuscmEditFieldLabel.HorizontalAlignment = 'right';
            app.RadiuscmEditFieldLabel.Enable = 'off';
            app.RadiuscmEditFieldLabel.Position = [7 9 69 22];
            app.RadiuscmEditFieldLabel.Text = 'Radius (cm)';

            % Create RadiuscmEditField
            app.RadiuscmEditField = uieditfield(app.Panel, 'numeric');
            app.RadiuscmEditField.Editable = 'off';
            app.RadiuscmEditField.Enable = 'off';
            app.RadiuscmEditField.Position = [91 9 52 22];

            % Create Panel_2
            app.Panel_2 = uipanel(app.SelecttypeofmaskPanel);
            app.Panel_2.BorderType = 'none';
            app.Panel_2.Visible = 'off';
            app.Panel_2.Position = [98 46 138 62];

            % Create HeightcmEditFieldLabel
            app.HeightcmEditFieldLabel = uilabel(app.Panel_2);
            app.HeightcmEditFieldLabel.HorizontalAlignment = 'right';
            app.HeightcmEditFieldLabel.Position = [1 33 67 22];
            app.HeightcmEditFieldLabel.Text = 'Height (cm)';

            % Create HeightcmEditField
            app.HeightcmEditField = uieditfield(app.Panel_2, 'numeric');
            app.HeightcmEditField.Editable = 'off';
            app.HeightcmEditField.Position = [82 33 50 22];

            % Create WidthcmEditFieldLabel
            app.WidthcmEditFieldLabel = uilabel(app.Panel_2);
            app.WidthcmEditFieldLabel.HorizontalAlignment = 'right';
            app.WidthcmEditFieldLabel.Position = [4 1 63 22];
            app.WidthcmEditFieldLabel.Text = 'Width (cm)';

            % Create WidthcmEditField
            app.WidthcmEditField = uieditfield(app.Panel_2, 'numeric');
            app.WidthcmEditField.Editable = 'off';
            app.WidthcmEditField.Position = [81 1 51 22];

            % Create CropthedataPanel
            app.CropthedataPanel = uipanel(app.UIFigure);
            app.CropthedataPanel.Title = '4. Crop the data';
            app.CropthedataPanel.Position = [10 258 240 93];

            % Create CropButton
            app.CropButton = uibutton(app.CropthedataPanel, 'push');
            app.CropButton.ButtonPushedFcn = createCallbackFcn(app, @CropButtonPushed, true);
            app.CropButton.Enable = 'off';
            app.CropButton.Position = [63 40 100 22];
            app.CropButton.Text = 'Crop';

            % Create ViewallslicesButton
            app.ViewallslicesButton = uibutton(app.CropthedataPanel, 'push');
            app.ViewallslicesButton.ButtonPushedFcn = createCallbackFcn(app, @ViewallslicesButtonPushed, true);
            app.ViewallslicesButton.Enable = 'off';
            app.ViewallslicesButton.Position = [13 9 100 22];
            app.ViewallslicesButton.Text = 'View all slices';

            % Create VolumeViewerButton
            app.VolumeViewerButton = uibutton(app.CropthedataPanel, 'push');
            app.VolumeViewerButton.ButtonPushedFcn = createCallbackFcn(app, @VolumeViewerButtonPushed, true);
            app.VolumeViewerButton.Enable = 'off';
            app.VolumeViewerButton.Position = [127 9 98 23];
            app.VolumeViewerButton.Text = 'Volume Viewer';

            % Create AdjustifneededPanel
            app.AdjustifneededPanel = uipanel(app.UIFigure);
            app.AdjustifneededPanel.Title = '3. Adjust if needed';
            app.AdjustifneededPanel.Position = [10 359 240 190];

            % Create XangleSliderLabel
            app.XangleSliderLabel = uilabel(app.AdjustifneededPanel);
            app.XangleSliderLabel.HorizontalAlignment = 'right';
            app.XangleSliderLabel.Enable = 'off';
            app.XangleSliderLabel.Position = [18 137 50 22];
            app.XangleSliderLabel.Text = 'X angle°';

            % Create XangleSlider
            app.XangleSlider = uislider(app.AdjustifneededPanel);
            app.XangleSlider.Limits = [-10 10];
            app.XangleSlider.ValueChangingFcn = createCallbackFcn(app, @XangleSliderValueChanging, true);
            app.XangleSlider.Enable = 'off';
            app.XangleSlider.Position = [14 130 212 3];

            % Create YangleSliderLabel
            app.YangleSliderLabel = uilabel(app.AdjustifneededPanel);
            app.YangleSliderLabel.HorizontalAlignment = 'right';
            app.YangleSliderLabel.Enable = 'off';
            app.YangleSliderLabel.Position = [15 57 50 22];
            app.YangleSliderLabel.Text = 'Y angle°';

            % Create YangleSlider
            app.YangleSlider = uislider(app.AdjustifneededPanel);
            app.YangleSlider.Limits = [-10 10];
            app.YangleSlider.ValueChangingFcn = createCallbackFcn(app, @YangleSliderValueChanging, true);
            app.YangleSlider.Enable = 'off';
            app.YangleSlider.Position = [14 49 213 3];

            % Create PlotthedataPanel
            app.PlotthedataPanel = uipanel(app.UIFigure);
            app.PlotthedataPanel.Title = '5. Plot the data';
            app.PlotthedataPanel.Position = [10 137 240 115];

            % Create PlotButton
            app.PlotButton = uibutton(app.PlotthedataPanel, 'push');
            app.PlotButton.ButtonPushedFcn = createCallbackFcn(app, @PlotButtonPushed, true);
            app.PlotButton.Enable = 'off';
            app.PlotButton.Position = [59 9 100 22];
            app.PlotButton.Text = 'Plot';

            % Create TopmmEditFieldLabel
            app.TopmmEditFieldLabel = uilabel(app.PlotthedataPanel);
            app.TopmmEditFieldLabel.HorizontalAlignment = 'right';
            app.TopmmEditFieldLabel.Enable = 'off';
            app.TopmmEditFieldLabel.Position = [30 66 56 22];
            app.TopmmEditFieldLabel.Text = 'Top (mm)';

            % Create TopmmEditField
            app.TopmmEditField = uieditfield(app.PlotthedataPanel, 'numeric');
            app.TopmmEditField.Enable = 'off';
            app.TopmmEditField.Position = [36 43 52 22];

            % Create EndmmEditFieldLabel
            app.EndmmEditFieldLabel = uilabel(app.PlotthedataPanel);
            app.EndmmEditFieldLabel.HorizontalAlignment = 'right';
            app.EndmmEditFieldLabel.Enable = 'off';
            app.EndmmEditFieldLabel.Position = [137 66 57 22];
            app.EndmmEditFieldLabel.Text = 'End (mm)';

            % Create EndmmEditField
            app.EndmmEditField = uieditfield(app.PlotthedataPanel, 'numeric');
            app.EndmmEditField.Enable = 'off';
            app.EndmmEditField.Position = [137 43 50 22];

            % Create ExportPanel
            app.ExportPanel = uipanel(app.UIFigure);
            app.ExportPanel.Title = '9. Export';
            app.ExportPanel.Position = [1052 9 234 204];

            % Create SaveDICOMfileasEditFieldLabel
            app.SaveDICOMfileasEditFieldLabel = uilabel(app.ExportPanel);
            app.SaveDICOMfileasEditFieldLabel.HorizontalAlignment = 'right';
            app.SaveDICOMfileasEditFieldLabel.Enable = 'off';
            app.SaveDICOMfileasEditFieldLabel.Position = [14 149 110 22];
            app.SaveDICOMfileasEditFieldLabel.Text = 'Save DICOM file as';

            % Create SaveDICOMfileasEditField
            app.SaveDICOMfileasEditField = uieditfield(app.ExportPanel, 'text');
            app.SaveDICOMfileasEditField.Enable = 'off';
            app.SaveDICOMfileasEditField.Position = [14 128 201 22];

            % Create ExportDICOMButton
            app.ExportDICOMButton = uibutton(app.ExportPanel, 'push');
            app.ExportDICOMButton.ButtonPushedFcn = createCallbackFcn(app, @ExportDICOMButtonPushed, true);
            app.ExportDICOMButton.Enable = 'off';
            app.ExportDICOMButton.Position = [63 93 100 22];
            app.ExportDICOMButton.Text = 'Export DICOM';

            % Create SavedataButton
            app.SavedataButton = uibutton(app.ExportPanel, 'push');
            app.SavedataButton.ButtonPushedFcn = createCallbackFcn(app, @SavedataButtonPushed, true);
            app.SavedataButton.Enable = 'off';
            app.SavedataButton.Position = [13 53 100 22];
            app.SavedataButton.Text = 'Save data';

            % Create SaveFigure1Button
            app.SaveFigure1Button = uibutton(app.ExportPanel, 'push');
            app.SaveFigure1Button.ButtonPushedFcn = createCallbackFcn(app, @SaveFigure1ButtonPushed, true);
            app.SaveFigure1Button.Enable = 'off';
            app.SaveFigure1Button.Position = [13 9 100 22];
            app.SaveFigure1Button.Text = 'Save Figure 1';

            % Create SaveFigure3Button
            app.SaveFigure3Button = uibutton(app.ExportPanel, 'push');
            app.SaveFigure3Button.ButtonPushedFcn = createCallbackFcn(app, @SaveFigure3ButtonPushed, true);
            app.SaveFigure3Button.Enable = 'off';
            app.SaveFigure3Button.Position = [120 9 100 22];
            app.SaveFigure3Button.Text = 'Save Figure 3';

            % Create GetpositionsButton
            app.GetpositionsButton = uibutton(app.ExportPanel, 'push');
            app.GetpositionsButton.ButtonPushedFcn = createCallbackFcn(app, @GetpositionsButtonPushed, true);
            app.GetpositionsButton.Enable = 'off';
            app.GetpositionsButton.Position = [120 53 100 22];
            app.GetpositionsButton.Text = 'Get positions';

            % Create Figure1Panel
            app.Figure1Panel = uipanel(app.UIFigure);
            app.Figure1Panel.Title = 'Figure 1';
            app.Figure1Panel.Position = [259 9 442 869];

            % Create Figure3Panel
            app.Figure3Panel = uipanel(app.UIFigure);
            app.Figure3Panel.Title = 'Figure 3';
            app.Figure3Panel.Position = [712 9 331 583];

            % Create AnalysisPanel
            app.AnalysisPanel = uipanel(app.UIFigure);
            app.AnalysisPanel.Title = '6. Analysis';
            app.AnalysisPanel.Position = [10 9 240 120];

            % Create DetectButton
            app.DetectButton = uibutton(app.AnalysisPanel, 'push');
            app.DetectButton.ButtonPushedFcn = createCallbackFcn(app, @DetectButtonPushed, true);
            app.DetectButton.Enable = 'off';
            app.DetectButton.Position = [69 43 100 22];
            app.DetectButton.Text = 'Detect';

            % Create NumberofcontactsEditFieldLabel
            app.NumberofcontactsEditFieldLabel = uilabel(app.AnalysisPanel);
            app.NumberofcontactsEditFieldLabel.HorizontalAlignment = 'right';
            app.NumberofcontactsEditFieldLabel.Enable = 'off';
            app.NumberofcontactsEditFieldLabel.Position = [30 10 112 22];
            app.NumberofcontactsEditFieldLabel.Text = 'Number of contacts';

            % Create NumberofcontactsEditField
            app.NumberofcontactsEditField = uieditfield(app.AnalysisPanel, 'numeric');
            app.NumberofcontactsEditField.Editable = 'off';
            app.NumberofcontactsEditField.Enable = 'off';
            app.NumberofcontactsEditField.Position = [151 10 42 22];

            % Create ThresholdEditFieldLabel
            app.ThresholdEditFieldLabel = uilabel(app.AnalysisPanel);
            app.ThresholdEditFieldLabel.HorizontalAlignment = 'right';
            app.ThresholdEditFieldLabel.Enable = 'off';
            app.ThresholdEditFieldLabel.Position = [9 71 59 22];
            app.ThresholdEditFieldLabel.Text = 'Threshold';

            % Create ThresholdEditField
            app.ThresholdEditField = uieditfield(app.AnalysisPanel, 'numeric');
            app.ThresholdEditField.Enable = 'off';
            app.ThresholdEditField.Position = [83 71 51 22];

            % Create IdentifyingpeaksPanel
            app.IdentifyingpeaksPanel = uipanel(app.UIFigure);
            app.IdentifyingpeaksPanel.Title = '8. Identifying peaks';
            app.IdentifyingpeaksPanel.Position = [1052 220 234 214];

            % Create NoofPeaksEditFieldLabel
            app.NoofPeaksEditFieldLabel = uilabel(app.IdentifyingpeaksPanel);
            app.NoofPeaksEditFieldLabel.HorizontalAlignment = 'right';
            app.NoofPeaksEditFieldLabel.Enable = 'off';
            app.NoofPeaksEditFieldLabel.Position = [44 55 74 22];
            app.NoofPeaksEditFieldLabel.Text = 'No. of Peaks';

            % Create NoofPeaksEditField
            app.NoofPeaksEditField = uieditfield(app.IdentifyingpeaksPanel, 'numeric');
            app.NoofPeaksEditField.Editable = 'off';
            app.NoofPeaksEditField.Enable = 'off';
            app.NoofPeaksEditField.Position = [133 59 49 22];

            % Create IdentifypeaksButton
            app.IdentifypeaksButton = uibutton(app.IdentifyingpeaksPanel, 'push');
            app.IdentifypeaksButton.ButtonPushedFcn = createCallbackFcn(app, @IdentifypeaksButtonPushed, true);
            app.IdentifypeaksButton.Enable = 'off';
            app.IdentifypeaksButton.Position = [13 103 100 22];
            app.IdentifypeaksButton.Text = 'Identify peaks';

            % Create MinDistancecmEditFieldLabel
            app.MinDistancecmEditFieldLabel = uilabel(app.IdentifyingpeaksPanel);
            app.MinDistancecmEditFieldLabel.HorizontalAlignment = 'right';
            app.MinDistancecmEditFieldLabel.Enable = 'off';
            app.MinDistancecmEditFieldLabel.Position = [20 148 105 22];
            app.MinDistancecmEditFieldLabel.Text = 'Min. Distance (cm)';

            % Create MinDistancecmEditField
            app.MinDistancecmEditField = uieditfield(app.IdentifyingpeaksPanel, 'numeric');
            app.MinDistancecmEditField.Enable = 'off';
            app.MinDistancecmEditField.Position = [135 148 63 22];

            % Create CalculatedistanceButton
            app.CalculatedistanceButton = uibutton(app.IdentifyingpeaksPanel, 'push');
            app.CalculatedistanceButton.ButtonPushedFcn = createCallbackFcn(app, @CalculatedistanceButtonPushed, true);
            app.CalculatedistanceButton.Enable = 'off';
            app.CalculatedistanceButton.Position = [55 17 115 22];
            app.CalculatedistanceButton.Text = 'Calculate distance';

            % Create ShowpeaksButton
            app.ShowpeaksButton = uibutton(app.IdentifyingpeaksPanel, 'push');
            app.ShowpeaksButton.ButtonPushedFcn = createCallbackFcn(app, @ShowpeaksButtonPushed, true);
            app.ShowpeaksButton.Enable = 'off';
            app.ShowpeaksButton.Position = [117 103 100 22];
            app.ShowpeaksButton.Text = 'Show peaks';

            % Create SelectplanePanel
            app.SelectplanePanel = uipanel(app.UIFigure);
            app.SelectplanePanel.Title = '7. Select plane';
            app.SelectplanePanel.Position = [1052 441 234 151];

            % Create ButtonGroup_2
            app.ButtonGroup_2 = uibuttongroup(app.SelectplanePanel);
            app.ButtonGroup_2.SelectionChangedFcn = createCallbackFcn(app, @ButtonGroup_2SelectionChanged, true);
            app.ButtonGroup_2.BorderType = 'none';
            app.ButtonGroup_2.Position = [21 75 100 53];

            % Create CoronalButton
            app.CoronalButton = uiradiobutton(app.ButtonGroup_2);
            app.CoronalButton.Enable = 'off';
            app.CoronalButton.Text = 'Coronal';
            app.CoronalButton.Position = [13 24 64 22];
            app.CoronalButton.Value = true;

            % Create SagitalButton
            app.SagitalButton = uiradiobutton(app.ButtonGroup_2);
            app.SagitalButton.Enable = 'off';
            app.SagitalButton.Text = 'Sagital';
            app.SagitalButton.Position = [13 3 59 22];

            % Create DrawlineButton
            app.DrawlineButton = uibutton(app.SelectplanePanel, 'push');
            app.DrawlineButton.ButtonPushedFcn = createCallbackFcn(app, @DrawlineButtonPushed, true);
            app.DrawlineButton.Enable = 'off';
            app.DrawlineButton.Position = [67 11 100 22];
            app.DrawlineButton.Text = 'Draw line';

            % Create SlicenumberSpinnerLabel
            app.SlicenumberSpinnerLabel = uilabel(app.SelectplanePanel);
            app.SlicenumberSpinnerLabel.HorizontalAlignment = 'right';
            app.SlicenumberSpinnerLabel.Enable = 'off';
            app.SlicenumberSpinnerLabel.Position = [40 45 76 22];
            app.SlicenumberSpinnerLabel.Text = 'Slice number';

            % Create SlicenumberSpinner
            app.SlicenumberSpinner = uispinner(app.SelectplanePanel);
            app.SlicenumberSpinner.Limits = [0 Inf];
            app.SlicenumberSpinner.ValueChangedFcn = createCallbackFcn(app, @SlicenumberSpinnerValueChanged, true);
            app.SlicenumberSpinner.Enable = 'off';
            app.SlicenumberSpinner.Position = [128 45 70 22];

            % Create ColorCheckBox
            app.ColorCheckBox = uicheckbox(app.SelectplanePanel);
            app.ColorCheckBox.ValueChangedFcn = createCallbackFcn(app, @ColorCheckBoxValueChanged, true);
            app.ColorCheckBox.Enable = 'off';
            app.ColorCheckBox.Text = 'Color';
            app.ColorCheckBox.Position = [130 86 51 22];

            % Create Figure2Panel_2
            app.Figure2Panel_2 = uipanel(app.UIFigure);
            app.Figure2Panel_2.Title = 'Figure 2';
            app.Figure2Panel_2.Position = [712 605 574 273];

            % Create VoidareaEditFieldLabel
            app.VoidareaEditFieldLabel = uilabel(app.Figure2Panel_2);
            app.VoidareaEditFieldLabel.HorizontalAlignment = 'right';
            app.VoidareaEditFieldLabel.Enable = 'off';
            app.VoidareaEditFieldLabel.Position = [448 164 58 22];
            app.VoidareaEditFieldLabel.Text = 'Void area ';

            % Create VoidareaEditField
            app.VoidareaEditField = uieditfield(app.Figure2Panel_2, 'numeric');
            app.VoidareaEditField.Enable = 'off';
            app.VoidareaEditField.Position = [505 164 59 22];

            % Create AxialsliceSpinnerLabel
            app.AxialsliceSpinnerLabel = uilabel(app.Figure2Panel_2);
            app.AxialsliceSpinnerLabel.HorizontalAlignment = 'right';
            app.AxialsliceSpinnerLabel.Enable = 'off';
            app.AxialsliceSpinnerLabel.Position = [433 226 60 22];
            app.AxialsliceSpinnerLabel.Text = 'Axial slice';

            % Create AxialsliceSpinner
            app.AxialsliceSpinner = uispinner(app.Figure2Panel_2);
            app.AxialsliceSpinner.ValueChangedFcn = createCallbackFcn(app, @AxialsliceSpinnerValueChanged, true);
            app.AxialsliceSpinner.Enable = 'off';
            app.AxialsliceSpinner.Position = [498 226 68 22];

            % Create IdentifyareaButton
            app.IdentifyareaButton = uibutton(app.Figure2Panel_2, 'push');
            app.IdentifyareaButton.ButtonPushedFcn = createCallbackFcn(app, @IdentifyareaButtonPushed, true);
            app.IdentifyareaButton.Enable = 'off';
            app.IdentifyareaButton.Position = [454 196 100 22];
            app.IdentifyareaButton.Text = 'Identify area';

            % Create CalculateVolumeButton
            app.CalculateVolumeButton = uibutton(app.Figure2Panel_2, 'push');
            app.CalculateVolumeButton.ButtonPushedFcn = createCallbackFcn(app, @CalculateVolumeButtonPushed, true);
            app.CalculateVolumeButton.Enable = 'off';
            app.CalculateVolumeButton.Position = [453 68 108 22];
            app.CalculateVolumeButton.Text = 'Calculate Volume';

            % Create VoidvolumeEditFieldLabel
            app.VoidvolumeEditFieldLabel = uilabel(app.Figure2Panel_2);
            app.VoidvolumeEditFieldLabel.HorizontalAlignment = 'right';
            app.VoidvolumeEditFieldLabel.Enable = 'off';
            app.VoidvolumeEditFieldLabel.Position = [443 37 71 22];
            app.VoidvolumeEditFieldLabel.Text = 'Void volume';

            % Create VoidvolumeEditField
            app.VoidvolumeEditField = uieditfield(app.Figure2Panel_2, 'numeric');
            app.VoidvolumeEditField.Editable = 'off';
            app.VoidvolumeEditField.Enable = 'off';
            app.VoidvolumeEditField.Position = [521 37 43 22];

            % Create TopsliceEditFieldLabel
            app.TopsliceEditFieldLabel = uilabel(app.Figure2Panel_2);
            app.TopsliceEditFieldLabel.HorizontalAlignment = 'right';
            app.TopsliceEditFieldLabel.Enable = 'off';
            app.TopsliceEditFieldLabel.Position = [460 126 53 22];
            app.TopsliceEditFieldLabel.Text = 'Top slice';

            % Create TopsliceEditField
            app.TopsliceEditField = uieditfield(app.Figure2Panel_2, 'numeric');
            app.TopsliceEditField.Enable = 'off';
            app.TopsliceEditField.Position = [518 126 46 22];

            % Create BottomsliceEditFieldLabel
            app.BottomsliceEditFieldLabel = uilabel(app.Figure2Panel_2);
            app.BottomsliceEditFieldLabel.HorizontalAlignment = 'right';
            app.BottomsliceEditFieldLabel.Enable = 'off';
            app.BottomsliceEditFieldLabel.Position = [442 100 72 22];
            app.BottomsliceEditFieldLabel.Text = 'Bottom slice';

            % Create BottomsliceEditField
            app.BottomsliceEditField = uieditfield(app.Figure2Panel_2, 'numeric');
            app.BottomsliceEditField.Enable = 'off';
            app.BottomsliceEditField.Position = [518 100 46 22];

            % Create ViewVolumeButton
            app.ViewVolumeButton = uibutton(app.Figure2Panel_2, 'push');
            app.ViewVolumeButton.ButtonPushedFcn = createCallbackFcn(app, @ViewVolumeButtonPushed, true);
            app.ViewVolumeButton.Enable = 'off';
            app.ViewVolumeButton.Position = [461 6 100 22];
            app.ViewVolumeButton.Text = 'View Volume';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = CoreCT_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end