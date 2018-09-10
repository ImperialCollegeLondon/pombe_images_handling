classdef ExploThreshContours

    % same as Exploration except we can extract contours from thresholding

    methods(Static)

        % starting
        function start(fov,output_file)
            % initialize view state and main display
            state.fov = fov;
            state.z = 1; state.t = 1; state.ch = 1; state.rect = [1 1 fov.Width fov.Height];
            im = fov.loadImage(1,state.z,fov.Ch_short_names{state.ch});
            state.h_im = Operations.imshowfit( Operations.normalize(im,'normalize_contrast') ); hold on;
            state.threshold = 0.5; state.threshold_dir = 1;
            text_label = strjoin(['ch = ' state.fov.Ch_short_names(state.ch) '/ z = ' num2str(state.z) '/ t = ' num2str(state.t)] );
            state.h_label = text('String',text_label,'Position',[0.6*state.rect(3) 0.05*state.rect(4)],'Color','w','FontSize',20);
            % for output
            state.output_file = output_file;
            % to store the segments
            if exist(state.output_file,'file') %load if exist already
                contours = load(state.output_file);
                state.contours = contours.contours;
            else
                state.contours = struct('pts',{},'rect',{},'t',{},'z',{},'ch',{});
            end
            state.h_contours = [];
            setappdata(gcf,'state',state);
            % update once to show existing rectangles
            ExploThreshContours.update();
            % bind display to callback
            set(gcf,'KeyPressFcn',@(obj,event)ExploThreshContours.callback(obj,event.Key));
            waitfor(gcf);
        end

        % callback entry
        function callback(~,key)
            % first check for basic change of view (use Exploration code)
            if Exploration.callbackView(key)
                ExploThreshContours.update();
                return;
            end
            % check for rectangle click
            if ExploThreshContours.callbackPointClick(key)
                return;
            end
            % quitting
            if strcmp(key,'space')
                close; return;
            end
        end

        % callback check: rectangle click
        function valid_key = callbackPointClick(key)
            valid_key = false;
            state = getappdata(gcf,'state');
            % click new rect
            if strcmp(key,'x')
                valid_key = true;
                % compute the masks of objects
                bw = bwareaopen(imfill(imclearborder(state.im > state.threshold),'holes'),30); % hardcoded min pixels, to improve
                boundaries = bwboundaries(bw);
                % add each of them
                for i_cc = 1:length(boundaries)
                    boundary = boundaries{i_cc};
                    boundary(:,2) = boundary(:,2) + state.rect(1) - 1;
                    boundary(:,1) = boundary(:,1) + state.rect(2) - 1;
                    % add contour
                    state.contours(end+1).pts = [boundary(:,2),boundary(:,1)];
                    state.contours(end).rect = state.rect;
                    state.contours(end).t = state.t;
                    state.contours(end).z = state.z;
                    state.contours(end).ch = state.ch;
                end
                setappdata(gcf,'state',state);
                ExploThreshContours.save2disk();
                ExploThreshContours.update();
                return;
            end
            % delete last
            if strcmp(key,'backspace')
                state.contours(end) = [];
                setappdata(gcf,'state',state);
                ExploThreshContours.save2disk();
                ExploThreshContours.update();
                return;
            end
        end

        % update
        function update()
            % delete previously drawn contours
            state = getappdata(gcf,'state');
            for i=1:length(state.h_contours)
                delete(state.h_contours(i));
            end
            state.h_contours = [];
            % draw the contours
            colors = hsv(100);
            for i=1:length(state.contours)
                contour = state.contours(i);
                if (contour.t == state.t)
                    if (contour.z == state.z)
                        state.h_contours(end+1) = plot( contour.pts(:,1) - state.rect(1) + 1 , contour.pts(:,2) - state.rect(2) + 1 , 'Color' , colors(ceil(rand()*100),:) );
                    else
                        state.h_contours(end+1) = plot( contour.pts(:,1) - state.rect(1) + 1 , contour.pts(:,2) - state.rect(2) + 1 , '--' , 'Color' , colors(ceil(rand()*100),:) );
                    end
                end
            end
            setappdata(gcf,'state',state);
        end

        % save to disk
        function save2disk()
            state = getappdata(gcf,'state');
            contours = state.contours;
            fov = state.fov;
            save(state.output_file,'contours','fov');
        end

        % methods to extract data from contour
        function extractCellData(matlab_data_file,output_data_file,output_cells_folder,pixel_size_um,do_plot)
            if nargin < 5
                do_plot = true;
            end
            if do_plot; mkdir(output_cells_folder); end
            contours = load(matlab_data_file);
            fov = contours.fov;
            contours = contours.contours;
            l_um = []; CS_um2 = []; V_um3 = []; w_avg_um = []; SA_um2 = []; n_pixels = []; barycenter_x = []; barycenter_y = [];
            for i_ch = 1:fov.N_ch
                eval(['avg_pixel_intensity_' fov.Ch_short_names{i_ch} ' = [];']);
                eval(['nhood_median_pixel_intensity_' fov.Ch_short_names{i_ch} ' = [];']);
                eval(['nhood_lowest_pixel_intensities_' fov.Ch_short_names{i_ch} ' = [];']);
            end
            for i=1:length(contours)
                disp([num2str(i) '/' num2str(length(contours))]);
                contour = contours(i);
                % geometry
                geom = Operations.contour2geom(contour.pts);
                l_um(end+1,1) = geom.l * pixel_size_um;
                V_um3(end+1,1) = geom.V * pixel_size_um^3;
                CS_um2(end+1,1) = geom.CS * pixel_size_um^2;
                n_pixels(end+1,1) = geom.CS;
                w_avg_um(end+1,1) = geom.width_avg * pixel_size_um;
                SA_um2(end+1,1) = geom.SA * pixel_size_um^2;
                % location of barycenter
                barycenter_x(end+1,1) = mean(contour.pts(:,1));
                barycenter_y(end+1,1) = mean(contour.pts(:,2));
                % channels
                im_to_plot = []; contours_to_plot = {};
                for i_ch = 1:fov.N_ch
                    im = fov.loadImage(contour.t,contour.z,fov.Ch_short_names{i_ch});
                    cell_mask = poly2mask(contour.pts(:,1),contour.pts(:,2),size(im,1),size(im,2));
                    pixel_vals = im(cell_mask);
                    rect_cell = Operations.boundary2rect([contour.pts(:,2),contour.pts(:,1)],size(im,2),size(im,1));
                    % pad by 100 pixels (50 left/right/top/bottom) to get neighbourhood
                    rect_cell = Operations.padRect(rect_cell,100,size(im,2),size(im,1));
                    im_cropped = imcrop(im,rect_cell);
                    region_pixel_vals = sort(im_cropped(:));
                    eval(['avg_pixel_intensity_' fov.Ch_short_names{i_ch} '(end+1,1) = mean(pixel_vals);']);
                    eval(['nhood_median_pixel_intensity_' fov.Ch_short_names{i_ch} '(end+1,1) = median(region_pixel_vals);']);
                    eval(['nhood_lowest_pixel_intensities_' fov.Ch_short_names{i_ch} '(end+1,1) = mean(region_pixel_vals(1:10));']);
                    % for plotting
                    if do_plot
                        im_to_plot = [ im_to_plot , Operations.normalize(im_cropped,'normalize_contrast')];
                        contours_to_plot{end+1,1} = [contour.pts(:,1)-rect_cell(1)+1+(i_ch-1)*size(im_cropped,2),contour.pts(:,2)-rect_cell(2)+1];
                    end
                end
                % plotting
                if do_plot
                    Operations.imshowfit(im_to_plot); hold on;
                    for i_ch = 1:fov.N_ch
                        plot(contours_to_plot{i_ch}(:,1),contours_to_plot{i_ch}(:,2),'r'); hold on;
                    end
                    saveas(gcf,[output_cells_folder '/cell-' num2str(i) '.png']); close;
                end
            end
            table_creation_str = 'data = table(l_um,w_avg_um,CS_um2,SA_um2,V_um3,n_pixels,barycenter_x,barycenter_y';
            for i_ch = 1:fov.N_ch
                table_creation_str = [table_creation_str ',avg_pixel_intensity_' fov.Ch_short_names{i_ch}];
                table_creation_str = [table_creation_str ',nhood_median_pixel_intensity_' fov.Ch_short_names{i_ch}];
                table_creation_str = [table_creation_str ',nhood_lowest_pixel_intensities_' fov.Ch_short_names{i_ch}];
            end
            table_creation_str = [table_creation_str ');'];
            eval(table_creation_str);
            writetable(data,output_data_file);
        end
    end

end
