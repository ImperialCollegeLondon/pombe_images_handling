classdef myFieldOfView
    
    properties
        Ch_names % names of channels (as in image files saved my Micromanager)
        Ch_short_names % user-given shorthand names for the channels (e.g. BF for Brightfield, etc)
        Z_max % number of z slices
        T_max % timepoint number
        Folder_path % where the images are
        N_ch % number of channels
        Width % image width
        Height % image height
    end
    
    methods
        % constructor
        % you need to provide channel names (and shorthand names) and the
        % folder path
        % Z-max and T-max are computed automatically by inspecting files in
        % the folder
        function obj = myFieldOfView(ch_names,ch_short_names,folder_path)
            obj.Ch_names = ch_names;
            obj.Ch_short_names = ch_short_names;
            obj.Folder_path = cell2mat(folder_path);
            obj.N_ch = length(ch_names);
            im = obj.loadImage(1,1,obj.Ch_short_names{1});
            obj.Width = size(im,2);
            obj.Height = size(im,1);
            obj.Z_max = 1;
            obj.T_max = 1;
        end
        
        % load single image
        function im = loadImage(obj,t,z,ch_short_name,normalize_opt)
            if nargin<5
                normalize_opt = 'none';
            end
            for ch=1:obj.N_ch
                if strcmp(ch_short_name,obj.Ch_short_names{ch})
                    im = imread(strcat(obj.Folder_path,obj.Ch_names{ch},'.ome.tif'));
                end
            end
            im = Operations.normalize(im,normalize_opt);
        end
        
        % load single image and crop
        function im = loadCroppedImage(obj,t,z,ch_short_name,rect,normalize_opt)
            if nargin<6
                normalize_opt = 'none';
            end
            im = obj.loadImage(t,z,ch_short_name);
            im = imcrop(im,rect);
            im = Operations.normalize(im,normalize_opt);
        end
        
        % get Z-MAX image and crop
        function im = loadZMaxCroppedImage(obj,t,ch_short_name,rect,normalize_opt)
            im_stack = obj.loadAs3D(t,ch_short_name,rect);
            im = max(im_stack,[],3);
            im = Operations.normalize(im,normalize_opt);
        end
        
        % get Z-MIN image and crop
        function im = loadZMinCroppedImage(obj,t,ch_short_name,rect,normalize_opt)
            im_stack = obj.loadAs3D(t,ch_short_name,rect);
            im = min(im_stack,[],3);
            im = Operations.normalize(im,normalize_opt);
        end
        
        % get Z-MIN x Z-MAX
        function im = loadZMinZMaxProductCroppedImage(obj,t,ch_short_name,rect)
            im_zmin = obj.loadZMinCroppedImage(t,ch_short_name,rect,'normalize_contrast');
            im_zmax = obj.loadZMaxCroppedImage(t,ch_short_name,rect,'normalize_contrast');
            im = immultiply( im_zmax , im_zmin );
            im = Operations.normalize( im , 'normalize_contrast' );
        end
        
        % compute CV along z and crop
        function im = loadZCvCroppedImage(obj,t,ch_short_name,rect)
            im_stack = obj.loadAs3D(t,ch_short_name,rect);
            im = std(im_stack,[],3) ./ mean(im_stack,3);
        end
        
        % side2side ch X z and crop
        function im = loadSide2SideCroppedImage(obj,t,rect)
            im_size = size( imcrop( obj.loadImage(t,1,obj.Ch_short_names{1}) , rect ) );
            im = zeros( obj.N_ch * im_size(1) , obj.Z_max * im_size(2));
            for c=1:obj.N_ch
                for z=1:obj.Z_max
                    im( 1+im_size(1)*(c-1):im_size(1)*c , 1+im_size(2)*(z-1):im_size(2)*z ) =  ...
                        obj.loadCroppedImage(t,z,obj.Ch_short_names{c},rect,'normalize_contrast');
                end
            end
        end
        
        % to 3D data
        function im_stack = loadAs3D(obj,t,ch_short_name,rect)
            im_size = size( imcrop( obj.loadImage(t,1,ch_short_name) , rect ) );
            im_stack=zeros(im_size(1),im_size(2),obj.Z_max);
            for z=1:obj.Z_max
                im = obj.loadImage(t,z,ch_short_name);
                im = imcrop(im,rect);
                im_stack(:,:,z) = im;
            end
        end
        
        % computeBestZ
        function z_best = computeBestZ(obj,t,ch_short_name,rect,score_fun)
            M = -1e300;
            for z=1:obj.Z_max
                im = obj.loadCroppedImage(t,z,ch_short_name,rect);
                score = score_fun(im);
                if score>M
                    M = score;
                    z_best = z;
                end
            end
        end
        
        % exploration
        function explore(obj)
            Exploration.start(obj);
        end
        
        % filter and copy dataset
        function filterAndCopyDataset(obj,ts,zs,channels,rect,output_folder)
            mkdir(output_folder);
            new_t = 1;
            for t=ts
                new_z = 1;
                for z=zs
                    for ch=channels
                        im = obj.loadCroppedImage(t,z,obj.Ch_short_names{ch},rect);
                        imwrite(im,[output_folder '/img_' sprintf('%09i',new_t-1) '_' obj.Ch_names{ch} '_' sprintf('%03i',new_z-1) '.tif']);
                    end
                    new_z = new_z + 1;
                end
                new_t = new_t + 1;
            end
        end
    end
    
end
