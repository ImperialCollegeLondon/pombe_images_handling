classdef Operations

    methods(Static)

        % normalize images to [.. 1] or [0 1]
        function im_n = normalize(im,normalize_opt)
            im_n = im;
            if strcmp(normalize_opt,'normalize')
                im_n = double(im);
                im_n = im_n./max(im_n(:));
            elseif strcmp(normalize_opt,'normalize_contrast')
                im_n = double(im_n);
                im_n = (im_n-min(im_n(:)))./(max(im_n(:))-min(im_n(:)));
            end
        end

        % correct a rect to fit in full size
        function rect = makeRectFit(rect,im_w,im_h)
            rect(1) = min([im_w , rect(1)]);
            rect(2) = min([im_h , rect(2)]);
            rect(1) = max([1 , rect(1)]);
            rect(2) = max([1 , rect(2)]);
            rect(3) = min([im_w-rect(1),rect(3)]);
            rect(4) = min([im_h-rect(2),rect(4)]);
            rect(3) = max([1,rect(3)]);
            rect(4) = max([1,rect(4)]);
        end

        %%% add margins to a rect
        % padding: num pixels to add
        % im_w and im_h: dimensions of FULL image
        function p_rect = padRect(rect,padding,im_w,im_h)
            % ensure padding/2 is integer
            if (round(padding/2) ~= padding/2)
                padding = padding + 1;
            end
            p_rect(1) = max( [1 , rect(1) - padding/2] );
            p_rect(2) = max( [1 , rect(2) - padding/2] );
            p_rect(3) = max( [1 , rect(3) + padding] );
            p_rect(4) = max( [1 , rect(4) + padding] );
            if p_rect(1)+p_rect(3) > im_w
                p_rect(3) = im_w - p_rect(1);
            end
            if p_rect(2)+p_rect(4) > im_h
                p_rect(4) = im_h - p_rect(2);
            end
        end

        % user selection of ROI
        function [rect,h] = selectRect(im_w,im_h)
            [x1,y1] = ginput(1); x1 = round(x1); y1 = round(y1);
            h = plot(x1,y1,'g+'); hold on;
            [x2,y2] = ginput(1); x2 = round(x2); y2 = round(y2);
            delete(h);
            h = plot([x1 x1 x2 x2 x1],[y1 y2 y2 y1 y1],'g'); hold on; pause(0.1);
            rect = [ min([x1 x2]) , min([y1 y2]), ...
                max([x1 x2])-min([x1 x2]), max([y1 y2])-min([y1 y2])];
            rect = Operations.makeRectFit(rect,im_w,im_h);
        end

        % user selection of segment
        function [x1,x2,y1,y2,h] = selectSegment()
            [x1,y1] = ginput(1);
            h = plot(x1,y1,'gx'); hold on;
            [x2,y2] = ginput(1);
            delete(h);
            h = plot([x1 x2],[y1 y2],'-gx'); hold on;
            pause(0.1);
        end

        % boundary to rect (i.e. bounding box)
        % careful: here boundary as from bwboundaries, i.e. x is column 2
        function rect = boundary2rect(boundary,im_w,im_h)
            x_min = floor(min(boundary(:,2))) - 1;
            x_max = ceil(max(boundary(:,2))) + 1;
            y_min = floor(min(boundary(:,1))) - 1;
            y_max = ceil(max(boundary(:,1))) + 1;
            rect = [ x_min , y_min  , x_max - x_min , y_max - y_min ];
            rect = Operations.makeRectFit(rect,im_w,im_h);
        end

        % get the most central central component
        function bw_central = getCentralCC(bw)
            regions = regionprops( bw , 'Centroid' , 'PixelIdxList' );
            min_d2 = 1e300; r_min = 0;
            for r=1:length(regions)
                d2 = ( size(bw,2)/2 - regions(r).Centroid(1) )^2 + ( size(bw,1)/2 - regions(r).Centroid(2) )^2 ;
                if d2 < min_d2
                    r_min = r; min_d2 = d2;
                end
            end
            bw_central = false(size(bw));
            bw_central( regions(r_min).PixelIdxList ) = true;
        end

        % get the biggest central components
        function bw_biggest = getBiggestCC( bw )
            regions = regionprops( bw , 'Area' , 'PixelIdxList' );
            areas = [regions.Area];
            [~,r_max] = max(areas);
            bw_biggest = false(size(bw));
            bw_biggest( regions(r_max).PixelIdxList ) = true;
        end

        % imshow but with fitting fullscreen
        function h = imshowfit(I)
            h = imshow(I, 'Border', 'tight', 'InitialMagnification', 'fit');
        end
        
        % convert a contour into geometry stats
        function geom = contour2geom(pts,n_slices,show)
            pts(:,1) = pts(:,1) - min(pts(:,1)) + 1;
            pts(:,2) = pts(:,2) - min(pts(:,2)) + 1;
            if nargin < 3
                show = false;
                if nargin < 2
                    n_slices = 50;
                end
            end
            mask = poly2mask(pts(:,1),pts(:,2),round(max(pts(:)))+10,round(max(pts(:)))+10);
            geom = Operations.mask2geom(mask,n_slices,show);
        end

        % convert a mask into geometry stats like volume (SA inaccurate, to correct)
        function geom = mask2geom(mask,n_slices,show)
            if nargin < 3
                show = false;
                if nargin < 2
                    n_slices = 50;
                end
            end
            props = regionprops(mask,'FilledArea','MajorAxisLength','MinorAxisLength','Centroid','Orientation','BoundingBox');
            if length(props)>1; Operations.imshowfit(mask); error('more than 1 cc ?'); end
            props.BoundingBox(1) = max([1 props.BoundingBox(1)-15]);
            props.BoundingBox(2) = max([1 props.BoundingBox(2)-15]);
            props.BoundingBox(3) = min([size(mask,2) props.BoundingBox(3)+30]);
            props.BoundingBox(4) = min([size(mask,1) props.BoundingBox(4)+30]);
            % find 'tips'
            x1=props.Centroid(1)+0.7*props.MajorAxisLength*cos(-props.Orientation/180*pi);
            x2=props.Centroid(1)-0.7*props.MajorAxisLength*cos(-props.Orientation/180*pi);
            y1=props.Centroid(2)+0.7*props.MajorAxisLength*sin(-props.Orientation/180*pi);
            y2=props.Centroid(2)-0.7*props.MajorAxisLength*sin(-props.Orientation/180*pi);
            [x,y,val] = improfile(mask,[x1 x2],[y1 y2]);
            i1 = find(val>0,1,'first');
            i2 = find(val>0,1,'last');
            x1 = x(i1); x2 = x(i2);
            y1 = y(i1); y2 = y(i2);
            % slicing analysis
            dl = norm([x1 y1] - [x2 y2]) / n_slices;
            V_slices = 0; l_slices = 0; SA_slices = 0;
            last_xm = x1; last_ym = y1;
            if show; Operations.imshowfit(mask); hold on; end
            wsVec=zeros(n_slices,1);
            for s=1:n_slices
                % ref point
                xl = x1 + (dl/2 + (s-1) * dl ) * (x2-x1)/norm([x1 y1] - [x2 y2]);
                yl = y1 + (dl/2 + (s-1) * dl ) * (y2-y1)/norm([x1 y1] - [x2 y2]);
                xl1 = xl + (y1-y2);
                xl2 = xl - (y1-y2);
                yl1 = yl + (x2-x1);
                yl2 = yl - (x2-x1);
                % find border points from the ref point
                [x,y,val] = improfile(mask,[xl1 xl2],[yl1 yl2]);
                i1 = find(val>0,1,'first');
                i2 = find(val>0,1,'last');
                if isempty(i1) || isempty(i2)
                    geom.success = false;
                    return;
                end
                xs1 = x(i1); xs2 = x(i2);
                ys1 = y(i1); ys2 = y(i2);
                if show
                    plot([xs1 xs2],[ys1 ys2],'Color',[1 1 1].*0.5); hold on;
                end
                ws = norm([xs1 ys1]-[xs2 ys2]);
                wsVec(s) = ws;
                V_slices = V_slices + dl * pi * ws^2 / 4;
                xm = (xs1+xs2)/2;
                ym = (ys1+ys2)/2;
                l_slices = l_slices + sqrt( (xm-last_xm)^2 + (ym-last_ym)^2 );
                last_xm = xm; last_ym = ym;
                if (s==1 || s==n_slices)
                    SA_slices = pi * ws^2/4 + dl * pi * ws;
                else
                    SA_slices = SA_slices + dl * pi * ws;
                end
            end
            % finish l_slices
            l_slices = l_slices + sqrt( (x2-last_xm)^2 + (y2-last_ym)^2 );
            % wrap up results in struct
            geom.l = l_slices;
            geom.V = V_slices;
            geom.SA = SA_slices;
            geom.CS = sum(mask(:));
            geom.width_avg = sum(wsVec)/n_slices;
            geom.width_vec = wsVec;
            geom.success = true;
        end

        % convert a mask into geometry stats (dye only!)
        function geom = mask2geomDye(mask,mask2,~)
            % find 'tips'
            props = regionprops(mask,'MajorAxisLength','Centroid','Orientation');
            x1=props.Centroid(1)+0.7*props.MajorAxisLength*cos(-props.Orientation/180*pi);
            x2=props.Centroid(1)-0.7*props.MajorAxisLength*cos(-props.Orientation/180*pi);
            y1=props.Centroid(2)+0.7*props.MajorAxisLength*sin(-props.Orientation/180*pi);
            y2=props.Centroid(2)-0.7*props.MajorAxisLength*sin(-props.Orientation/180*pi);
            [x,y,val] = improfile(mask2,[x1 x2],[y1 y2]);
            i1 = find(val>0,1,'first');
            i2 = find(val>0,1,'last');
            x1 = x(i1); x2 = x(i2);
            y1 = y(i1); y2 = y(i2);
            % slicing analysis
            Lsimple = norm([x1 y1] - [x2 y2]);  n_slices = round(Lsimple/0.3);
            dl = Lsimple/n_slices;
            l_slices = 0;
            last_xm = x1; last_ym = y1;
            xl1 = size(mask2,2);
            for s=1:n_slices
                % ref point
                x1 = x1 + (dl/2 + (s-1) * dl ) * (x2-x1)/norm([x1 y1] - [x2 y2]);
                yl = y1 + (dl/2 + (s-1) * dl ) * (y2-y1)/norm([x1 y1] - [x2 y2]);
                yl1 = yl + (x2-x1);
                yl2 = yl - (x2-x1);
                % find border points from the ref point, dye
                [x,y,val] = improfile(mask2,[xl1 1],[yl1 yl2]);
                i1 = find(val>0,1,'first');
                i2 = find(val>0,1,'last');
                if isempty(i1) || isempty(i2)
                    geom.success = false;
                    return;
                end
                xm = (x(i1)+x(i2))/2;
                ym = (y(i1)+y(i2))/2;
                l_slices = l_slices + sqrt( (xm-last_xm)^2 + (ym-last_ym)^2 );
                last_xm = xm; last_ym = ym;
            end
            l_slices = l_slices + sqrt( (x2-last_xm)^2 + (y2-last_ym)^2 );
            % wrap up results in struct
            geom.l = l_slices;
            geom.success = true;
        end

        % convert a mask into geometry stats (dye only!)
        function geom = mask2geomDye_extended(mask,mask2,n_slices)
            if nargin < 2;  n_slices = 50;  end
            props = regionprops(mask,'FilledArea','MajorAxisLength','MinorAxisLength','Centroid','Orientation');
            % find 'tips'
            x1=props.Centroid(1)+0.7*props.MajorAxisLength*cos(-props.Orientation/180*pi);
            x2=props.Centroid(1)-0.7*props.MajorAxisLength*cos(-props.Orientation/180*pi);
            y1=props.Centroid(2)+0.7*props.MajorAxisLength*sin(-props.Orientation/180*pi);
            y2=props.Centroid(2)-0.7*props.MajorAxisLength*sin(-props.Orientation/180*pi);
            [x,y,val] = improfile(mask2,[x1 x2],[y1 y2]);
            i1 = find(val>0,1,'first');
            i2 = find(val>0,1,'last');
            x1 = x(i1); x2 = x(i2);
            y1 = y(i1); y2 = y(i2);
            % slicing analysis
            dl = norm([x1 y1] - [x2 y2]) / n_slices;
            l_slices = 0;
            last_xm = x1; last_ym = y1;
            wsVec=zeros(n_slices,1); BFwsVec=zeros(n_slices,1);
            xl1 = size(mask2,2); %inside loop: xl1 = xl + (y1-y2);
            xl2 = 1;             %inside loop: xl2 - (y1-y2);
            for s=1:n_slices
                % ref point
                x1 = x1 + (dl/2 + (s-1) * dl ) * (x2-x1)/norm([x1 y1] - [x2 y2]);
                yl = y1 + (dl/2 + (s-1) * dl ) * (y2-y1)/norm([x1 y1] - [x2 y2]);
                yl1 = yl + (x2-x1);
                yl2 = yl - (x2-x1);
                % find border points from the ref point, bf
                [x,y,val] = improfile(mask,[xl1 xl2],[yl1 yl2]);
                i1 = find(val>0,1,'first');
                i2 = find(val>0,1,'last');
                xs1 = x(i1); xs2 = x(i2);
                ys1 = y(i1); ys2 = y(i2);
                BFwsVec(s) = norm([xs1 ys1]-[xs2 ys2]);

                % find border points from the ref point, dye
                [x,y,val] = improfile(mask2,[xl1 xl2],[yl1 yl2]);
                i1 = find(val>0,1,'first');
                i2 = find(val>0,1,'last');
                if isempty(i1) || isempty(i2)
                    geom.success = false;
                    return;
                end
                xs1 = x(i1); xs2 = x(i2);
                ys1 = y(i1); ys2 = y(i2);
                wsVec(s) = norm([xs1 ys1]-[xs2 ys2]);
                xm = (xs1+xs2)/2;
                ym = (ys1+ys2)/2;
                l_slices = l_slices + sqrt( (xm-last_xm)^2 + (ym-last_ym)^2 );
                last_xm = xm; last_ym = ym;
            end
            l_slices = l_slices + sqrt( (x2-last_xm)^2 + (y2-last_ym)^2 );
            % wrap up results in struct
            geom.l = l_slices;
            geom.dye_width_avg = sum(wsVec)/n_slices;
            geom.bf_width_avg  = sum(BFwsVec)/n_slices;
            geom.success = true;
        end
    end

end
