classdef Skeleton
    
    % this class represents geometrically a pombe cell contour using a
    % skeleton, with a medial (potentially curved) line joigning tips and
    % orthogonal segments with potentially different radiuses
    
    properties
        Tips;
        Midline;
        Im;
    end
    
    
    methods
        
        % constructor
        function obj = Skeleton(tips,midline,im)
            obj.Tips = tips;
            obj.Midline = midline;
            obj.Im = im;
        end
        
        % display
        function hs = plotSkeleton(obj,i_active)
            if nargin < 2
                i_active = 0;
            end
            xs = obj.Tips(1).x; ys = obj.Tips(1).y;
            xys_rad = struct('xs',{},'ys',{});
            for i=1:length(obj.Midline)
                xs(end+1) = obj.Midline(i).x; ys(end+1) = obj.Midline(i).y;
            end
            xs(end+1) = obj.Tips(2).x; ys(end+1) = obj.Tips(2).y;
            for i=1:length(obj.Midline)
                ul = [ xs(i+2)-xs(i) , ys(i+2) - ys(i) ];
                ul = ul ./ norm(ul); un = [-ul(2) ul(1)];
                M = [ xs(i+1) , ys(i+1) ];
                N1 = M + un .* obj.Midline(i).width/2;
                N2 = M - un .* obj.Midline(i).width/2;
                xys_rad(end+1).xs = [ N1(1) N2(1) ];
                xys_rad(end).ys = [ N1(2) N2(2) ];
            end
            hs = plot(xs,ys,'-ro'); hold on;
            for i=1:length(xys_rad)
                if i==i_active
                    hs(end+1) = plot(xys_rad(i).xs,xys_rad(i).ys,'-bx'); hold on;
                else
                    hs(end+1) = plot(xys_rad(i).xs,xys_rad(i).ys,'-gx'); hold on;
                end
            end
        end
        
        % adjust manually
        function obj = adjustManual(obj)
            % start
            state.skel = obj;
            state.current = 1;
            state.h_im = Operations.imshowfit(state.skel.Im); hold on;
            state.hs_plots = obj.plotSkeleton(state.current);
            setappdata(gcf,'state',state);
            set(gcf,'KeyPressFcn',@(obj,event)Skeleton.callbackAdjustManual(obj,event.Key));
            waitfor(gcf);
            state = getappdata(gcf,'state');
            close;
            obj = state.skel;
        end
        
        % compute the ul and un vectors for a midline point
        function [ul,un,l] = getMidlinePoint(obj,i)
            if i==1
                A = [ obj.Tips(1).x , obj.Tips(1).y ];
            else
                A = [ obj.Midline(i-1).x , obj.Midline(i-1).y ];
            end
            if i==length(obj.Midline)
                B = [ obj.Tips(2).x , obj.Tips(2).y ];
            else
                B = [ obj.Midline(i+1).x , obj.Midline(i+1).y ];
            end
            ul = (B-A) ./ norm(B-A);
            un = [-ul(2) ul(1)];
            l = norm(B-A);
        end
        
        % add midline points
        function obj = addMidlinePoints(obj)
            new_midline = struct('x',{},'y',{},'width',{});
            for i=1:length(obj.Midline)
                [ul,~,l] = obj.getMidlinePoint(i);
                new_midline(end+1).width = obj.Midline(i).width;
                new_midline(end).x = obj.Midline(i).x - ul(1) * l / 4;
                new_midline(end).y = obj.Midline(i).y - ul(2) * l / 4;
                new_midline(end+1).width = obj.Midline(i).width;
                new_midline(end).x = obj.Midline(i).x;
                new_midline(end).y = obj.Midline(i).y;
            end
            new_midline(end+1).width = obj.Midline(i).width;
            new_midline(end).x = obj.Midline(i).x + ul(1) * l / 4;
            new_midline(end).y = obj.Midline(i).y + ul(2) * l / 4;            
            obj.Midline = new_midline;
        end
        
        % delete midline points
        function obj = deleteMidlinePoints(obj)
            obj.Midline = obj.Midline(2:2:end);
        end
        
        % compute the length
        function l = computeLength(obj)
           l = 0;
            for i=1:length(obj.Midline)
                [~,~,delta_l] = obj.getMidlinePoint(i);
                l = l + delta_l / 2;
            end
            l = l + delta_l / 2;
        end
        % compute the average width
        function w = computeAvgWidth(obj)
            w = mean([obj.Midline.width]); 
        end
        % compute median width
        function w_med = computeMedianWidth(obj)
            w_med = median([obj.Midline.width]);
        end
        % compute simple volume
        function V = computeSimpleVolume(obj)
            V = pi * obj.computeAvgWidth()^2 / 4 * obj.computeLength();
        end
        % compute volume
        function V = computeVolume(obj)
            mask = obj.computeContour;
            geom = Operations.mask2geom(mask);
            V = geom.V;
        end
        % compute SA
        function SA = computeSA(obj)
            mask = obj.computeContour;
            geom = Operations.mask2geom(mask);
            SA = geom.SA;
        end        
        
        % get the polygon
        function [xs,ys] = computePolygon(obj)
            xs = []; ys = [];
            % first side
            for i=1:length(obj.Midline)
                [~,un,~] = obj.getMidlinePoint(i);
                xs(end+1) = obj.Midline(i).x  + un(1) .* obj.Midline(i).width/2;
                ys(end+1) = obj.Midline(i).y  + un(2) .* obj.Midline(i).width/2;
            end
            % first cap
            n_cap_pts = 8;
            ul = [ obj.Midline(end).x - obj.Midline(end-1).x , obj.Midline(end).y - obj.Midline(end-1).y ];
            ul = ul ./ norm(ul);
            d = sqrt((obj.Tips(2).x-obj.Midline(end).x)^2+(obj.Tips(2).y-obj.Midline(end).y)^2);
            if d > obj.Midline(end).width/2
                delta = d - obj.Midline(end).width/2;
            else
                delta = 0;
            end
            for i=1:n_cap_pts
                alpha = (n_cap_pts-i+1)*pi/(n_cap_pts+1) - pi/2;
                xs(end+1) = obj.Midline(end).x + delta*ul(1) + obj.Midline(end).width/2 * ( cos(alpha) * ul(1) - sin(alpha) * ul(2) );
                ys(end+1) = obj.Midline(end).y + delta*ul(2) + obj.Midline(end).width/2 * ( sin(alpha) * ul(1) + cos(alpha) * ul(2) );
            end
            % second side
            for i=1:length(obj.Midline)
                [~,un,~] = obj.getMidlinePoint(length(obj.Midline)-i+1);
                xs(end+1) = obj.Midline(length(obj.Midline)-i+1).x  - un(1) .* obj.Midline(length(obj.Midline)-i+1).width/2;
                ys(end+1) = obj.Midline(length(obj.Midline)-i+1).y  - un(2) .* obj.Midline(length(obj.Midline)-i+1).width/2;
            end          
            % second cap
            n_cap_pts = 8;
            ul = [ obj.Midline(2).x - obj.Midline(1).x , obj.Midline(2).y - obj.Midline(1).y ];
            ul = ul ./ norm(ul);
            d = sqrt((obj.Tips(1).x-obj.Midline(1).x)^2+(obj.Tips(1).y-obj.Midline(1).y)^2);
            if d > obj.Midline(1).width/2
                delta = d - obj.Midline(1).width/2;
            else
                delta = 0;
            end
            for i=1:n_cap_pts
                alpha = (n_cap_pts-i+1)*pi/(n_cap_pts+1) - pi/2;
                xs(end+1) = obj.Midline(1).x - delta*ul(1) - obj.Midline(1).width/2 * ( cos(alpha) * ul(1) - sin(alpha) * ul(2) );
                ys(end+1) = obj.Midline(1).y - delta*ul(2) - obj.Midline(1).width/2 * ( sin(alpha) * ul(1) + cos(alpha) * ul(2) );
            end
            % close
            xs(end+1) = xs(1); ys(end+1) = ys(1);
        end
                
        % compute the contour
        function mask = computeContour(obj)
            [xs,ys] = obj.computePolygon();
            mask = poly2mask(xs,ys,size(obj.Im,1),size(obj.Im,2));
        end
        
        % compute the cross section
        function S = computeCrossSectionSurface(obj)
            mask = obj.computeContour();
            S = sum(mask(:));
        end
        
        % 
        function obj = improveWidths(obj,delta)
            new_midline = obj.Midline;
            for i=1:length(obj.Midline)
                [~,un,~] = obj.getMidlinePoint(i);
                % point 1
                N1 = [obj.Midline(i).x + obj.Midline(i).width/2 * un(1) , obj.Midline(i).y + obj.Midline(i).width/2 * un(2)] ;
                P1 = N1 + delta .* un; P2 = N1 - delta .* un;
                [x,y,v] = improfile(obj.Im,[P1(1) P2(1)],[P1(2) P2(2)]);
                [~,j_min] = min(v); N1_new = [x(j_min),y(j_min)];
%                 plot([P1(1) P2(1)],[P1(2) P2(2)],'m'); hold on;
                % point 2
                N2 = [obj.Midline(i).x - obj.Midline(i).width/2 * un(1) , obj.Midline(i).y - obj.Midline(i).width/2 * un(2)] ;
                P1 = N2 + delta .* un; P2 = N2 - delta .* un;
                [x,y,v] = improfile(obj.Im,[P1(1) P2(1)],[P1(2) P2(2)]);
                [~,j_min] = min(v); N2_new = [x(j_min),y(j_min)];
%                 plot([P1(1) P2(1)],[P1(2) P2(2)],'m'); hold on;
                % update midline point if not too big change
                old_width = obj.Midline(i).width;
                new_width = norm(N1_new-N2_new);
                if abs(new_width-old_width) < 0.15 * old_width
                    new_midline(i).x = (N1_new(1) + N2_new(1))/2;
                    new_midline(i).y = (N1_new(2) + N2_new(2))/2;
                    new_midline(i).width = norm(N1_new-N2_new);
                end
            end
            obj.Midline = new_midline;
        end
        
        % iterative improvements
        function obj = iterativeImprovement(obj,delta,n_iter)
            for i=1:n_iter
                obj = obj.improveWidths(delta);
                obj = obj.addMidlinePoints();
            end
            obj = obj.improveWidths(delta);
        end
        
        
    end
    
    
    
    
    
    
    
    methods(Static)
        
        % build a straight skeleton from tips
        function skel = buildStraightSkeleton(tips,n_midline,width)
            midline = struct('x',{},'y',{},'width',{});
            for i=1:n_midline
                midline(end+1).x = tips(1).x + (tips(2).x - tips(1).x) * i / (n_midline+1);
                midline(end).y = tips(1).y + (tips(2).y - tips(1).y) * i / (n_midline+1);
                midline(end).width = width;
            end
            skel = Skeleton(tips,midline,[]);
        end
        
        
        % create skeleton straigth from clicked tips
        function skel = clickTipsBuildStraightSkeleton(im,n_midline,width)
            Operations.imshowfit(im); hold on;
            [x1,x2,y1,y2,~] = Operations.selectSegment(); close;
            tips=struct('x',{},'y',{});
            tips(1).x = x1; tips(1).y = y1; tips(2).x = x2; tips(2).y = y2;
            skel = Skeleton.buildStraightSkeleton(tips,n_midline,width);
            skel.Im = im;
        end
        
        
        function testing_2()
            tips=struct('x',{},'y',{});
            tips(1).x = -2; tips(1).y = 0;
            tips(2).x = 2; tips(2).y = 0;
            skel = Skeleton.buildStraightSkeleton(tips,5,0.5);
            skel.Midline(3).width = 0.6;
            skel.plotSkeleton(1);
        end
        
        
        function skel = testing_1()
            %
            im = imread('long-cell.tif'); im = Operations.normalize(im,'normalize_contrast');
            %
            skel = Skeleton.clickTipsBuildStraightSkeleton(im,1,15);
            %
            skel = skel.adjustManual();
        end
        
        
        function callbackAdjustManual(~,key)
            state = getappdata(gcf,'state');
            state.delta = 0.5;
            i = state.current; state.skel;
            if strcmp(key,'w') % increase width
                state.skel.Midline(i).width = state.skel.Midline(i).width + state.delta; setappdata(gcf,'state',state); Skeleton.updateAdjustManual(); return;
            end
            if strcmp(key,'s') % decrease width
                state.skel.Midline(i).width = state.skel.Midline(i).width - state.delta; setappdata(gcf,'state',state); Skeleton.updateAdjustManual(); return;
            end
            if strcmp(key,'a') || strcmp(key,'d') % move center 'left' / 'right'
                [~,un] = state.skel.getMidlinePoint(state.current);
                mult = 1; if strcmp(key,'d'); mult = -1; end
                state.skel.Midline(i).x = state.skel.Midline(i).x + mult * un(1) * state.delta;
                state.skel.Midline(i).y = state.skel.Midline(i).y + mult * un(2) * state.delta;
                setappdata(gcf,'state',state); Skeleton.updateAdjustManual(); return;
            end
            if strcmp(key,'equal')
                state.current = state.current + 1; 
                if state.current == length(state.skel.Midline) + 1
                    state.current = 1;
                end
                setappdata(gcf,'state',state); Skeleton.updateAdjustManual(); return;
            end
            if strcmp(key,'hyphen')
                state.current = state.current - 1; 
                if state.current == 0
                    state.current = length(state.skel.Midline);
                end
                setappdata(gcf,'state',state); Skeleton.updateAdjustManual(); return;
            end
            if strcmp(key,'rightbracket')
                state.skel = state.skel.addMidlinePoints();
                setappdata(gcf,'state',state); Skeleton.updateAdjustManual(); return;
            end
            if strcmp(key,'leftbracket')
                state.skel = state.skel.deleteMidlinePoints();
                setappdata(gcf,'state',state); Skeleton.updateAdjustManual(); return;
            end            
            if strcmp(key,'1')
                state.delta = max([0.5,state.delta - 0.5]);
                setappdata(gcf,'state',state); Skeleton.updateAdjustManual(); return;
            end
            if strcmp(key,'2')
                state.delta = state.delta + 0.5;
                setappdata(gcf,'state',state); Skeleton.updateAdjustManual(); return;
            end
            if strcmp(key,'space')
                state = getappdata(gcf,'state');
                close; figure('visible','off');
                setappdata(gcf,'state',state);
                return;
            end
            
        end
        
        function updateAdjustManual()
            state = getappdata(gcf,'state');
            if isfield(state,'hs_plots')
                delete(state.hs_plots);
            end
            state.hs_plots = state.skel.plotSkeleton(state.current);
            setappdata(gcf,'state',state);
        end
        
    end
    
end












%         function userBisectionCallback(~,key)
%             state = getappdata(gcf,'state');
%             if strcmp(key,'w') % increase width
%                 state.width = state.width + 1; setappdata(gcf,'state',state); Skeleton.update(); return;
%             end
%             if strcmp(key,'s') % decrease width
%                 state.width = state.width - 1; setappdata(gcf,'state',state); Skeleton.update(); return;
%             end
%             if strcmp(key,'a') % move center 'left'
%                 state.delta_center = state.delta_center - sign(state.un(1)); setappdata(gcf,'state',state); Skeleton.update(); return;
%             end
%             if strcmp(key,'d') % move center 'right'
%                 state.delta_center = state.delta_center + sign(state.un(1)); setappdata(gcf,'state',state); Skeleton.update(); return;
%             end
%         end
%
%         function update()
%             state = getappdata(gcf,'state');
%             state.ul = (state.B-state.A) ./ norm(state.B-state.A);
%             state.un = [-state.ul(2) state.ul(1)];
%             M = state.A + state.ul .* norm(state.B-state.A) ./ 2 + state.un .* state.delta_center;
%             N1 = M + state.un .* state.width/2;
%             N2 = M - state.un .* state.width/2;
%             if isfield(state,'h_plot')
%                 delete(state.h_plot);
%             end
%             state.h_plot = plot( [N1(1) N2(1)] , [N1(2) N2(2)] , '--gx'); hold on;
%             setappdata(gcf,'state',state);
%         end