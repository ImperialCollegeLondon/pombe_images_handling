classdef PtsContour
    
    % this class represents a contour made of joined points
    % can be interactively modified to correct a segmentation
    
    properties
        Points; % one point per row, colums x , y
        Im;
    end
    
    
    methods
        
        % constructor
        function obj = PtsContour(points,im)
            [points(:,1),points(:,2)] = poly2cw(points(:,1),points(:,2)); % ensure clockwise
            obj.Points = points;
            obj.Im = im;
        end
        
        % display
        function hs = plot(obj,i_active)
            if nargin < 2
                i_active = 0; i_next = 0;
            else
                i_next = i_active + 1;
                if i_next > length(obj.Points)
                    i_next = 1;
                end
            end
            pts = [ obj.Points ; obj.Points(1,:) ];
            hs = plot(pts(:,1),pts(:,2),'r'); hold on;
            for i=1:length(obj.Points)
                if i==i_active
                    hs(end+1) = plot(obj.Points(i,1),obj.Points(i,2),'g+'); hold on;
                elseif i==i_next
                    hs(end+1) = plot(obj.Points(i,1),obj.Points(i,2),'b+'); hold on;
                else
                    hs(end+1) = plot(obj.Points(i,1),obj.Points(i,2),'r+'); hold on;
                end
            end
        end
        
        % adjust manually
        function obj = adjustManual(obj)
            % start
            state.contour = obj;
            state.current = 1;
            state.h_im = Operations.imshowfit(state.contour.Im); hold on;
            state.hs_plots = obj.plot(state.current);
            setappdata(gcf,'state',state);
            set(gcf,'KeyPressFcn',@(obj,event)PtsContour.callbackAdjustManual(obj,event.Key));
            waitfor(gcf);
            state = getappdata(gcf,'state');
            close;
            obj = state.contour;
        end

        % compute the mask of the contour
        function [mask,im] = toMask(obj,scale)
            if nargin < 2
                scale = 1;
            end
            im = imresize(obj.Im,scale);
            pts = obj.Points .* scale;
            mask = poly2mask(pts(:,1),pts(:,2),size(im,1),size(im,2));
        end
        
        % compute geometry
        function geom = toGeometry(obj,scale,show)
            if nargin < 3
                show = false;
            end
            [mask,im] = obj.toMask(scale);
            Operations.imshowfit(im); hold on;
            geom = Operations.mask2geom(mask,50,show);
            geom.V = geom.V / scale^3;
            geom.CS = geom.CS / scale^2;
        end
        
        
        % add point btw current and next
        function obj = addPoint(obj,i_current)
            i_next = i_current + 1;
            if i_next > length(obj.Points); i_next = 1; end
            new_pt = (obj.Points(i_current,:) + obj.Points(i_next,:)) ./ 2;
            obj.Points = [obj.Points(1:i_current,:) ; new_pt ; obj.Points(i_current+1:end,:) ];
        end
        
        % delete point
        function obj = deletePoint(obj,i_current)
            obj.Points(i_current,:) = [];
        end
        
    end
    
    
    
    
    
    
    
    methods(Static)

        % testing
        function testing_1()
            im = rand(80,60) ./ 5;
            points = [ 10,10 ; 10,15 ; 10,20 ; 10,25 ; 20,25 ; 20,20 ; 20,15 ; 20,10 ];
            contour = PtsContour(points,im);
            contour = contour.adjustManual();
            contour.adjustManual();
        end
        
        % build a contour looking like a rod from tips and radius
        function contour = buildRod(tips,n_body,n_cap,width)
            error('not implemented yet');
        end
        
        % for manual adjustment of contour (callback)
        function callbackAdjustManual(~,key)
            state = getappdata(gcf,'state');
            state.delta = 0.5;
            i = state.current;
            if strcmp(key,'w') % move up
                state.contour.Points(i,2) = state.contour.Points(i,2) - state.delta;
                setappdata(gcf,'state',state); PtsContour.updateAdjustManual(); return;
            end
            if strcmp(key,'s') % move down
                state.contour.Points(i,2) = state.contour.Points(i,2) + state.delta;
                setappdata(gcf,'state',state); PtsContour.updateAdjustManual(); return;
            end
            if strcmp(key,'a') % move left
                state.contour.Points(i,1) = state.contour.Points(i,1) - state.delta;
                setappdata(gcf,'state',state); PtsContour.updateAdjustManual(); return;
            end
            if strcmp(key,'d') % move right
                state.contour.Points(i,1) = state.contour.Points(i,1) + state.delta;
                setappdata(gcf,'state',state); PtsContour.updateAdjustManual(); return;
            end                        
            if strcmp(key,'hyphen')
                state.current = state.current + 1; 
                if state.current == length(state.contour.Points) + 1
                    state.current = 1;
                end
                setappdata(gcf,'state',state); PtsContour.updateAdjustManual(); return;
            end
            if strcmp(key,'equal')
                state.current = state.current - 1; 
                if state.current == 0
                    state.current = length(state.contour.Points);
                end
                setappdata(gcf,'state',state); PtsContour.updateAdjustManual(); return;
            end
            if strcmp(key,'rightbracket')
                state.contour = state.contour.addPoint(state.current); state.current = state.current + 1;
                setappdata(gcf,'state',state); PtsContour.updateAdjustManual(); return;
            end
            if strcmp(key,'backspace')
                state.contour = state.contour.deletePoint(state.current); 
                state.current = state.current - 1; 
                if state.current == 0; state.current = length(state.contour.Points); end
                setappdata(gcf,'state',state); PtsContour.updateAdjustManual(); return;
            end                     
            if strcmp(key,'1')
                state.delta = max([0.5,state.delta - 0.5]);
                setappdata(gcf,'state',state); PtsContour.updateAdjustManual(); return;
            end
            if strcmp(key,'2')
                state.delta = state.delta + 0.5;
                setappdata(gcf,'state',state); PtsContour.updateAdjustManual(); return;
            end
            if strcmp(key,'space')
                state = getappdata(gcf,'state');
                close; figure('visible','off');
                setappdata(gcf,'state',state);
                return;
            end
            
        end
        
        % for manual adjustment of contour (update)
        function updateAdjustManual()
            state = getappdata(gcf,'state');
            if isfield(state,'hs_plots')
                delete(state.hs_plots);
            end
            state.hs_plots = state.contour.plot(state.current);
            setappdata(gcf,'state',state);
        end
        
    end
    
end




