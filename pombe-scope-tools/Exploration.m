classdef Exploration
    
    % those methods allow to do a basic exploration of a field of view:
    % move, zoom, change channel...
    % other methods build on those to also allow to click points, click
    % ROIs, segment cells...
    % importantly, no need to re-write code for the basic exploration, as
    % those functions can be-reused in the other 'richer' exploration
    % scenario (cf ExploClickPoints)
    
    methods(Static)
        
        % starting an exploration
        function start(fov)
            state.fov = fov;
            state.z = 1; state.t = 1; state.ch = 1; state.rect = [1 1 fov.Width fov.Height];
            im = fov.loadImage(1,state.z,fov.Ch_short_names{state.ch});
            state.h_im = Operations.imshowfit( Operations.normalize(im,'normalize_contrast') ); hold on;
            text_label = strjoin(['ch = ' state.fov.Ch_short_names(state.ch) '/ z = ' num2str(state.z) '/ t = ' num2str(state.t)] );
            state.h_label = text('String',text_label,'Position',[0.6*state.rect(3) 0.05*state.rect(4)],'Color','w','FontSize',20);
            state.threshold = 0.5; state.threshold_dir = 1;
            setappdata(gcf,'state',state);
            set(gcf,'KeyPressFcn',@(obj,event)Exploration.callback(obj,event.Key));
            waitfor(gcf);
        end
        
        % callback entry
        function callback(~,key)
            if Exploration.callbackView(key)
                return;
            end
            if strcmp(key,'space')
                close; return;
            end      
        end
        
        % callback check: view changes
        function valid_key = callbackView(key)
            valid_key = false;
            state = getappdata(gcf,'state');
            % change of z-slice
            if strcmp(key,'equal')
                valid_key = true;
                state.z = min([state.z + 1,state.fov.Z_max]);
                setappdata(gcf,'state',state); Exploration.updateView(); return;
            end
            if strcmp(key,'hyphen')
                valid_key = true;
                state.z = max([state.z - 1,-3]);
                setappdata(gcf,'state',state); Exploration.updateView(); return;
            end
            % change channel with c
            if strcmp(key,'c')
                valid_key = true;
                state.ch = state.ch + 1;
                if (state.ch > state.fov.N_ch); state.ch = 1; end
                setappdata(gcf,'state',state); Exploration.updateView(); return;
            end
            % change time with [ and ]
            if strcmp(key,'rightbracket')
                valid_key = true;
                state.t = min([state.t+1 , state.fov.T_max]);
                setappdata(gcf,'state',state); Exploration.updateView(); return;
            end
            if strcmp(key,'leftbracket')
                valid_key = true;
                state.t = max([state.t-1 , 1]);
                setappdata(gcf,'state',state); Exploration.updateView(); return;
            end            
            % change of zoom
            if strcmp(key,'p') || strcmp(key,'o')
                valid_key = true;
                delta = state.rect(3);
                xc = state.rect(1)+delta/2; yc = state.rect(2)+delta/2;
                if strcmp(key,'p')
                    delta = delta * 0.5;
                else
                    delta = delta * 2;
                end
                state.rect(1) = xc-delta/2; state.rect(2) = yc-delta/2;
                state.rect(3) = delta; state.rect(4) = delta;
                state.rect(1) = min([state.fov.Width , state.rect(1)]);
                state.rect(2) = min([state.fov.Height , state.rect(2)]);
                state.rect(1) = max([1 , state.rect(1)]);
                state.rect(2) = max([1 , state.rect(2)]);
                state.rect(3) = min([state.fov.Width-state.rect(1),state.rect(3)]);
                state.rect(4) = min([state.fov.Height-state.rect(2),state.rect(4)]);
                state.rect(3) = max([10,state.rect(3)]);
                state.rect(4) = max([10,state.rect(4)]);
                setappdata(gcf,'state',state); Exploration.updateView(); return;
            end    
            % zoom via roi selection
            if strcmp(key,'i')
                valid_key = true;
                [rect,h] = Operations.selectRect(state.rect(3),state.rect(4));
                state.rect(1:2) = rect(1:2) + state.rect(1:2) - 1;
                state.rect(3:4) = rect(3:4);
                delete(h);
                setappdata(gcf,'state',state); Exploration.updateView(); return;
            end
            % displacement
            if strcmp(key,'a') || strcmp(key,'d') || strcmp(key,'w') || strcmp(key,'s')
                valid_key = true;
                delta = state.rect(3);
                if strcmp(key,'a'); state.rect(1) = max([1 , state.rect(1)-0.1*delta]); end
                if strcmp(key,'d'); state.rect(1) = min([state.fov.Width , state.rect(1)+0.1*delta]); end
                if strcmp(key,'w'); state.rect(2) = max([1 , state.rect(2)-0.1*delta]); end
                if strcmp(key,'s'); state.rect(2) = min([state.fov.Height , state.rect(2)+0.1*delta]); end
                setappdata(gcf,'state',state); Exploration.updateView(); return;
            end
            % change of threshold
            if strcmp(key,'t')
                state.threshold = max([0,state.threshold-0.015]);
                setappdata(gcf,'state',state); Exploration.updateView(true); return;
            end
            if strcmp(key,'y')
                state.threshold = min([1,state.threshold+0.015]);
                setappdata(gcf,'state',state); Exploration.updateView(true); return;
            end
            if strcmp(key,'u')
                state.threshold_dir = 1 - state.threshold_dir;
                setappdata(gcf,'state',state); Exploration.updateView(true); return;
            end            
        end
        
        % update view
        function updateView(show_threshold)
            % by default no show of threshold
            if nargin<1
                show_threshold = false;
            end 
            % get state
            state = getappdata(gcf,'state');
            % load image
            im = Operations.normalize(imcrop(state.fov.loadImage(state.t,state.z,state.fov.Ch_short_names{state.ch}),state.rect),'normalize_contrast');
            % save image
            state.im = im;
            % apply threshold if show threshold
            if (show_threshold && state.threshold<1)
                if state.threshold_dir==1
                    im(im > state.threshold) = 1;
                else
                    im(im < state.threshold) = 1;
                end
            end
            % set image
            set(state.h_im,'CData',im);
            set(gca,'XLim',[1 size(im,2)] + [-0.5 0.5]);
            set(gca,'YLim',[1 size(im,1)] + [-0.5 0.5]);
            % draw labels
            delete(state.h_label);
            z_label = num2str(state.z);
            special_labels = { 'CV' , 'min' , 'max' , 'min x max' };
            if (state.z<1)
                z_label = special_labels{-state.z+1};
            end
            text_label = strjoin(['ch = ' state.fov.Ch_short_names(state.ch) '/ z = ' z_label '/ t = ' num2str(state.t)] );
            state.h_label = text('String',text_label,'Position',[0.6*state.rect(3) 0.05*state.rect(4)],'Color','w','FontSize',20);
            % save state
            setappdata(gcf,'state',state);
        end
        
    end
    
end