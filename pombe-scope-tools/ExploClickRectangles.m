classdef ExploClickRectangles
    
    % same as Exploration except we can click rectangle segments and store them
    % note how we do not rewrite the logic of the basic exploration, we use
    % callback and update functions from the basic exploration methods
    
    methods(Static)
        
        % starting
        function start(fov,output_file)
            % initialize view state and main display
            state.fov = fov;
            state.z = 1; state.t = 1; state.ch = 1; state.rect = [1 1 fov.Width fov.Height];
            im = fov.loadImage(1,state.z,fov.Ch_short_names{state.ch});
            state.h_im = Operations.imshowfit( Operations.normalize(im,'normalize_contrast') ); hold on;
            text_label = strjoin(['ch = ' state.fov.Ch_short_names(state.ch) '/ z = ' num2str(state.z) '/ t = ' num2str(state.t)] );
            state.h_label = text('String',text_label,'Position',[0.6*state.rect(3) 0.05*state.rect(4)],'Color','w','FontSize',20);
            % for output
            state.output_file = output_file;
            % to store the segments
            if exist(state.output_file,'file') %load if exist already
                rectangles = load(state.output_file);
                state.rectangles = rectangles.rectangles;
            else
                state.rectangles = struct('rect',{},'t',{},'z',{},'ch',{});
            end
            state.h_rectangles = [];
            setappdata(gcf,'state',state);
            % update once to show existing rectangles
            ExploClickRectangles.update();
            % bind display to callback
            set(gcf,'KeyPressFcn',@(obj,event)ExploClickRectangles.callback(obj,event.Key));
            waitfor(gcf);
        end
        
        % callback entry
        function callback(~,key)
            % first check for basic change of view (use Exploration code)
            if Exploration.callbackView(key)
                ExploClickRectangles.update();
                return;
            end
            % check for rectangle click
            if ExploClickRectangles.callbackPointClick(key)
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
                [rect,h] = Operations.selectRect(state.rect(3),state.rect(4)); delete(h);
                state.rectangles(end+1).rect = rect;
                state.rectangles(end).rect(1:2) = state.rectangles(end).rect(1:2) + state.rect(1:2) - 1;
                state.rectangles(end).t = state.t;
                state.rectangles(end).z = state.z;
                state.rectangles(end).ch = state.ch;
                setappdata(gcf,'state',state);
                ExploClickRectangles.save2disk();
                ExploClickRectangles.update();
                return;
            end
            % delete last
            if strcmp(key,'backspace')
                state.rectangles(end) = [];
                setappdata(gcf,'state',state);
                ExploClickRectangles.save2disk();
                ExploClickRectangles.update();
                return;
            end
        end
        
        % update
        function update()
            % delete previous rectangles
            state = getappdata(gcf,'state');
            for i=1:length(state.h_rectangles)
                delete(state.h_rectangles(i));
            end
            state.h_rectangles = [];
            % add the rectangles
            for i=1:length(state.rectangles)
                rectangle = state.rectangles(i);
                if (rectangle.t == state.t)
                    if (rectangle.z == state.z)
                        state.h_rectangles(end+1) = ...
                            plot( [rectangle.rect(1),rectangle.rect(1)+rectangle.rect(3),rectangle.rect(1)+rectangle.rect(3),rectangle.rect(1),rectangle.rect(1)] - state.rect(1) + 1 , ...
                                  [rectangle.rect(2),rectangle.rect(2),rectangle.rect(2)+rectangle.rect(4),rectangle.rect(2)+rectangle.rect(4),rectangle.rect(2)] - state.rect(2) + 1 , 'g'); hold on;
                    else
                        state.h_rectangles(end+1) = ...
                            plot( [rectangle.rect(1),rectangle.rect(1)+rectangle.rect(3),rectangle.rect(1)+rectangle.rect(3),rectangle.rect(1),rectangle.rect(1)] - state.rect(1) + 1 , ...
                                  [rectangle.rect(2),rectangle.rect(2),rectangle.rect(2)+rectangle.rect(4),rectangle.rect(2)+rectangle.rect(4),rectangle.rect(2)] - state.rect(2) + 1 , 'r'); hold on;
                    end
                end
            end
            setappdata(gcf,'state',state);
        end
        
        % save to disk
        function save2disk()
            state = getappdata(gcf,'state');
            rectangles = state.rectangles;
            fov = state.fov;
            save(state.output_file,'rectangles','fov');
        end
        
    end
    
end