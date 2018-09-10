classdef ExploClickPoints
    
    % same as Exploration except we can click points and store them
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
            % to store the points
            if exist(state.output_file,'file') %load if exist already
                points = load(state.output_file);
                state.points = points.points;
            else
                state.points = struct('x',{},'y',{},'t',{},'z',{},'ch',{});
            end
            state.h_points = [];
            setappdata(gcf,'state',state);
            % update once to show existing points
            ExploClickPoints.update();
            % bind display to callback
            set(gcf,'KeyPressFcn',@(obj,event)ExploClickPoints.callback(obj,event.Key));
            waitfor(gcf);
        end
        
        % callback entry
        function callback(~,key)
            % first check for basic change of view (use Exploration code)
            if Exploration.callbackView(key)
                ExploClickPoints.update();
                return;
            end
            % check for point click
            if ExploClickPoints.callbackPointClick(key)
                return;
            end
            % quitting
            if strcmp(key,'space')
                close; return;
            end
        end
        
        % callback check: point click
        function valid_key = callbackPointClick(key)
            valid_key = false;
            state = getappdata(gcf,'state');
            % new point
            if strcmp(key,'x')
                valid_key = true;
                pt = ginput(1);
                state.points(end+1).x = pt(1) + state.rect(1) - 1;
                state.points(end).y = pt(2)  + state.rect(2) - 1;
                state.points(end).t = state.t;
                state.points(end).z = state.z;
                state.points(end).ch = state.ch;
                setappdata(gcf,'state',state);
                ExploClickPoints.save2disk();
                ExploClickPoints.update();
                return;
            end
            % delete last clicked point
            if strcmp(key,'backspace')
                state.points(end) = [];
                setappdata(gcf,'state',state);
                ExploClickPoints.save2disk();
                ExploClickPoints.update();                
            end
        end
        
        % update
        function update()
            % delete previous points
            state = getappdata(gcf,'state');
            for i=1:length(state.h_points)
                delete(state.h_points(i));
            end
            state.h_points = [];
            % add the points
            for i=1:length(state.points)
                point = state.points(i);
                if (point.t == state.t)
                    if (point.z == state.z)
                        state.h_points(end+1) = plot( point.x - state.rect(1) + 1 , point.y - state.rect(2) + 1 , 'go'); hold on;
                    else
                        state.h_points(end+1) = plot( point.x - state.rect(1) + 1 , point.y - state.rect(2) + 1 , 'ro'); hold on;
                    end
                end
            end
            setappdata(gcf,'state',state);
        end
        
        % save to disk
        function save2disk()
            state = getappdata(gcf,'state');
            points = state.points;
            fov = state.fov;
            save(state.output_file,'points','fov');
        end
        
    end
    
end