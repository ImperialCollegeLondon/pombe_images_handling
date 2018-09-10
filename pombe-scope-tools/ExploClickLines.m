classdef ExploClickLines
    
    % same as Exploration except we can click line segments and store them
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
                lines = load(state.output_file);
                state.lines = lines.lines;
            else
                state.lines = struct('x1',{},'y1',{},'x2',{},'y2',{},'t',{},'z',{},'ch',{});
            end
            state.h_lines = [];
            setappdata(gcf,'state',state);
            % update once to show existing lines
            ExploClickLines.update();
            % bind display to callback
            set(gcf,'KeyPressFcn',@(obj,event)ExploClickLines.callback(obj,event.Key));
            waitfor(gcf);
        end
        
        % callback entry
        function callback(~,key)
            % first check for basic change of view (use Exploration code)
            if Exploration.callbackView(key)
                ExploClickLines.update();
                return;
            end
            % check for line click
            if ExploClickLines.callbackPointClick(key)
                return;
            end
            % quitting
            if strcmp(key,'space')
                close; return;
            end
        end
        
        % callback check: line click
        function valid_key = callbackPointClick(key)
            valid_key = false;
            state = getappdata(gcf,'state');
            % click segment
            if strcmp(key,'x')
                valid_key = true;
                [x1,x2,y1,y2,h] = Operations.selectSegment(); delete(h);
                state.lines(end+1).x1 = x1 + state.rect(1) - 1;
                state.lines(end).y1 = y1 + state.rect(2) - 1;
                state.lines(end).x2 = x2 + state.rect(1) - 1;
                state.lines(end).y2 = y2 + state.rect(2) - 1;
                state.lines(end).t = state.t;
                state.lines(end).z = state.z;
                state.lines(end).ch = state.ch;
                setappdata(gcf,'state',state);
                ExploClickLines.save2disk();
                ExploClickLines.update();
                return;
            end
            % delete last one
            if strcmp(key,'backspace')
                state.lines(end) = [];
                setappdata(gcf,'state',state);
                ExploClickLines.save2disk();
                ExploClickLines.update();
                return;                
            end
        end
        
        % update
        function update()
            % delete previous lines
            state = getappdata(gcf,'state');
            for i=1:length(state.h_lines)
                delete(state.h_lines(i));
            end
            state.h_lines = [];
            % add the lines
            for i=1:length(state.lines)
                line = state.lines(i);
                if (line.t == state.t)
                    if (line.z == state.z)
                        state.h_lines(end+1) = plot( [line.x1 line.x2] - state.rect(1) + 1 , [line.y1 line.y2] - state.rect(2) + 1 , '-g+'); hold on;
                    else
                        state.h_lines(end+1) = plot( [line.x1 line.x2] - state.rect(1) + 1 , [line.y1 line.y2] - state.rect(2) + 1 , '-r+'); hold on;
                    end
                end
            end
            setappdata(gcf,'state',state);
        end
        
        % save to disk
        function save2disk()
            state = getappdata(gcf,'state');
            lines = state.lines;
            fov = state.fov;
            save(state.output_file,'lines','fov');
        end
        
    end
    
end