classdef Grid
    
    properties
        Base_path
        Grid_size_x
        Grid_size_y
        Ch_names
        Ch_short_names
        FOVs
    end
    
    methods
        % constructor
        function obj = Grid(base_path,grid_size_x,grid_size_y,ch_names,ch_short_names)
            obj.Base_path = base_path;
            obj.Grid_size_x = grid_size_x;
            obj.Grid_size_y = grid_size_y;
            obj.Ch_names = ch_names;
            obj.Ch_short_names = ch_short_names;
            obj.FOVs = cell(grid_size_x,grid_size_y);
            for i=1:obj.Grid_size_x
                for j=1:obj.Grid_size_y
                    obj.FOVs{i,j} = FieldOfView(obj.Ch_names,obj.Ch_short_names,[base_path '-Pos_' sprintf('%03i',i-1) '_' sprintf('%03i',j-1)]);
                end
            end
        end
        % explore on of the FOV
        function explore(obj,grid_x,grid_y)
            obj.FOVs{grid_x,grid_y}.explore();
        end
    end
    
end

