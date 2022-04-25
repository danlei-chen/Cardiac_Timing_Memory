classdef DynamicTable < handle
    %DynamicTable To receive and send input to a GUI table
    % Ian Kleckner
    % Interdisciplinary Affective Science Lab (IASL)
    % Continuous Flash Suppression (CFS)
    %
    % 2011/04/22 Start coding
    % 2011/09/22 Add getDataStructure
    % 2011/09/27 Allow STRING type
    % 2011/10/04 Write to log file
    % 2012/02/22 More informative error message for out of bounds error
    % 2012/04/20 More informative error messages
    %
    % TODO
    %
    
    properties (SetAccess = private)
        % These can only be set via functions (but read directly via ".")
        modes             = '';   % Mode for the row (NUMERIC, LIST, BOOLEAN)
        names            = {};
        values           = {};
        descriptions     = {};
        Nrows            = 0;
        
        isVisible       = {};
        
        values_Min      = {};
        values_Max     = {};
        
        table_handle	= NaN;  % Handle to table for input/display
    end
    
    methods
        function obj = DynamicTable( table_handle )
            %% Default constructor
            obj.table_handle = table_handle;
            
            columnNames = {'Name', 'Value', 'Description'};
            set(table_handle, 'ColumnName', columnNames);
            set(table_handle, 'RowName', {});
            set(table_handle, 'ColumnEditable', [false,true,false]);
            set(table_handle, 'ColumnWidth', {150, 60, 450});
        end
        
        function addRow(obj, mode, name, value, description, value_Min, value_Max)
            %% Add the row to the table
            
            % Check all inputs
           	switch upper(mode)
                case {'NUMERIC', 'HEADING', 'BOOLEAN', 'FOLDER', 'FILE', 'STRING'}
                    % Input is OK, commit settings
                    r = obj.Nrows+1;
                    obj.modes{r}            = mode;
                    obj.names{r}            = name;
                    obj.values{r}           = value;
                    obj.descriptions{r}     = description;
                    obj.values_Min{r}       = value_Min;
                    obj.values_Max{r}       = value_Max;
                    obj.isVisible{r}        = true;
                    obj.Nrows               = r;            
                    
                otherwise
                    error('Bad input');            
            end
            
            
        end
        
        function displayTable( obj )
            %% Prepare table_data cell array and display to table
            table_data = cell(1, 3);
            
            rDisplay = 0;
            for r = 1:obj.Nrows                
                if( obj.isVisible{r} )
                    rDisplay = rDisplay+1;
                    switch obj.modes{r}
                        case 'HEADING'
                            table_data{rDisplay,1}     = sprintf('<html><i><b>%s</b></i></html>', obj.names{r});
                            table_data{rDisplay,2}     = sprintf('<html><i><b>- - - -</b></i></html>');
                            table_data{rDisplay,3}     = sprintf('<html><i><b>%s</b></i></html>', obj.names{r});

                        otherwise
                            table_data{rDisplay,1}     = obj.names{r};
                            table_data{rDisplay,2}     = obj.values{r};
                            table_data{rDisplay,3}     = sprintf('<html><i>%s</i></html>', obj.descriptions{r});
                    end
                end
            end
            
            % Set it in the GUI
            set(obj.table_handle, 'Data', table_data);
        end
        
        function mode = getMode( obj, name )
            %% Return the value of the specified parameter
            for r = 1:obj.Nrows
                if(strcmp(name,obj.names{r}))
                    mode = obj.modes{r};
                    return
                end                
            end
            error('Parameter %s name not found', name);
        end
        
        function value = getValue( obj, name )
            %% Return the value of the specified parameter
            for r = 1:obj.Nrows
                if(strcmp(name,obj.names{r}))
                    value = obj.values{r};
                    return
                end                
            end
            error('Parameter %s name not found', name);
        end
        
        function table_dataStructure = getDataStructure(obj)
            %% Build a data structure of all of the elements of the table
            % Ignore the HEADER types since those don't contain information
            % of interest
            %  Each row from the table is turned into a parameter in this structure
            %  Spaces ( ) in row names are converted to undersctores (_) in paramter
            %  names
            
            table_dataStructure = [];
            
            for r = 1:obj.Nrows
                % Save the information as long as the row is not HEADER
                if( ~strcmpi( obj.modes{r}, 'HEADER') )
                   
                    % Assign the value to a field with the name of the row
                    variableName = obj.names{r};
                    
                    % Convert spaces to undersctore
                    variableName = strrep(variableName, ' ', '_');
                    
                    % Execute the assignment
                    try
                        eval( sprintf('table_dataStructure.%s = obj.values{r};', variableName) );
                    catch err
                        fprintf('\n\nERROR: table_dataStructure.%s = obj.values{r};', variableName);
                        fprintf('\n\nERROR: make sure your variable name (%s) does not have invalid characters: . - + / \ ^ and others' ,variableName);
                        fprintf('\n\n');
                        rethrow(err);
                    end
                        
                end
            end
        end
        
        function readTable( obj )
            %% Read data from table into properties of the class            
            table_data = get(obj.table_handle, 'Data');
            
            rDisplay = 0;
            for r = 1:obj.Nrows
                if( obj.isVisible{r} )
                    rDisplay = rDisplay+1;
                    value = table_data{rDisplay,2};
                    
                    switch obj.modes{r}
                        case 'NUMERIC'
                            if( value < obj.values_Min{r} || value > obj.values_Max{r} )
                                message_string = sprintf('Bad input: %s = %f (valid input is within %f and %f', ...
                                    obj.names{r}, value, obj.values_Min{r}, obj.values_Max{r});
                                msgbox(message_string, 'Bad input', 'error');                                
                                return
                            end
                    end
                
                    % Set parameter
                    obj.values{r} = value;
                end
            end
        end
        
        function setVisibility(obj, name, isVisible)
            %% Set the visibility of a row
            for r = 1:obj.Nrows
                if(strcmp(name,obj.names{r}))
                    obj.isVisible{r} = isVisible;
                    return
                end                
            end
            error('Parameter %s name not found', name);            
        end
        
        function writeLogFile( obj, OUTPUTSTREAM )
            %% Append the contents and values of the table to a log file
            %  Format this nicely for easy viewing
            % The OUTPUTSTREAM must be already open
            
            for r = 1:obj.Nrows
                switch obj.modes{r}
                    case 'HEADING'
                        % Add an extra newline and format at a string %s
                        fprintf(OUTPUTSTREAM, '\n');
                        fprintf(OUTPUTSTREAM, '\n%s', obj.names{r});
                        
                    case {'NUMERIC','BOOLEAN'}
                        % Format these as a number %f
                        fprintf(OUTPUTSTREAM, '\n%s\t%f\t%s', ...
                            obj.names{r}, obj.values{r}, obj.descriptions{r});
                        
                        
                    case {'FOLDER', 'FILE', 'STRING'}
                        % Format these as a string %s
                        fprintf(OUTPUTSTREAM, '\n%s\t%s\t%s', ...
                            obj.names{r}, obj.values{r}, obj.descriptions{r});
                        
                    otherwise, 
                        error('Invalid format "%s". See code', obj.modes{r});
                        
                end
            end                
        end
            
        
        function removeRow(obj, r)
            
        end
    end    
end