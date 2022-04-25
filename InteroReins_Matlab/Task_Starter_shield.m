% visdiff('/Users/chendanlei/Google Drive/Interoceptive Rinstatement/IanKleckner_Study_Tasks-2019.08.28/Task_Starter.m', '/Users/chendanlei/Google Drive/Interoceptive Rinstatement/IanKleckener_Study_Tasks-2019.05.15/InteroAffect_TaskStarter.m')
% see here for gui creation from MathWorks: https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/5506/versions/1/previews/gui_mainfcn.m/index.html


% Ian Kleckner
% ian.kleckner@gmail.com
% Cancer Control Division
% University of Rochester Medical Center
%
% Interoception / Affect Study launcher program
%
% Goal: Specify basic settings, and call one of several programs 
%
% 2012/01/30 Start coding
%            Copy code from Ian Kleckner's CFS Initializer [2011/03/13 ->]
% 2012/02/07 Update to start one of many Intero/Affect tasks
% 2012/05/08 Update PPnumber to PPnumber_in_group (each group is
%               counter-balanced within itself)
%            Change HBD_training_expt -> HeartbeatDetection_expt
% 2012/07/04 Update supporing program - Picture Rating Form
% 2012/07/06 Update trigger lists for Affective Judgment task
%            Include Move_if_you_Want, and Rate-Val and Rate-Aro
% 2012/07/24 Update instructions for piture rating form
% 2013/05/20 Update for new Interoception/Affect study
% 2016/12/23 Update for Kleckner R25 Pilot at URMC
% 2019/06/19 Update for Kleckner K07 study at URMC


function varargout = Task_Starter(varargin)
% Task_Starter M-file for Task_Starter.fig
%      Task_Starter, by itself, creates a new Task_Starter or raises the existing
%      singleton*.
%
%      H = Task_Starter returns the handle to a new Task_Starter or the handle to
%      the existing singleton*.
%
%      Task_Starter('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in Task_Starter.M with the given input arguments.
%
%      Task_Starter('Property','Value',...) creates a new Task_Starter or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Task_Starter_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Task_Starter_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Task_Starter

% Last Modified by GUIDE v2.5 19-Jun-2019 18:31:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Task_Starter_OpeningFcn, ...
                   'gui_OutputFcn',  @Task_Starter_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before Task_Starter is made visible.
function Task_Starter_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Task_Starter (see VARARGIN)

% Choose default command line output for Task_Starter
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Task_Starter wait for user response (see UIRESUME)
% uiwait(handles.GUI_task_starter);

%% Set up the tasks that can be started
set(handles.popupmenu_task, 'String', ...
    { 'View_and_Record_Shimmer', ...
      'HeartbeatTracking_expt_Shimmer', ...
      'HeartbeatDetection_expt_Shimmer'} ...
    );

%% Input specs table
SoftwareTitle   = 'Physiological Tasks Starter';
SoftwareVersion = 'v.2019.06.25';

% Commit the software title to the GUI
set(handles.text_title, 'String', SoftwareTitle);
set(handles.GUI_task_starter, 'Name', sprintf('[%s]',SoftwareTitle));
setappdata(handles.GUI_task_starter, 'SoftwareTitle', SoftwareTitle);

% Commit the software version to the GUI and save it for future use
setappdata(handles.GUI_task_starter, 'SoftwareVersion', SoftwareVersion);
set(handles.text_version, 'String', sprintf('Ian Kleckner, URCC\n%s ', SoftwareVersion));

specsTable = DynamicTable(handles.table_input_specs);
% Syntax for table
% specsTable.addRow(mode, name, value, description, value_Min, value_Max)

specsTable.addRow('HEADING', 'Participant Information', NaN, '', NaN, NaN);    
    specsTable.addRow('NUMERIC', 'PP_Number', 999, '** Number in the study (1, 2, 3, ...)', 1, Inf);
    specsTable.addRow('NUMERIC', 'Timepoint', 1, '** 1, 2, or 3', 1, 3);
    %specsTable.addRow('FILE', 'StimOrder_PriorSession', 'SELECT FILE', '** AJ ONLY - StimOrder.mat file from session ONE for this PP and this task', NaN, NaN);

specsTable.addRow('HEADING', 'Shimmer Connection', NaN, '', NaN, NaN);
    specsTable.addRow('BOOLEAN', 'Use_ECG', true, '** Acquire ECG data', false, true);    
    specsTable.addRow('NUMERIC', 'COM_ECG', 8, 'COM port for Shimmer with ECG', 0, Inf);
    specsTable.addRow('NUMERIC', 'Sampling_rate_ECG_Hz', 256, 'Sampling rate for ECG in Hz for NON-HBD tasks', 0, Inf);
    
    specsTable.addRow('BOOLEAN', 'Use_EDA', false, '** Acquire EDA data', false, true);    
    specsTable.addRow('NUMERIC', 'COM_EDA', 9, 'COM port for Shimmer with EDA', 0, Inf);
    specsTable.addRow('NUMERIC', 'Sampling_rate_EDA_Hz', 32, 'Sampling rate for EDA in Hz', 0, Inf);

specsTable.addRow('HEADING', 'ECG R Spike Detection', NaN, '', NaN, NaN);
    specsTable.addRow('NUMERIC', 'Minimum_RR_Interval_sec', 0.5, 'Minimum duration between consecutive ECG R spikes (sec)', 0, Inf);
    specsTable.addRow('NUMERIC', 'Minimum_R_Prominence_mV', 0.2, '** Minimum height of ECG R spike(mV)', 0, Inf);

specsTable.addRow('HEADING', 'View and Record Shimmer', NaN, '', NaN, NaN);
    specsTable.addRow('STRING', 'Event_name', 'Test', '** Name of the event (Test, 6MWT, Biodex)', NaN, NaN);
    specsTable.addRow('NUMERIC', 'Plot_update_period_sec', 0.5, 'Time between plot updates (sec)', 1e-3, Inf);
    specsTable.addRow('NUMERIC', 'Plot_viewable_duration_sec', 15, 'Visible region of plot (sec)', 1, Inf);
    specsTable.addRow('NUMERIC', 'Test_duration_sec', Inf, 'Duration of test (sec)', 0, Inf);
    specsTable.addRow('BOOLEAN', 'Write_data_to_file', true, 'Write Shimmer data to text file', 0, Inf);
  
specsTable.addRow('HEADING', 'Heartbeat Detection', NaN, '', NaN, NaN);
    specsTable.addRow('NUMERIC', 'Number_of_trials_HBD', 60, 'Number of trials', 1, Inf);
    specsTable.addRow('NUMERIC', 'Sampling_rate_ECG_HBD_Hz', 1024, 'Sampling rate for ECG in Hz for HBD only', 0, Inf);    
    specsTable.addRow('BOOLEAN', 'TrainingMode', false, 'Present Knowledge of Results (KOR) after each trial', false, true);    
    specsTable.addRow('BOOLEAN', 'ShowProgress', true, 'Show PP their progress at 25%, 50%, and 75% done', false, true);
    
specsTable.addRow('HEADING', 'Advanced Display', NaN, '', NaN, NaN);
    specsTable.addRow('BOOLEAN', 'SpeedMode', false, 'Skip program instructions for rapid execution (debugging)', false, true);

specsTable.addRow('HEADING', 'Display', NaN, '', NaN, NaN);    
    specsTable.addRow('BOOLEAN', 'FullScreenMode', true, 'Display in full-screen mode', false, true);
        specsTable.addRow('NUMERIC', 'WindowPixels_Width', 1200, 'Width of screen (if not fullscreen; Pixels)', 1, Inf);
        specsTable.addRow('NUMERIC', 'WindowPixels_Height', 800, 'Height of screen (if not fullscreen; Pixels)', 1, Inf);    
    specsTable.addRow('BOOLEAN', 'HideMousePointer', true, 'Hide the mouse pointer', false, true);
    specsTable.addRow('NUMERIC', 'MinBorderPercent', 1, 'Percent window width or height (smaller of the two) for IMAGE border (1-49%)', 1, 49);    
    
if( ismac )
    SYNCTEST_DEFAULT = true;
else
    %SYNCTEST_DEFAULT = false;
    SYNCTEST_DEFAULT = true;
end

specsTable.addRow('BOOLEAN', 'SetSyncTest_StDev', SYNCTEST_DEFAULT, 'MacOS=ON, Manually set VBL SyncTest max standard deviation (below)', false, true);
    specsTable.addRow('NUMERIC', 'VBL_MaxStd_ms', 5, 'Maximum allowable stdev in monitor refresh interval (5 ms for "relaxed" to fix issue on Mac)', 0.0001, Inf);

specsTable.displayTable();
setappdata(handles.GUI_task_starter, 'specsTable', specsTable);

% Current directory (to address a bug)
setappdata(handles.GUI_task_starter, 'currentDirectory', pwd);

clc;
refresh_display( handles )


% --- Outputs from this function are returned to the command line.
function varargout = Task_Starter_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function refresh_display( handles )
%% Update graphical elements

% Specifications table
specsTable = getappdata(handles.GUI_task_starter, 'specsTable');


% Optional display elements

% For loading the prior stimulus order
%if( specsTable.getValue('SessionNumber') ~= 1 )
%    specsTable.setVisibility('StimOrder_PriorSession', true);
%else
%    specsTable.setVisibility('StimOrder_PriorSession', false);
%end

if( specsTable.getValue('FullScreenMode') == 1 )
    specsTable.setVisibility('WindowPixels_Width', false);
    specsTable.setVisibility('WindowPixels_Height', false);
else
    specsTable.setVisibility('WindowPixels_Width', true);
    specsTable.setVisibility('WindowPixels_Height', true);
end


% % Optional display elements
% if( specsTable.getValue('CalibrateContrast') )
% %    specsTable.setVisibility('NumPracticeTrials', true);
%     specsTable.setVisibility('NumCalibTrials', true);
%     specsTable.setVisibility('NumFeedbackTrials_Calib', true);
%     specsTable.setVisibility('NumTrialsPerReset', true);    
%     specsTable.setVisibility('CalibrationFilename', false);
% else
% %    specsTable.setVisibility('NumPracticeTrials', false);
%     specsTable.setVisibility('NumCalibTrials', false);
%     specsTable.setVisibility('NumFeedbackTrials_Calib', false);
%     specsTable.setVisibility('NumTrialsPerReset', false);    
%     specsTable.setVisibility('CalibrationFilename', true);
% end


% % Optional display elements
% if( specsTable.getValue('CalculateTrialContrast') )    
%     specsTable.setVisibility('DesiredProbCorrect', true);
%     specsTable.setVisibility('LogTrialContrast', false);
% else    
%     specsTable.setVisibility('DesiredProbCorrect', false);
%     specsTable.setVisibility('LogTrialContrast', true);
% end


% For advanced display options
if( specsTable.getValue('SetSyncTest_StDev') )
    specsTable.setVisibility('VBL_MaxStd_ms', true);
else
    specsTable.setVisibility('VBL_MaxStd_ms', false);
end


% Display the contents of the table
specsTable.displayTable();


% --- Executes on button press in button_GO.
function button_GO_Callback(hObject, eventdata, handles)
% hObject    handle to button_GO (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Start the Task_Starter program

% Read a special data structure that can be passed to the experimental
% program
%  Each row from the table is turned into a parameter in this structure
%  Spaces ( ) in row names are converted to undersctores (_) in paramter
%  names
specsTable          = getappdata(handles.GUI_task_starter, 'specsTable');
table_dataStructure = specsTable.getDataStructure();

table_dataStructure.currentDirectory          = getappdata(handles.GUI_task_starter, 'currentDirectory');
table_dataStructure.SoftwareVersion           = getappdata(handles.GUI_task_starter, 'SoftwareVersion');
table_dataStructure.SoftwareTitle             = getappdata(handles.GUI_task_starter, 'SoftwareTitle');

% To pass the contents of the table to the experiment
%  This allows for the log file to be written dynamically
table_dataStructure.dynamicTable = specsTable;

% Disable the GO button so that it is not accidentally activated during
% runtime -- this usually happens if the user hits the SPACE bar
set(handles.button_GO, 'Enable', 'off');
set(handles.button_closefigs, 'Enable', 'off');

% Determine which task should be started
taskName_array  = get(handles.popupmenu_task, 'String');
taskNumber      = get(handles.popupmenu_task, 'Value');
taskName        = taskName_array{ taskNumber };

% Task names are something like this
%  HeartbeatDetection_expt
%  EvocativePhotos_expt
%  AffectiveMisattribution_expt

% The eval command will call the name of the task and provide it with the
% information from the dynamic table
eval( sprintf( '%s( table_dataStructure )', taskName ) )

% This is an example of an explicit call: HeartbeatDetection_expt( table_dataStructure )

% Turn the GO button back on again
set(handles.button_GO, 'Enable', 'on');
set(handles.button_closefigs, 'Enable', 'on');



% --- Executes when entered data in editable cell(s) in table_input_specs.
function table_input_specs_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to table_input_specs (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
%% Data in the table has changed
specsTable = getappdata(handles.GUI_task_starter, 'specsTable');
specsTable.readTable();


refresh_display(handles);

% --- Executes when selected cell(s) is changed in table_input_specs.
function table_input_specs_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to table_input_specs (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
%% Data in the table has been clicked
% 2011/09/22 Update for getting FILE (not just FOLDER)

% Check if something has been clicked
if( ~isempty(eventdata.Indices) )
    
    % Read what has been clicked in the table
    row         = eventdata.Indices(1);
    table_data  = get(handles.table_input_specs, 'Data');
    paramName   = table_data{row,1};
    
    % Check what type of parameter this is, from the dynamic table
    specsTable  = getappdata(handles.GUI_task_starter, 'specsTable');
    
    try
        mode        = specsTable.getMode(paramName);
    catch
        % Parameter not found
        return
    end

    % If it is a folder, then select via GUI window
    if( strcmp(mode, 'FOLDER') )
        directoryName = uigetdir('.', 'Select folder...');
        

        % If the path WAS selected...
        if( ~isequal(directoryName,0) )
            % ...then commit the change
            table_data{row, 2} = directoryName;
            set(handles.table_input_specs, 'Data', table_data);
            specsTable.readTable();
            refresh_display(handles);
        end
        
    elseif( strcmp(mode, 'FILE') )
        [filename, filepath] = uigetfile('*.*', 'Select file...');
       
        % If the file WAS selected...
        if( ~isequal(filename,0) )
            % ...then commit the change            
            table_data{row, 2} = sprintf('%s/%s', filepath, filename);
            set(handles.table_input_specs, 'Data', table_data);
            specsTable.readTable();
            refresh_display(handles);
        end
    end
end


% --- Executes on button press in button_clear.
function button_clear_Callback(hObject, eventdata, handles)
% hObject    handle to button_clear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Turn the GO button on in case it was off
set(handles.button_GO, 'Enable', 'on');
set(handles.button_closefigs, 'Enable', 'on');

% Clear the workspace and variables
clc;
clear;
clear all;
%clear classes;


% --- Executes on selection change in popupmenu_task.
function popupmenu_task_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_task (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu_task contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_task


% --- Executes during object creation, after setting all properties.
function popupmenu_task_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_task (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in button_closefigs.
function button_closefigs_Callback(hObject, eventdata, handles)
% hObject    handle to button_closefigs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figHandles = findobj('Type', 'figure');
Nfigs_closed = 0;

for f = 1:length(figHandles)    
    %figHandles(f).Name    
    if( ~contains(figHandles(f).Name, 'Starter') )
        close(figHandles(f));
        Nfigs_closed = Nfigs_closed + 1;
    end    
end
fprintf('\nClosed %d figures\n', Nfigs_closed);
