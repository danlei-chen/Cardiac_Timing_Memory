% Ian Kleckner
% Interdisciplinary Affective Science Lab (IASL)
%
% Sandbox for testing various program functions
%
% 2012/01/21 Heartbeat Detection Trainingtesting
%

%% ECG basleine correct [2019/05/02]

% h = findobj(gca,'Type','line')
% x=get(h,'Xdata')
% y=get(h,'Ydata')
% 
% X = x{2};
% Y = y{2};
% 
% X_peaks = x{1};
% Y_peaks = y{2};
load('data/data_ECG_baseline-2019.05.03.mat', 'X', 'Y', 'X_peaks', 'Y_peaks');

sampling_period_ECG_sec = X(2)-X(1);


Minimum_RR_Interval_sec = 0.5;
Minimum_R_Prominence_mV = 1;

ISOELECTRIC_TIME_BEFORE_RSPIKE_SEC = 66e-3;
ISOELECTRIC_SAMPLES_BEFORE_RSPIKE = floor(ISOELECTRIC_TIME_BEFORE_RSPIKE_SEC / sampling_period_ECG_sec);



[peak_Y_array, peak_X_array, peak_width_array, peak_prom_array] = ...
                findpeaks(Y, X, ...
                'MinPeakDistance', Minimum_RR_Interval_sec, ...
                'MinPeakProminence', Minimum_R_Prominence_mV );
            
X_for_spline = peak_X_array - ISOELECTRIC_TIME_BEFORE_RSPIKE_SEC;
Y_for_spline = NaN * X_for_spline;

for p = 1:length(X_for_spline)
    k_this_peak     = find( X == peak_X_array(p) );
    Y_for_spline(p) =  Y( k_this_peak - ISOELECTRIC_SAMPLES_BEFORE_RSPIKE );
end

X_spline = X;
Y_spline = spline(X_for_spline, Y_for_spline, X);

plot(X, Y);
hold all;
plot( X_spline, Y_spline, '-r');

Y_adj = Y - Y_spline;

plot(X, Y_adj, '-b');

            
            




%% ECG filter [2019/04/30]
clc;
clear all;

load('data/data_noisy_ECG-2019.04.30.mat');

%--------------------------------------------------------------------------
% Power spectrum
plot(X_ECG, Y_ECG)

Y_FFT = fft(Y_ECG);

% Sampling frequency
fs = 1 / (X_ECG(2) - X_ECG(1));

%hist(Y_FFT);
% 
% L = length(Y_ECG);
% 
% P2 = abs(Y_FFT/L);
% P1 = P2(1:L)/2+1;
% P1(2:end-1) = 2*P1(2:end-1);
% 
% f = sampling_frequency*(0:(L/2))/L;
% plot(f,P1) 
% 
% 

n = length(Y_ECG);          % number of samples
f = (0:n-1)*(fs/n);     % frequency range
power = abs(Y_FFT).^2/n;    % power of the DFT

semilogy(f,power)

%--------------------------------------------------------------------------
% Attempt to filter






%% Names of accelerometer channels [2017/05/10]


comPort_ECG = '18';
sampling_rate_ECG_Hz = 512;

% Connect to Shimmer at specified COM port
fprintf('\n\nAttempting to connect to Shimmer, this may take 10 sec...\n');
shimmer = ShimmerHandleClass(comPort_ECG);

% assign user friendly macros for setenabledsensors
SensorMacros = SetEnabledSensorsMacrosClass;                               

% Ensure connection
if (shimmer.connect)

    fprintf('\n\nSuccessfully connected to Shimmer device!\n');

    % Select sampling rate
    shimmer.setsamplingrate(sampling_rate_ECG_Hz);                                           

    % Select internal expansion board; select 'ECG' to enable both SENSOR_EXG1 and SENSOR_EXG2 
    shimmer.setinternalboard('ECG');                                       

    % Disable other sensors
    shimmer.disableallsensors;                                             

    % Enable SENSOR_EXG1 and SENSOR_EXG2
    shimmer.setenabledsensors(SensorMacros.ECG,1, SensorMacros.ACCEL, 1);

else
    error('Could not connect to Shimmer device');
end

%%

shimmer.start;
WaitSecs(3);

[newData, signalNameArray, signalFormatArray, signalUnitArray] = shimmer.getdata('c');

shimmer.stop;


%%

signalNameArray_ECG = signalNameArray;

k_accelX_ECG = find(strcmp('Low Noise Accelerometer X', signalNameArray_ECG))
k_accelY_ECG = find(strcmp('Low Noise Accelerometer Y', signalNameArray_ECG))
k_accelZ_ECG = find(strcmp('Low Noise Accelerometer Z', signalNameArray_ECG))





%% Part I: TIming of shimmer start and stop 2017/01/01

comPort_ECG = '18';
sampling_rate_ECG_Hz = 512;

% Connect to Shimmer at specified COM port
fprintf('\n\nAttempting to connect to Shimmer, this may take 10 sec...\n');
shimmer = ShimmerHandleClass(comPort_ECG);

% assign user friendly macros for setenabledsensors
SensorMacros = SetEnabledSensorsMacrosClass;                               

% Ensure connection
if (shimmer.connect)

    fprintf('\n\nSuccessfully connected to Shimmer device!\n');

    % Select sampling rate
    shimmer.setsamplingrate(sampling_rate_ECG_Hz);                                           

    % Select internal expansion board; select 'ECG' to enable both SENSOR_EXG1 and SENSOR_EXG2 
    shimmer.setinternalboard('ECG');                                       

    % Disable other sensors
    shimmer.disableallsensors;                                             

    % Enable SENSOR_EXG1 and SENSOR_EXG2
    shimmer.setenabledsensors(SensorMacros.ECG,1);

else
    error('Could not connect to Shimmer device');
end

%% Part II: TIming of shimmer start and stop 2017/01/01

% First need to run code above to set up shimmer device

duration_desired_sec = 10;

%--------------------------------------------------------------------------
shimmer.stop;

time_before_start_sec = GetSecs();
shimmer.start;
time_after_start_sec = GetSecs();

WaitSecs(duration_desired_sec);

duration_before_start_to_before_pull_sec = GetSecs() - time_before_start_sec;
duration_after_start_to_before_pull_sec = GetSecs() - time_after_start_sec;

[newData, signalNameArray, signalFormatArray, signalUnitArray] = shimmer.getdata('c');

duration_before_start_to_after_pull_sec = GetSecs() - time_before_start_sec;
duration_after_start_to_after_pull_sec = GetSecs() - time_after_start_sec;

% Calculate elapsed time in actual data
k_time = find(strcmp('Time Stamp', signalNameArray));
sampling_period_sec = (newData(2, k_time) - newData(1, k_time)) / 1e3;
duration_data_sec = length(newData(:,1)) * sampling_period_sec;

% Doesn't matter as much when this occurs
shimmer.stop;

fprintf('\n\n%f\tDesired duration of data', duration_desired_sec);
fprintf('\n%f\tActual duration of data', duration_data_sec);
fprintf('\n%f\tBefore start to before pull', duration_before_start_to_before_pull_sec);
fprintf('\n%f\tBefore start to after pull', duration_before_start_to_after_pull_sec);
fprintf('\n%f\tAfter start to before pull', duration_after_start_to_before_pull_sec);
fprintf('\n%f\tAfter start to after pull', duration_after_start_to_after_pull_sec);


%% Test PTB keyboard input [2017/01/01]

%reply=Ask(window,'Who are you?',[],[],'GetChar',RectLeft,RectTop); % Accept keyboard input, echo it to screen.


% Create a DisplaySession
DS = DisplaySession;

% The minimum border percentage for scaling displayed images (<50)
DS.setMinBorderPercent( 0.1 );

X_left = 0;
Y_upper = 0;
X_width = 1200;
Y_height = 800;

rect_windowed = [X_left, Y_upper, X_left+X_width, Y_upper+Y_height];

rect_windowed = [];


% Stereo-mode as per the PTB documentation (see below)
%  0 -> no stereo, 4 -> split left/right
stereoMode = 0;

% Create the PTB window for display
[window, rect] = DS.createWindow(rect_windowed, stereoMode);


text_string = 'How many heartbeats\ndid you count?\n\nType number, hit Enter:';
text_x = (rect(3) - rect(1)) / 8;
text_y = 0.2 * rect(4);

question_string = '';
question_x = (rect(3) - rect(1)) / 2.5;
question_y = 0.7 * rect(4);

quesstion_fontsize = 250;
oldTextSize = Screen('TextSize', window, quesstion_fontsize);

[nx, ny, textbounds] = DrawFormattedText(window, text_string, 'center', text_y, 255*[1 1 1]);

[string,terminatorChar] = GetEchoString(window, question_string, question_x, question_y, 255*[1 1 1], [0 0 0] );
%[string,terminatorChar] = GetEchoString(window, msg, x, y);

string
terminatorChar

% Accept keyboard input, echo it to screen.
%reply = Ask(window, 'How many heartbeats did you count?', 255*[1 1 1], [.5 .5 .5], 'GetChar', 'top', 'bottom', 80);
%reply = Ask(window, 'How many heartbeats did you count? ', 255*[1 1 1], [.5 .5 .5], 'GetChar', [], [], 40);
%]reply

Screen('CloseAll')





%% Test randomizing faces [2012/05/02]

clc

for pp = 1:10
    fprintf('\n\n==========================================================');
    fprintf('\n==========================================================');
    fprintf('\nPPNUMBER = %d', pp);
    fprintf('\n\n');

    experiment_specs.USE_CROPPED_FACE_FOR_SUPPRESSED_SIDE = true;
    experiment_specs.PPnumber = pp;
    experiment_specs.SessionNumber = 1;

    % Constraints for this section
    %  * Images are named {M,F}-###-{NEUTRAL,ANGRY,HAPPY}-{CROP,FULL}.jpg
    %
    %  * For a given PP, each face is shown with the same suppressed emotion
    %  across all judgments
    %
    %  * Across PPs, each face is assigned with each suppressed emotion
    %
    %  * Each PP has two sessions
        
    %----------------------------------------------------------------------    
    % Options
    
    % Folder name for the face images
    foldername_faces = 'data/AffectiveMisattribution/faces';
    
    % This NEEDS to match the contents of the folder  
    %  Num of IDs must be divisible by 3 (for expressions) and by 4 (to
    %  divide into two sessions and then into halves of each session)
    Nfemale_IDs     = 24;
    Nmale_IDs       = 24;
    
    if( rem(Nfemale_IDs,3) ~= 0 || rem(Nfemale_IDs,4) ~= 0 || ...
        rem(Nmale_IDs,3) ~= 0 || rem(Nmale_IDs,4) ~= 0  || ...
        Nfemale_IDs ~= Nmale_IDs )
        error('Num of IDs must be divisible by 3 (for expressions) and by 4 AND Nfemale = Nmale');
    end
    
    % Number of trials (each trial is a unique face ID)
    %  Half the males and half the females into each of 2 sessions
    Ntrials_per_session = Nfemale_IDs/2 + Nmale_IDs/2;
    
    % Number of sessions
    Nsessions = 2;
    
    % Facial experessions per identity (ANGRY HAPPY NEUTRAL)
    faceExpressionString_array = {'ANGRY', 'HAPPY', 'NEUTRAL'};    
    NfaceExpressions = length( faceExpressionString_array );
    
    % Use a CROPPED face for suppressed image
    %  if false, use a FULL face for suppressed image
    USE_CROPPED_FACE_FOR_SUPPRESSED_SIDE = experiment_specs.USE_CROPPED_FACE_FOR_SUPPRESSED_SIDE;
    %----------------------------------------------------------------------    
    
    % Browse the directory where the face images reside
    dirStruct_female    = dir(sprintf('%s/%s', foldername_faces, 'F-*-NEUTRAL-FULL.jpg'));
    dirStruct_male      = dir(sprintf('%s/%s', foldername_faces, 'M-*-NEUTRAL-FULL.jpg'));
    
    % Check to see how many identities are in the directory
    Nfemale_files       = length(dirStruct_female);
    Nmale_files         = length(dirStruct_male);
    
    % Throw an error if there is not the expected file structure
    if( Nfemale_files ~= Nfemale_IDs )
        error('File number mismatch');
    end
    
    if( Nmale_files ~= Nmale_IDs )
        error('File number mismatch');
    end
    
    
    % Get face ID numbers from the filename
    for f = 1:Nfemale_files
        filename            = dirStruct_female(f).name;
        femaleID_array(f)   = sscanf(filename, 'F-%03d-NEUTRAL-FULL.jpg');
    end
    
    for f = 1:Nmale_files
        filename            = dirStruct_male(f).name;
        maleID_array(f)   = sscanf(filename, 'M-%03d-NEUTRAL-FULL.jpg');
    end
    
    % Get session number from GUI
    sessionNum = experiment_specs.SessionNumber;
    
    if( sessionNum == 1 )
        % Generate a new trial order for this PP
        %
        % GOALS
        % * Across PPs, each ID is paired with each expression the SAME
        % number of times
        %   E.g., F-001 is ANGRY, HAPPY, and NEUTRAL an equal number of
        %   times across PPs
        %   I.e., F-001 is ANGRY for PPs 1, 4, 7, 10, ...
        %         F-001 is HAPPY for PPs 2, 5, 8, 11, ...
        %         F-001 is NEUTRAL for PPs 3, 6, 9, 12, ...
        %
        % * For each PP, each session has the same number of each
        % expression
        %    E.g., each session contains 10 ANGRY, 10 HAPPY, and 10 NEUTRAL
        %    faces
        %
        % * For each PP, randomize order of overt gender
        %    E.g., M F F M M M F F F  ...
        %
        % * For each PP, the faces in session one and session two are
        % different
        %
        % STEPS
        % * Assign each ID an expression (counter-balanced across PPs)
        % * Divide each ID-Expression for that PP into their session 1 or 2
        % (random for each PP)
        % * Randomize the order of each session
        %
        % TRIAL STRUCTURE
        %  Overt face
        %   Gender:(Male or Female) with Face:(Neutral)
        %
        %  Suppressed face
        %   Gender:(Opposite of overt) with Face:(Smile or Scowl or Neutral)
        
        % Number of face IDs that go to each expression
        %  E.g., 30 faces and 3 expresions => 10 IDs that are HAPPY, or
        %  ANGRY, or NEUTRAL
        NIDs_per_expression = Nfemale_IDs / NfaceExpressions;
        
        % Get participant number from GUI
        PPnumber = experiment_specs.PPnumber;
        
        % PPtrialType is 1, 2, or 3
        %  It determines whether each face ID gets ANGRY, HAPPY, or NEUTRAL
        %PPtrialType = mod(PPnumber-1, NfaceExpressions)+1;
        
        %------------------------------------------------------------------
        % Determine whether each Expression should go to session 1 or
        % session 2
        %
        %  Each ROW is for a face expression (1, 2, 3) = (ANGRY, HAPPY, NEUTRAL)
        %  Each COLUMN is for a face ID that would make that expression (1, 2, 3, 4, ...)
        %  Each VALUE is 1 or 2 for the session number that the
        %   ID/Expression is displayed at
        expression_to_session = ones(NfaceExpressions,NIDs_per_expression);
        
        for ex = 1:NfaceExpressions
            % Half of the expressions set to session 2
            % (the other half are already set to session 1)
            expression_to_session(ex,1:NIDs_per_expression/2) = 2;
            
            % Shuffle the order
            expression_to_session(ex,:) = expression_to_session(ex, randperm(NIDs_per_expression));
        end
        %------------------------------------------------------------------
        
        %------------------------------------------------------------------
        % Assign face expression to each female
        faceFilename_female_set_session{1} = Set('Female face filenames session 1');
        faceFilename_female_set_session{2} = Set('Female face filenames session 2');
        
        % Counts number of times each expression was used so that each
        % expression can be divided equally among the sessions
        expression_counter_array = zeros(3,1);
        
        for id = 1:Nfemale_IDs
            % Get the expression of this face: 1, 2, or 3
            %  ACross PPnumber, this makes a list 1, 2, 3, 1, 2, 3, 1, 2, 3, ...            
            expressionIndex = mod(id-1 + PPnumber-1, NfaceExpressions)+1;
            
            % The filename for the FULL expression
            filename = sprintf('%s-%03d-%s-FULL.jpg', 'F', femaleID_array(id), faceExpressionString_array{expressionIndex});
            
            % Count that this expression has been used
            %  E.g., this is the 3rd HAPPY face
            expression_counter_array(expressionIndex) = expression_counter_array(expressionIndex) + 1;
            
            % Determine which session it goes to (1 or 2)
            %  E.g., the 3rd HAPPY face goes to session 1
            sessionNum = expression_to_session(expressionIndex, expression_counter_array(expressionIndex));
            
            % Add it to the Set for the session
            faceFilename_female_set_session{sessionNum}.addElement(filename);
        end
        %------------------------------------------------------------------
        
        %------------------------------------------------------------------
        % Repeat this process for males
        
        % Determine whether each Expression should go to session 1 or
        % session 2
        expression_to_session = ones(NfaceExpressions,NIDs_per_expression);
        
        for ex = 1:NfaceExpressions
            % Half of the expressions set to 2 (the other half set to 1)
            expression_to_session(ex,1:NIDs_per_expression/2)    = 2;
            
            % Shuffle the order
            expression_to_session(ex,:) = expression_to_session(ex, randperm(NIDs_per_expression));
        end
        
        
        % Assign face expression to each male
        faceFilename_male_set_session{1} = Set('Male face filenames session 1');
        faceFilename_male_set_session{2} = Set('Male face filenames session 2');
        
        % Counts number of times each expression was used so that each
        % expression can be divided equally among the sessions
        expression_counter_array = zeros(3,1);
        
        for id = 1:Nmale_IDs
            % Get the expression 1, 2, 3, 1, 2, 3, 1, 2, 3, ...
            %  This is offset by PP trial type
            expressionIndex = mod(id-1 + PPnumber-1, NfaceExpressions)+1;
            
            % The filename for the FULL expression
            filename = sprintf('%s-%03d-%s-FULL.jpg', 'M', maleID_array(id), faceExpressionString_array{expressionIndex});
            
            % Count that this expression has been used
            expression_counter_array(expressionIndex) = expression_counter_array(expressionIndex) + 1;
            
            % Determine which session it goes to (1 or 2)
            sessionNum = expression_to_session(expressionIndex, expression_counter_array(expressionIndex));
            
            % Add it to the Set for the session
            faceFilename_male_set_session{sessionNum}.addElement(filename);
        end
        
        % Scamble the order of each set
        faceFilename_female_set_session{1}.resetRandomSample();
        faceFilename_female_set_session{2}.resetRandomSample();
        faceFilename_male_set_session{1}.resetRandomSample();
        faceFilename_male_set_session{2}.resetRandomSample();
        
        %------------------------------------------------------------------
        % Finally, prepare a GENDER order
        % GOAL - Show each face overt and each face suppressed
        %  Each face pair will be fixed, with identify A and B always shown
        %  together
        %
        %  First half of trials A is overt and B is suppressed with some
        %  expression
        %
        %  Second half of trials, B is overt and A is suppressed with some
        %  expression
        %
        %  This will have half 1's and half 0's for the first half of trials
        %  Then, that array is duplicated, inverted, and concatenated to
        %  the first array
        %
        % E.g., 1 0 1 0 is used to make 1 0 1 0 0 1 0 1
        
        
        % Prepare a random order for each session "s" (1 and 2)
        
        % Number of trials across all sessions = Nfemale_IDs + Nmale_IDs
        % Number of trials per session = Nfemale_IDs/2 + Nmale_IDs/2
        %  Each MALE and each FEMALE is shown as overt
        % But each ID is shown twice (once OVERT and once SUPPRESSED)
        
        Nids_per_session = Ntrials_per_session / 2;
        
        for s = 1:2
            % Prepare the first half of the trials
            overtFace_is_female_array_session{s} = ones(1, Nids_per_session);

            % Half of the first half of trials are male (i.e., 0)
            overtFace_is_female_array_session{s}(1:Nids_per_session/2) = 0;

            % Randomize the order of these trials
            overtFace_is_female_array_session{s} = overtFace_is_female_array_session{s}( randperm( Nids_per_session ) );            

            % Duplicate and invert this vector for the second half of all trials
            overtFace_is_female_array_part_2 = ~overtFace_is_female_array_session{s};

            % Add these trials to the end of the first set
            overtFace_is_female_array_session{s} = [overtFace_is_female_array_session{s} overtFace_is_female_array_part_2];
        end        
        %------------------------------------------------------------------
        
        
    elseif( sessionNum == 2 )
        
        % Load the image order
        fileOrder_struct = load( experiment_specs.StimOrder_PriorSession );
        
        % Report the loaded file order to the screen
        fprintf('\n\n\tLoaded file: %s', experiment_specs.StimOrder_PriorSession);
        
        % Get variables so they can be saved below
        faceFilename_set_session = fileOrder_struct.faceFilename_set_session;        
    else
        error('Bad session number');
    end
    
    %------------------------------------------------------------------
    % Show the presentation order for session 1 and for session 2
    for s = 1:2        
        fprintf('\n\n---------------------------------------------------');
        fprintf('\n\nPresentation stimuli for session %d', s);
        fprintf('\nTrialNum\tOvertFace\tSuppressedFace');
        
        % Assign the set to be used this session
        faceFilename_female_set     = faceFilename_female_set_session{s};
        faceFilename_male_set       = faceFilename_male_set_session{s};
        overtFace_is_female_array   = overtFace_is_female_array_session{s};
        
        % Determine total number of trials
        Ntrials = faceFilename_female_set.Nelements + faceFilename_male_set.Nelements;
        
        for t = 1:Ntrials

            % Get the stimuli for this presentation
            %  E.g., M-006-NEUTRAL-FULL.jpg, F-027-HAPPY-FULL.jpg                
            filename_full_female    = faceFilename_female_set.getElementRandomNoReplacementCircularReset();
            filename_full_male      = faceFilename_male_set.getElementRandomNoReplacementCircularReset();

            % Determine whether the overt face is female (or male)
            OVERT_FACE_IS_FEMALE    = overtFace_is_female_array(t);

            if( OVERT_FACE_IS_FEMALE )
                filenameOvert_full      = filename_full_female;
                filenameSuppressed_full = filename_full_male;

            else
                filenameOvert_full      = filename_full_male;
                filenameSuppressed_full = filename_full_female;                    
            end

            if( USE_CROPPED_FACE_FOR_SUPPRESSED_SIDE )
                % Get the cropped facename for the suppressed side
                filenameSuppressed = strrep(filenameSuppressed_full, 'FULL', 'CROP');  
            else
                filenameSuppressed = filenameSuppressed_full;
            end

            % Get the neutral facename for the overt side
            %  and make sure the face is NEUTRAL (find/replace any
            %  instance of HAPPY or ANGRY)
            filenameOvert = filenameOvert_full;
            filenameOvert = strrep(filenameOvert, 'HAPPY', 'NEUTRAL');
            filenameOvert = strrep(filenameOvert, 'ANGRY', 'NEUTRAL');

            % Save the filename without the directory on it
            filenameOvert_nodir         = filenameOvert;
            filenameSuppressed_nodir    = filenameSuppressed;
            
            % Report the trial details
            fprintf('\n%d\t%s\t%s', t, filenameOvert_nodir, filenameSuppressed_nodir);
        end
    end
    %------------------------------------------------------------------
end

return

%% Testing array randomizing stuff [2012/04/12]

clc
clear classes
clear all

Nfemale_IDs = 6;
s = 1;

% Prepare the first half of the trials
overtFace_is_female_array_session{s} = ones(1, Nfemale_IDs);
overtFace_is_female_array_session{s}

% Half of the first half of trials are male (i.e., 0)
overtFace_is_female_array_session{s}(1:Nfemale_IDs/2) = 0;
overtFace_is_female_array_session{s}

% Randomize the order of these trials
overtFace_is_female_array_session{s} = overtFace_is_female_array_session{s}( randperm( Nfemale_IDs ) );
overtFace_is_female_array_session{s}

% Duplicate and invert this vector for the second half of all trials
overtFace_is_female_array_part_2 = ~overtFace_is_female_array_session{s}

% Add these trials to the end of the first set
overtFace_is_female_array_session{s} = [overtFace_is_female_array_session{s} overtFace_is_female_array_part_2];
overtFace_is_female_array_session{s}



%% Testing SliderScale with stereo view [2012/04/10], [2012/04/19]

clc

clear all;
clear classes;

% Enable unified mode of KbName, so KbName accepts identical key names on
% all operating systems (KbDemo.m)
KbName('UnifyKeyNames');

% This script calls Psychtoolbox commands available only in OpenGL-based
% versions of the Psychtoolbox.
AssertOpenGL;

Screen('CloseAll');

% Create a DisplaySession
DS = DisplaySession;

% The minimum border percentage for scaling displayed images (<50)
DS.setMinBorderPercent( 30 );

DS.setEyeDominance( 'LEFT' )

X_left = 0;
Y_upper = 0;
X_width = 1280;
Y_height = 1024;

%rect_windowed = [X_left, Y_upper, X_left+X_width, Y_upper+Y_height];
rect_windowed = [];

% Allow the screen refresh rate to vary with stdev +/- maxStddev
%  This is to address an error in iMac screen setting
Screen('Preference', 'SyncTestSettings', 5/1000);

% Stereo-mode as per the PTB documentation (see below)
%  0 -> no stereo, 4 -> split left/right
stereoMode = 4;

% Create the PTB window for display
[window, rect] = DS.createWindow(rect_windowed, stereoMode);   

HideCursor();

% ---------------------------------------------------------------------
% Set the framing rectangle
DRAW_FRAMING_RECTANGLE          = true;

% (1) Add the default framing rectangle that is always on    
Frame_ColorRGB_border           = [255 255 255];
border_percent_X                = 25;
border_percent_Y                = 35;
PIXEL_WIDTH                     = 8;

frNum_Border = [];
frNum_Border = DS.setFramingRectangle( frNum_Border, ...
    DRAW_FRAMING_RECTANGLE, Frame_ColorRGB_border, border_percent_X, border_percent_Y, PIXEL_WIDTH );
% ---------------------------------------------------------------------

sliderScale_Judgment = SliderScale(DS);

% Show scale    
filename_rate_valence = 'data/AffectiveJudgment/SAM-V-9-Happy_on_Left.png';
sliderScale_Judgment.setDisplayImage(true, filename_rate_valence, false);            

% Do NOT show text on slide
sliderScale_Judgment.setDisplayText(false);    

% Set the resposne keys specified above
responseKey_commit  = 44;
responseKey_abort   = 41;
sliderScale_Judgment.setResponseKeys(responseKey_commit, responseKey_abort);
%----------------------------------------------------------------------

% Only add the following instructions if the program is not in speed
% mode
DS.addStep('NAVIGATE', 'STEREO', 'data/AffectiveJudgment/instructions/instructions-AJ-07.png', []);
DS.addNewSliderScale( sliderScale_Judgment );
DS.addStep('NAVIGATE', 'STEREO', 'data/AffectiveJudgment/instructions/instructions-AJ-12.png', []);        

DS.displayStep(1);

clear all;



%% More randomization [2012/03/20]

clc
clear all

Nfemale_IDs = 18;
NfaceExpressions = 3

% Determine whether each ID-Expression should go to session 1 or
% session 2        
for ex = 1:NfaceExpressions
    expression_to_session(ex,:)                  = zeros(1,Nfemale_IDs);
    expression_to_session(ex,1:Nfemale_IDs/2)    = 1;
    expression_to_session(ex,:)                  = expression_to_session(ex, randperm(Nfemale_IDs))+1
end

%% Test PP randomization [2012/03/19]

NfaceExpressions = 3;
Nidentities = 30;

PPnumber = 4;

clc



for id = 1:Nidentities
    % Get the expression 1, 2, 3, 1, 2, 3, 1, 2, 3, ...
    %
    %  This is offset by PP trial type
    expressionIndex = mod(id-1 + PPnumber-1, NfaceExpressions)+1;
    
    fprintf('\n%d\t%d', id, expressionIndex);

    % The filename for the FULL expression
    %filename = sprintf('%s-%03d-%s-FULL.jpg', 'F', femaleID_array(id), faceExpressionString_array{expressionIndex});

    % Add it to the Set
    %femaleFilename_set.addElement(filename);
end

%% Randomize the trial order with 50% of trials option1  [2012/03/19]

% The following is specified in the setup GUI
numTotalTrials_Calib    = 10;

% Generate list of trial types (1 or 2)
trialOption_array = ones(1,numTotalTrials_Calib)

if( mod(numTotalTrials_Calib,2) ~= 0 )
    error('Set an even number of trials');
end

% Change half of the elements into 2's 
trialOption_array(1:numTotalTrials_Calib/2) = 2

% Now shuffle the array of [2 2 2 2 2 ... 1 1 1 1 1]
random_index_order = randperm(numTotalTrials_Calib);

% Shuffle the array
trialOption_array = trialOption_array(random_index_order)
%------------------------------------------------------------------

%% Use Set class for blocks [2012/03/19]

clc
clear all
clear classes

%
% Constraints for this section
%  * Images need to be in folders with block names
%    E.g., data/AffectiveJudgment/Anchor/ where AffectiveJudgment is
%    the variable programName (above)
%  * There are only two sessions
%  * There must be exactly the number of images as required across both
%  sessions for each block
%    E.g., Nimages_in_block_array(2) = 10 requires taht block 2
%    (Negative images) must have 10*2 = 20 images to use across 2
%    sessions

%----------------------------------------------------------------------
% Options for this section

% Name of each block
block_name_array    = {'Anchor', ...            1
                       'Negative-Low', ...      2
                       'Negative-High', ...     3
                       'Neutral-Low', ...       4
                       'Positive-Low', ...      5
                       'Positive-High'}; ...    6

% Number of images per block
%  Anchor uses same images in both sessions
%  Other blocks must use half the total images in each session
%   i.e., 20 images for the Negative-Low block yields 10 images in that
%   block each session, and therefore write 10 below
Nimages_in_block_array  = [5 20 20 20 20 20];

% Total number of blocks
Nblocks             = length(block_name_array);

% The total number of sessions is not an option, but it's here FYI
Nsessions           = 2;

% Anchor, (pos/neg)-High, Neut-Low, (neg/pos)-High, (pos/neg)-Low, (neg/pos)-Low
block_order_options_session1 = { [1 6 4 3 5 2], ...
                                 [1 3 4 6 2 5] };

% Session 2 swaps the order of the two sets of blocks from session 1
block_order_options_session2 = { [1 3 4 6 2 5], ...
                                 [1 6 4 3 5 2] };

% Total number of possible orders (for counter-balancing PPs)
Nblock_order_options = size(block_order_options_session1, 2);
%----------------------------------------------------------------------

% Get session number from GUI
%sessionNum = experiment_specs.SessionNumber;

sessionNum = 1;
programName = 'AffectiveJudgment'
PPnumber = 1
input 'set session number. Ok?';

if( sessionNum == 1 )
    fprintf('\n\nReading files...');
    
    %----------------------------------------------------------------------
    % Create a set for each block
    % Read all files from each directory
    for b = 1:Nblocks
        blockSet_array{b} = Set( block_name_array{b} );        
        
        % Number of images in this block
        Nimages_in_block = Nimages_in_block_array(b);

        % Prepare the foldername
        foldername_block = sprintf('data/%s/%s/', programName, block_name_array{b});

        % Get the contents of the directory    
        fwf = FolderWithFiles(foldername_block);

        % Check if there are the proper number of files
        if( fwf.Nfiles ~= Nimages_in_block )
            error('There must be exactly enough images in the folder %s: %d (there are actually %d)', ...
                foldername_block, Nimages_in_block, fwf.Nfiles);
        end
        
        % Add each image to the block set
        fprintf('\n\n\tBlock %d: %s', b, foldername_block);            
        
        for f = 1:fwf.Nfiles
            % Read the filename
            filename    = fwf.getFilenameAbsolute( f );            
            
            % Add the filename to the set
            blockSet_array{b}.addElement(filename);
            
            fprintf('\n\tAdded file: %s', filename);
        end
    end
    %----------------------------------------------------------------------
    
    %----------------------------------------------------------------------
    % Divide the blocks into the two sessions
    
    for b = 1:Nblocks       
        
        if( strcmp(block_name_array{b}, 'Anchor') )
            % Anchor set is copied into each session
            blockSet_array_session1{b} = blockSet_array{b}.copy();
            blockSet_array_session2{b} = blockSet_array{b}.copy();
        
        else
            % Non-anchor sets are divided randomly into two parts
            dividedSet_array = blockSet_array{b}.divideSetRandomly(2);
            
            % Assign each division
            blockSet_array_session1{b} = dividedSet_array{1};
            blockSet_array_session2{b} = dividedSet_array{2};
        end
    end
    
    % Determine the index for the PP number
    %  Here, since there are 2 blocks, Odd PP#'s get block 1 first, 
    %  and even PP#s get block 2 first
    block_order_index = Set.getBlockNumber( PPnumber, Nblock_order_options );
    
    % Order of blocks for each session
    block_order_session1 = block_order_options_session1{ block_order_index };
    block_order_session2 = block_order_options_session2{ block_order_index };
    
    % Finally, prepare a blockSet for each of the two sessions
    %  so this can be used easily, in order
    blockSet_session1 = Set('Session 1');    
    
    % Go through each block in the order that they will appear in the task
    for b = block_order_session1
        blockSet_session1.addSubSet( blockSet_array_session1{b} );        
    end
    
    % Same for session 2
    blockSet_session2 = Set('Session 2');
    
    % Go through each block in the order that they will appear in the task
    for b = block_order_session2
        blockSet_session2.addSubSet( blockSet_array_session2{b} );
    end    
    %----------------------------------------------------------------------
    
    % Assign the set to be used this session
    blockSet_session = blockSet_session1;
    
elseif( sessionNum == 2 )
        % Load the image order
        fileOrder_struct = load( experiment_specs.StimOrder_PriorSession );
        
        % Set the current file order for session 2
        filename_array = fileOrder_struct.filename_array_session_2;
        
        % Report the loaded file order to the screen
        fprintf('\n\n\tLoaded file: %s', experiment_specs.StimOrder_PriorSession);
        
        % Get variables so they can be saved below
        blockSet_session1 = fileOrder_struct.blockSet_session1;
        blockSet_session2 = fileOrder_struct.blockSet_session2;
        
        % Assign the set to be used this session
        blockSet_session = blockSet_session2;
        
else
    error('Bad session number');
end

blockSet_session = blockSet_session2;
    
% Show the presentation order for session 1 and for session 2
Nblocks = blockSet_session.Nsubsets;

for b = 1:Nblocks

    % Get the current block in this session
    blockSet = blockSet_session.getSubSet(b);
    blockName = blockSet.name;

    fprintf('\nBlock %d / %d: %s', b, Nblocks, blockName);

    % Iterate through each image in the block
    Nimages = blockSet.Nelements;

    for im = 1:Nimages
        filename_image = blockSet.getElement(im);

        fprintf('\n\tImage %d / %d: %s', im, Nimages, filename_image);
    end
end


%% Test Set functions [2012/03/19]
clc
clear all
clear classes

set = Set('Test');
set.addElement('a');
set.addElement('b');
set.addElement('c');
set.addElement('d');
set.addElement('e');
set.addElement('f');

set.setPreventRepeatSample(true);

set

% Get subset

subsetArray = set.divideSetRandomly(2);

subsetArray{1}
subsetArray{2}


%% Testing PortIO [2012/03/12]

clear all;

% Options
filename_trigger_list = 'data/HeartbeatDetection/trigger_list-HBD-2012.03.12-MATLAB.txt';

% Load the file with the list of triggers
% portIO_out = PortIO('parallel', 'LPT1', 0:3, 'out');
% portIO_out.readTriggersFromFile(filename_trigger_list);
% portIO_out

% Create the communication to the parallel port
portIO_in = PortIO('parallel', 'LPT1', 4:7, 'in');


%% Wait for change [2012/03/12]
portIO_in.waitForPortChange();
portIO_out.findAndSendTrigger( 'HB--Tone' );    


%% Test reset [2012/03/12]

portIO_in.initialize()

%% Test triggers [2012/02/24]

clc
clear all
 
% Get into on installed adaptors
%hwinfo = daqhwinfo();
%hwinfo('parallel')
 
% Create a digital input/output object
dio = digitalio('parallel', 'LPT1');
%dio_in = digitalio('parallel', 'LPT1')
 
% Add a line group - one for input and one for output
% A line group consists of a mapping between hardware line IDs and MATLAB indices 
hwlines_out = addline(dio,0:7,'out')

% One-line sending a trigger
linevalue_array = [0 0 0 0 0 1 0 0]; putvalue(dio,linevalue_array); pause(1); putvalue(dio,0);

%% Check each value

for data = 0:255
    fprintf('\nWriting %03d', data);
	putvalue(dio,data);
    pause(0.25);
    putvalue(dio,0);
    pause(0.25);
end

%% Check each line

for bit = 0:7
    datum = 2^bit;
    linevalue_array = dec2binvec(datum, 8);
    
    fprintf('\nWriting value: ');
    disp(linevalue_array);    
    
	putvalue(dio,linevalue_array);
    pause(0.25);
    
    putvalue(dio,0);
    pause(0.25);
end


%% Testing Digital I/O [2012/02/20]

clc
clear all

% Get into on installed adaptors
%hwinfo = daqhwinfo();
%hwinfo('parallel')

% Create a digital input/output object
dio = digitalio('parallel', 'LPT1');
%dio_in = digitalio('parallel', 'LPT1')

% Add a line group - one for input and one for output
% A line group consists of a mapping between hardware line IDs and MATLAB indices 
hwlines_out = addline(dio,0:7,'in')
%hwlines_in = addline(dio_in,0:7,'in')

% Note: the "Port = 0" corresponds to the 8 data lines
% The other 17 lines on the parallel port are for other tasks
% 
% Index:  LineName:  HwLine:  Port:  Direction:  
%    1       'Pin2'     0        0      'In'        
%    2       'Pin3'     1        0      'In'        
%    3       'Pin4'     2        0      'In'        
%    4       'Pin5'     3        0      'In'        
%    5       'Pin6'     4        0      'In'        
%    6       'Pin7'     5        0      'In'        
%    7       'Pin8'     6        0      'In'        
%    8       'Pin9'     7        0      'In'   

% Check the properties of the digital input/output object
get(dio)

% Read the state of the parallel port
% http://www.mathworks.com/help/toolbox/daq/f11-20064.html

% Read value of the lines
linevalue_array = getvalue(dio);
last_value = binvec2dec(linevalue_array);
fprintf('\n\nInitial Value: %d', last_value);
datestr(now)

% Set value of the lines
% data = 23;
% putvalue(dio,data)

KEEP_GOING = true;

while KEEP_GOING
    
    % Read value of the lines to confirm change
    linevalue_array = getvalue(dio);
    current_value = binvec2dec(linevalue_array);
    
    if( current_value ~= last_value )    
        fprintf('\n\nNew Value: %d', current_value);
        datestr(now)
    end    
    
    last_value = current_value;
end


%% Testing for loop [2012/02/09]

FILE = fopen('00000-output.txt', 'w');

for OUTPUTSTREAM = [1, FILE]
    fprintf(OUTPUTSTREAM, '\nTesting output to %d', OUTPUTSTREAM);
end

%% Start HBD Training [2012/01/21]


% Set up the experiment specifications for the HBD program
experiment_specs.PPid = 'DEBUGGING';

experiment_specs.currentDirectory   = pwd;
experiment_specs.SoftwareVersion    = 'DEBUGGING-2012.01.21';
experiment_specs.SoftwareTitle      = 'Heartbeat Detection Training';

experiment_specs.HideMousePointer = false;
experiment_specs.MinBorderPercent = 5;
experiment_specs.FullScreenMode = false;


experiment_specs.filename_BioLab_Trials = '../data_03_heart_beat_detection_trials.txt';
experiment_specs.NumTrials_HBD = 10;
experiment_specs.SendTriggers = true;

experiment_specs.WindowPixels_Width = 600;
experiment_specs.WindowPixels_Height = 400;    

experiment_specs.SetSyncTest_StDev = true;
experiment_specs.VBL_MaxStd_ms = 5;

% Run the program
HBD_training_expt(experiment_specs);



%% Open files for editing



%% Test reading trigger structure [2011/11/16]

% Options
filename_trigger_list = 'data/trigger_list-Felt_Affect-2011.11.16.txt';

%--------------------------------------------------------------            
% Define event triggers for psychophysiology acquisition
% Uses Matab's Data Acquisition Toolbox and part of Spencer's Signals Approach Toolbox.
% See trigger_demo.m for more information. Email spencer.lynn@gmail.com
%
% Open a link to computer's parallel port for output of trigger/event codes.            
global parallel_port triggerinfo

clear parallel_port
clear triggerinfo

% Parallel port bit designations: pins 9 to 2 (of 25) are bits
% [7 6 5 4 3 2 1 0].
hlines = 0:7; 

% To create output lines on port.
direction = 'out'; 

%Create link to pararallel port for trigger output. The name of
%capturing variable must match Trigger_Port value in trigger
%list files.
parallel_port = init_dioport('parallel','LPT1',hlines,direction);

%Text file describes trigger names, codes, ports, durations.
%Used by fetch_signal, accrue_payoffs, etc.
triggerinfo = read_list( filename_trigger_list ); 

%{
% Create triggers
triggerNumber = 1;
% Human-readable information.
triggerinfo(triggerNumber).Event_Name       = {'Stimulus onset'};

% Port over which to send the trigger event. SAtb only works on
% parallel port at present (9/12/2011).
triggerinfo(triggerNumber).Trigger_Port     = {'parallel_port'}; 

% Decimal represenataion of the the trigger value to be
% received by remote hardware, sent over the port.
triggerinfo(triggerNumber).Trigger_Value    = 10;

% Duration of trigger event, in seconds. Not a stimulus
% duration, just the duration of the trigger itself. Smaller is
% better, to the limit of the acquisition hardware.
triggerinfo(triggerNumber).Trigger_Duration = 0.0100;

% Create triggers
% Human-readable information.
triggerNumber = triggerNumber + 1;
triggerinfo(triggerNumber).Event_Name       = {'Stimulus onset'};
triggerinfo(triggerNumber).Trigger_Port     = {'parallel_port'}; 
triggerinfo(triggerNumber).Trigger_Value    = 10;
triggerinfo(triggerNumber).Trigger_Duration = 0.0100;
%}

% For increased speed, dereferece the port ID strings in
% trigger, replace with the port uddobj/dioobj value.
for i=1:length(triggerinfo)
%                 fprintf('\n\nLOOP %d', i);
%                 triggerinfo(i)
%                 triggerinfo(i).Trigger_Port
%                 triggerinfo(i).Trigger_Port(1)

    triggerinfo(i).Trigger_Port=eval(triggerinfo(i).Trigger_Port{1});
end

% print to command-window
strcmpi( 'Baseline Start', triggerinfo(1).Event_Name )

% -- end trigger initialization parts --
%--------------------------------------------------------------
