% Ian Kleckner
% ian.kleckner@gmail.com
% Interdisciplinary Affective Science Lab (IASLab)
% Northeastern University
%
% Heartbeat Tracking Task
%
% Requirements
%  MindWare BioLab software for HBD
%
% 2012/05/22 Start coding
% 2017/01/01 Update to include Shimmer sensors
% 2017/05/18 Reduce rest delay from 30 sec to 1 sec
% 2017/05/26 Improve timing and instructions
%
% TODO

function HeartbeatTracking_expt_Shimmer( experiment_specs )
    %% Parse the input variable

    % For detecting R spikes
    Minimum_RR_Interval_sec = experiment_specs.Minimum_RR_Interval_sec;
    Minimum_R_Prominence_mV = experiment_specs.Minimum_R_Prominence_mV;
    
    % Read in COM port and format it as a string
    comPort_ECG = sprintf('%d', experiment_specs.COM_ECG);
    
    % Sample rate in [Hz]
    sampling_rate_ECG_Hz = experiment_specs.Sampling_rate_ECG_Hz;
    
    PPid                = sprintf('s%03d', experiment_specs.PP_Number);
    Timepoint           = experiment_specs.Timepoint;
        
    currentDirectory        = experiment_specs.currentDirectory;
    SoftwareVersion         = experiment_specs.SoftwareVersion;
    HideMousePointer        = experiment_specs.HideMousePointer;
    MinBorderPercent        = experiment_specs.MinBorderPercent;
    FullScreenMode          = experiment_specs.FullScreenMode;        
    
    %% OPTIONS    
    
    % Fontsize for asking how many HBs were counted
    % Fontsize for asking temperature when in fullscreen
    fontsize_questions = 75;
    
    if( ~FullScreenMode )
        fontsize_questions = 25;
    end
    
    % Order of the delay times listed below
    RANDOMIZE_TRIAL_ORDER = 0;
    
    if( ~experiment_specs.SpeedMode )
        % Duration for counting
        trial_duration_array    = [25, 35, 45];
        
        % Duration for rest between trials (if zero, the slide is skipped)
        duration_rest           = 0;     
        
        % Duration for get ready slide before counting starts
        delay_get_ready         = 1;
        
        % Duration for the stop sign before PP is asked how many HBs they
        % counted
        delay_stop_signal       = 3;
    else
        trial_duration_array    = [2.5, 3.5, 4.5];
        duration_rest           = 0;        
        delay_get_ready         = 0.5;
        delay_stop_signal       = 0.3;
    end
    
    delay_zero = 0;
    
    %----------------------------------------------------------------------
    % Inputs for detecting heartbeats

    % ECG channel used to detect R spikes
    string_ECG_signal_HBD = 'ECG LA-RA';
    %k_ECG_HBD = k_ECG_LA_RA;
    
    
    %% Program options    
    % Name of the program that goes on the output directory
    programName             = 'HBT';
    
    % Filenames of key slides
    filename_please_wait        = 'data/HeartbeatTracking/instructions/instructions-HBT-please_wait.jpg';
    filename_are_you_ready      = 'data/HeartbeatTracking/instructions/instructions-HBT-are_you_ready.jpg';
    filename_get_ready          = 'data/HeartbeatTracking/instructions/instructions-HBT-get_ready.jpg';        
    filename_start_counting     = 'data/HeartbeatTracking/instructions/instructions-HBT-start.jpg';
    filename_stop_counting      = 'data/HeartbeatTracking/instructions/instructions-HBT-stop.jpg';
    filename_task_is_over       = 'data/HeartbeatTracking/instructions/instructions-HBT-task_is_over.jpg';
    filename_how_many_HBs       = 'data/HeartbeatTracking/instructions/instructions-HBT-how_many_HBs.jpg';
    filename_rate_confidence    = 'data/HeartbeatTracking/instructions/instructions-HBT-how_confident.jpg';
    filename_saving_data        = 'data/HeartbeatTracking/instructions/saving_data.jpg';
    
    %% Ready user to start program
    userInput = questdlg('Click Start to connect to Shimmer and begin instruction slides', 'Ready to start?', ...
        'Start', 'Abort', 'Start');

    % Abort if the user wants to
    if( strcmp(userInput, 'Abort') )
        clear all;
        error('User aborted program instead of starting program');
    end

    % Start timing experiment
    timer_taskStart = tic();
     
    
    %% Initialize
    % Clear command window
    clc;
    
    % To address a bug
    cd(currentDirectory);    
    
    % Create an output directory to save results to files
    output_dir = sprintf('out-%s-%s-%s-T%d',datestr(now, 'yyyy.mm.dd-HHMM'), programName, PPid, Timepoint);
    mkdir(output_dir);
    
    % Create file for log
    outputfile = sprintf('%s/log.txt', output_dir);
    FILE = fopen(outputfile, 'a'); % 'a' for append
    
    % Write header
    writeSpecificationsHeader( FILE, experiment_specs )
    fprintf(FILE, '\n\n');
    
    % Backup code into this directory
    %copyCode( output_dir )        
    zip( sprintf('%s/code', output_dir),{'*.m','*.fig'});
    
    % Start saving contents of command window to log file
    diary( sprintf('%s/diary.txt', output_dir) )
    
    % Write header for program output
    OUTPUTSTREAM = 1; % 1 => command window
    writeSpecificationsHeader( OUTPUTSTREAM, experiment_specs );
    
    % Obtain new seed for random number generator [2011/11/15]
    %  This way, subsequent instances of MATLAB will be random with respect
    %  to one another
    %   On MATLAB R2011b: rng shuffle
    newStream = RandStream('mt19937ar', 'Seed', sum(100*clock))
    try
        RandStream.setDefaultStream(newStream);   
    catch exception
        RandStream.setGlobalStream(newStream);
    end
    
    fprintf('\n\n');
    fprintf('\n___________________________________________________________\n');
    
    % Enable unified mode of KbName, so KbName accepts identical key names on
    % all operating systems (KbDemo.m)
    KbName('UnifyKeyNames');

    % This script calls Psychtoolbox commands available only in OpenGL-based
    % versions of the Psychtoolbox.
    AssertOpenGL;
    
    Screen('CloseAll');
    
    
    %% Connect to Shimmer device

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
    
    %% Audio / Sound Setup    
    
    % The audio tone played by the computer
    tone_frequency_Hz_ready = 523.2511;
    tone_duration_sec_ready = 0.2;
    tone_sampling_rate_Hz_ready = 48000;

    %------------------------------------------------------------------
    % Open sound interface        
    InitializePsychSound(1);

    % Open audio device for low-latency output
    deviceid = [];
    audio_mode = 1; % Playback only
    reqlatencyclass = 0;
    nrchannels = 1;
    buffersize = 0;        
    repetitions = 1;

    % Hack to accomodate bad Windows systems or sound cards. By default,
    % the more aggressive default setting of something like 5 msecs can
    % cause sound artifacts on cheaper / less pro sound cards:
    suggestedLatencySecs = 0.1;
    warning('Using latency on Windows machine');

    pahandle_ready = PsychPortAudio('Open', ...
        deviceid, audio_mode, reqlatencyclass, ...
        tone_sampling_rate_Hz_ready, nrchannels, buffersize, suggestedLatencySecs);

    tone_sound_array = MakeBeep(tone_frequency_Hz_ready, tone_duration_sec_ready, tone_sampling_rate_Hz_ready);
    PsychPortAudio('FillBuffer', pahandle_ready, tone_sound_array);

    % Test and run first audio presentation to initialize the system         
    current_time_GetSecs = GetSecs();
    delay_until_play_sec = 0.1;
    time_beep_actual_sec = PsychPortAudio('Start', pahandle_ready, repetitions, current_time_GetSecs + delay_until_play_sec, 1);

    fprintf('\nTest presentation of delay %f msec was late by %f msec', 1e3*delay_until_play_sec, 1e3 * (time_beep_actual_sec - current_time_GetSecs - delay_until_play_sec))

    pause(0.5);
    
    %----------------------------------------------------------------------
    % For the TOUCH beep
    tone_frequency_Hz_touch = 1046.502;
    tone_duration_sec_touch = 0.3;
    tone_sampling_rate_Hz_touch = 48000;
    
    % Open audio device for low-latency output
    deviceid = [];
    audio_mode = 1; % Playback only
    reqlatencyclass = 0;
    nrchannels = 1;
    buffersize = 0;        
    repetitions = 1;

    % Hack to accomodate bad Windows systems or sound cards. By default,
    % the more aggressive default setting of something like 5 msecs can
    % cause sound artifacts on cheaper / less pro sound cards:
    suggestedLatencySecs = 0.1;
    warning('Using latency on Windows machine');

    pahandle_go = PsychPortAudio('Open', ...
        deviceid, audio_mode, reqlatencyclass, ...
        tone_sampling_rate_Hz_touch, nrchannels, buffersize, suggestedLatencySecs);

    tone_sound_array_dummy = MakeBeep(tone_frequency_Hz_touch, tone_duration_sec_touch, tone_sampling_rate_Hz_touch);
    PsychPortAudio('FillBuffer', pahandle_go, tone_sound_array_dummy);

    % Test and run first audio presentation to initialize the system         
    current_time_GetSecs = GetSecs();
    delay_until_play_sec = 0.5;
    time_beep_actual_sec = PsychPortAudio('Start', pahandle_go, repetitions, current_time_GetSecs + delay_until_play_sec, 1);


    %% Set up the display session
    %  See DisplaySession.m for more information

    % Create a DisplaySession
    DS = DisplaySession;

    % The minimum border percentage for scaling displayed images (<50)
    DS.setMinBorderPercent( MinBorderPercent );

    % Set the eye dominance
    %DS.setEyeDominance( DominantEyeString )

    if( FullScreenMode )
        rect_windowed = [];        
    else
        X_left = 0;
        Y_upper = 0;
        X_width = experiment_specs.WindowPixels_Width;
        Y_height = experiment_specs.WindowPixels_Height;        

        rect_windowed = [X_left, Y_upper, X_left+X_width, Y_upper+Y_height];
    end

    fprintf('\n\n_____________________________________');
    fprintf('\nSetting up display 1 [%s]', datestr(now));
    fprintf('\nAdding window rectangle as follows:');
    disp(rect_windowed);

    % Allow the screen refresh rate to vary with stdev +/- maxStddev
    %  This is to address an error in iMac screen setting
    if( experiment_specs.SetSyncTest_StDev )
        VBL_MaxStd_sec = experiment_specs.VBL_MaxStd_ms / 1000;

        Screen('Preference', 'SyncTestSettings', VBL_MaxStd_sec);

        fprintf('\n\n!! Sync test settings manually set to %f msec', VBL_MaxStd_sec * 1000);
        fprintf('\nThis determines time resolution of presentation');
    end

    % Stereo-mode as per the PTB documentation (see below)
    %  0 -> no stereo, 4 -> split left/right
    stereoMode = 0;

    % Create the PTB window for display
    [window, rect] = DS.createWindow(rect_windowed, stereoMode);

    % Hide the mouse cursor
    if( HideMousePointer )
        HideCursor();
    end
    
    %% Set up the display session
    %  See DisplaySession.m for more information
    
    % Create a DisplaySession
    DS = DisplaySession;

    % The minimum border percentage for scaling displayed images (<50)
    DS.setMinBorderPercent( MinBorderPercent );

    % Set the eye dominance
    %DS.setEyeDominance( DominantEyeString )

    if( FullScreenMode )
        rect_windowed = [];        
    else
        X_left = 0;
        Y_upper = 0;
        X_width = experiment_specs.WindowPixels_Width;
        Y_height = experiment_specs.WindowPixels_Height;        
        
        rect_windowed = [X_left, Y_upper, X_left+X_width, Y_upper+Y_height];
    end
    
    fprintf('\n\n_____________________________________');
    fprintf('\nSetting up display 1 [%s]', datestr(now));
    fprintf('\nAdding window rectangle as follows:');
    disp(rect_windowed);
    
    % Allow the screen refresh rate to vary with stdev +/- maxStddev
    %  This is to address an error in iMac screen setting
    if( experiment_specs.SetSyncTest_StDev )
        VBL_MaxStd_sec = experiment_specs.VBL_MaxStd_ms / 1000;
        
        Screen('Preference', 'SyncTestSettings', VBL_MaxStd_sec);
        
        fprintf('\n\n!! Sync test settings manually set to %f msec', VBL_MaxStd_sec * 1000);
        fprintf('\nThis determines time resolution of presentation');
    end

    % Stereo-mode as per the PTB documentation (see below)
    %  0 -> no stereo, 4 -> split left/right
    stereoMode = 0;
    
    % Create the PTB window for display
    [window, rect] = DS.createWindow(rect_windowed, stereoMode);    
    
    % Hide the mouse cursor
    if( HideMousePointer )
        HideCursor();
    end
    
    %% Using input
    
    %----------------------------------------------------------------------
    % Allow input of buttons 1, 2, 3, or 4
    if( ismac )
        % You pressed key 30 which is 1!
        % You pressed key 31 which is 2@
        % You pressed key 32 which is 3#
        % You pressed key 33 which is 4$
        % You pressed key 41 which is ESCAPE
        
        responseKey_commit = 44;
        responseKey_abort = 41;
        %error('Program keys for MAC');

    elseif( IsWin )  
        % You pressed key 49 which is 1!
        % You pressed key 50 which is 2@
        % You pressed key 51 which is 3#
        % You pressed key 52 which is 4$
        % You pressed key 27 which is ESCAPE
        
        responseKey_commit = 32;
        responseKey_abort = 27;
        %error('Program keys for WINDOWS');

    elseif( IsLinux )
        % You pressed key 11 which is 1!
        % You pressed key 12 which is 2@
        % You pressed key 13 which is 3#
        % You pressed key 14 which is 4$
        % You pressed key 10 which is ESCAPE
        
        responseKey_commit = 66;
        responseKey_abort = 10;
        %error('Program keys for LINUX');

    else
        error('Invalid OS');
    end
    
    %----------------------------------------------------------------------      
    % Prepare a slider scale for rating confidence
    sliderScale_Confidence = SliderScale(DS);
    sliderScale_Confidence.setDisplayImage(true, filename_rate_confidence);            
    sliderScale_Confidence.setDisplayText(false);
    sliderScale_Confidence.setResponseKeys(responseKey_commit, responseKey_abort);
    
    %% Instructions
    
    DS.addStep('NAVIGATE', 'STEREO', 'data/HeartbeatTracking/instructions/instructions-HBT-01.jpg', []);
    DS.addStep('NAVIGATE', 'STEREO', 'data/HeartbeatTracking/instructions/instructions-HBT-02.jpg', []);
    DS.addStep('NAVIGATE', 'STEREO', 'data/HeartbeatTracking/instructions/instructions-HBT-03.jpg', []);
    DS.addStep('NAVIGATE', 'STEREO', 'data/HeartbeatTracking/instructions/instructions-HBT-04.jpg', []);
    
    %DS.addStep('NAVIGATE', 'STEREO', 'data/HeartbeatTracking/instructions/instructions-HBT-05.jpg', []);
    %DS.addStep('NAVIGATE', 'STEREO', 'data/HeartbeatTracking/instructions/instructions-HBT-06.jpg', []);
    %DS.addNewSliderScale( sliderScale_Confidence );
    
    DS.addStep('NAVIGATE', 'STEREO', 'data/HeartbeatTracking/instructions/instructions-HBT-07.jpg', []);
    DS.addStep('NAVIGATE', 'STEREO', 'data/HeartbeatTracking/instructions/instructions-HBT-08.jpg', []);
    DS.addStep('NAVIGATE', 'STEREO', 'data/HeartbeatTracking/instructions/instructions-HBT-09.jpg', []);
    DS.addStep('NAVIGATE', 'STEREO', 'data/HeartbeatTracking/instructions/instructions-HBT-10.jpg', []);
    DS.addStep('NAVIGATE', 'STEREO', 'data/HeartbeatTracking/instructions/instructions-HBT-11.jpg', []);
    
    % Start displaying instructions
    DS.displayStep(1);
    
    %% Do the task
    
    Ntrials = length(trial_duration_array);
    
    if( RANDOMIZE_TRIAL_ORDER )
        k_random = randperm(Ntrials);
        trial_duration_array = trial_duration_array(k_random);
    end
    
    % Press a key to get ready for trial
    %DS.addStep('NAVIGATE', 'STEREO', filename_please_wait, []);
    
    % Header for log file
    fprintf(FILE, '\nTrialNum\tTrialDuration(sec)\tNum_HBs_Counted(reported)\tNum_HBs_Counted(convert_to_integer)\tNum_Rspikes_Detected\tHBT_Score\tConfidence(0_to_1)\tConfidenceRT(sec)');

    HBT_score_array = zeros(1,Ntrials);
    
    for t = 1:Ntrials
        trial_duration = trial_duration_array(t);
        
        fprintf('\nTrial %d / %d', t, Ntrials);
        fprintf('\n\tDelay\t%f', trial_duration);   
        
        % If this is not the last trial, use a rest period
        if( duration_rest > 0 )
            DS.addStep('DELAY', 'STEREO', filename_please_wait, [], delay_zero);
            DS.nextStepDisplay();
            WaitSecs(duration_rest); 
        end
                
        % Are you ready?
        DS.addStep('NAVIGATE', 'STEREO', filename_are_you_ready, []);
        
        % Get ready
        DS.addStep('DELAY', 'STEREO', filename_get_ready, [], delay_get_ready);        
        DS.nextStepDisplay();
        
        %-----------------------------
        % Start trial      
        
        
        % Show the start cue
        DS.addStep('DELAY', 'STEREO', filename_start_counting, [], delay_zero);
        
        % Get ready...
        current_time_GetSecs = GetSecs();
        delay_until_play_sec = 1.5;
        time_beep_actual_sec = PsychPortAudio('Start', pahandle_ready, repetitions, current_time_GetSecs + delay_until_play_sec, 1);
        PsychPortAudio('Stop', pahandle_ready);
        fprintf('\nReady...');

        % ...set...
        delay_until_play_sec = 2.5;
        time_beep_actual_sec = PsychPortAudio('Start', pahandle_ready, repetitions, current_time_GetSecs + delay_until_play_sec, 1);
        PsychPortAudio('Stop', pahandle_ready);
        fprintf('set...');

        % ...GO!
        delay_until_play_sec = 3.5;
        time_beep_actual_sec = PsychPortAudio('Start', pahandle_go, repetitions, current_time_GetSecs + delay_until_play_sec, 1);
        fprintf('GO');
        time_touch_on_GetSecs = GetSecs();
        time_touch_off_GetSecs = time_touch_on_GetSecs + 0.5;
        
        % Start displaying the steps
        DS.nextStepDisplay();
        
        % Start getting data from Shimmer
        shimmer.start;

        % Wait the duration of the slide
        WaitSecs(trial_duration);
        
        %------------------------------------------------------------------
        % Show the end cue
        DS.addStep('DELAY', 'STEREO', filename_stop_counting, [], delay_stop_signal);        
        %DS.addStep('NAVIGATE', 'STEREO', filename_how_many_HBs, []);
        
        % ...STOP!
        delay_until_play_sec = 0;
        time_beep_actual_sec = PsychPortAudio('Start', pahandle_go, repetitions, current_time_GetSecs + delay_until_play_sec, 1);
        fprintf('STOP');
        time_touch_on_GetSecs = GetSecs();
        time_touch_off_GetSecs = time_touch_on_GetSecs + 0.5;
        
        DS.nextStepDisplay();
        
        % Ask how many heartbeats were counted
        text_string = 'How many heartbeats\ndid you count?\n\nType number, hit Enter:';
        text_y = 0.2 * rect(4);

        question_string = '';
        question_x = (rect(3) - rect(1)) / 2.5;
        question_y = 0.7 * rect(4);
        
        % Draw the text
        oldTextSize = Screen('TextSize', window, fontsize_questions);
        DrawFormattedText(window, text_string, 'center', text_y, 255*[1 1 1]);
        [NHBs_counted_string, ~] = GetEchoString(window, question_string, question_x, question_y, 255*[1 1 1], [0 0 0] );

        % Ask for confidence
        DS.addNewSliderScale( sliderScale_Confidence );
        DS.nextStepDisplay();
        
        %------------------------------------------------------------------
        % Process data for analysis
        
        DS.addStep('DELAY', 'STEREO', filename_saving_data, [], 0.1);
        DS.nextStepDisplay();
        
        % Read the latest data from shimmer data buffer, signalFormatArray defines the format of the data and signalUnitArray the unit
        [newData, signalNameArray, signalFormatArray, signalUnitArray] = shimmer.getdata('c');
        
        % Stop recording data
        shimmer.stop;
        
        % Process data
        k_ECG_Signal_HBD = find(strcmp(string_ECG_signal_HBD, signalNameArray));
        k_time = find(strcmp('Time Stamp', signalNameArray));
        
        X = newData(:, k_time) / 1e3;
        Y = newData(:, k_ECG_Signal_HBD);
        sampling_period = X(2) - X(1);
        
        % Cut this to duration of trial because it will be too long
        Npoints_trial = trial_duration / sampling_period;
        
        X = X(1:Npoints_trial);
        Y = Y(1:Npoints_trial);
        Yraw = Y;
        
        for ADJUSTED_ECG = [0 1]
            
            if( ADJUSTED_ECG )                
                Y = ECG_adjust_baseline_spline( X, Yraw, Minimum_RR_Interval_sec, Minimum_R_Prominence_mV );
            end
            
            % Peak detection
            [peak_Y_array, peak_X_array, peak_width_array, peak_prom_array] = ...
                findpeaks(Y, X, 'MinPeakDistance', Minimum_RR_Interval_sec);

            % Narrow to the valid peaks based on
            % prominence
            k_valid_peaks = peak_prom_array > Minimum_R_Prominence_mV;
            peak_Y_array = peak_Y_array(k_valid_peaks);
            peak_X_array = peak_X_array(k_valid_peaks);
            peak_prom_array = peak_prom_array(k_valid_peaks);
            NRspikess_detected = length(peak_X_array);

            % Plot data and detected peaks
            hfig = figure('Visible', 'off');
            plot(X - X(1), Y, '-k');

            hold('all');
            plot(peak_X_array - X(1), peak_Y_array, 'or', 'LineWidth', 3);

            title(sprintf('%s, Trial %d, %0.1f sec, %d detected HBs Adjusted=%d', ...
                PPid, t, trial_duration, NRspikess_detected, ADJUSTED_ECG));
            xlabel('Time (sec)');
            ylabel('ECG mv');
            print(hfig, sprintf('%s/out-ECG-Trial_%02d-Adj%d.png', output_dir, t, ADJUSTED_ECG), '-dpng');
            close(hfig);

            %------------------------------------------------------------------
            % Save ECG data to file
            filename_out = sprintf('%s/out-ECG-Trial_%02d-Adj%d.csv', output_dir, t, ADJUSTED_ECG);

            % Write header
            FILE_OUT = fopen(filename_out, 'w');        
            fprintf(FILE_OUT, '%s,%s\n', 'Time_sec', 'ECG_mV');
            fclose(FILE_OUT);

            % Write data
            output_data = [X-X(1) Y];
            dlmwrite(filename_out, output_data, '-append', 'delimiter', ',');            
        end
        
        %------------------------------------------------------------------
        numQuestions = 1;
        
        % Get input from PP for the most recent questions
        for recencyIndex = 1:numQuestions
            % Get response from question number 1, 2, 3, 4, etc.
            [response_value_array(numQuestions-recencyIndex+1), ...
             response_time_array(numQuestions-recencyIndex+1)] = ...
             DS.getRecentResponse( recencyIndex );
        end
        
        confidence_array(t) = response_value_array(1);
        confidence_RT_array(t) = response_time_array(1);
        
        %------------------------------------------------------------------
        % Calculate score for this trial and save to array for later
        NHBs_counted = sscanf(NHBs_counted_string, '%d');
        
        if( isempty(NHBs_counted) )
            NHBs_counted = 0;
        end
        HBT_score_array(t) = 1 - (NRspikess_detected - NHBs_counted) / NRspikess_detected;
        
        %------------------------------------------------------------------        
        % Write to log file
        fprintf('\n\tConfidence\t%f\tin\t%f sec', confidence_array(t), confidence_RT_array(t));        
        fprintf(FILE, '\n%d\t%f\t%s\t%d\t%d\t%f\t%f\t%f', t, trial_duration, NHBs_counted_string, NHBs_counted, NRspikess_detected, HBT_score_array(t), confidence_array(t), confidence_RT_array(t));
        
        % User wants to exit
        if( any(response_value_array == -1) )
            fprintf('\n\nUser exits');            
            break;
        end
    end
    
    % Write average score to file
    fprintf(FILE, '\n\nTotal HBT score\t%f', mean(HBT_score_array));
    
    %% It's over!!
    shimmer.disconnect;
    
    % The task is over
    DS.addStep('NAVIGATE', 'STEREO', filename_task_is_over, filename_task_is_over);
    DS.nextStepDisplay();   
    
    % Close the log file
    fclose(FILE);
    
    % Close the onscreen window
    Priority(0);
    Screen('CloseAll')
    
    fprintf('\nAll done!\n');
    diary('off');    
end