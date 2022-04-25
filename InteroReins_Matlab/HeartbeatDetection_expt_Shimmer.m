% Ian Kleckner
% Univeristy of Rochester Medical Center
% ian.kleckner@gmail.com
%
% Testing heartbeat detection task using Shimmer ECG hardware
% Some code copied from Shimmer plotandwriteecgexample.m ( (C) Shimmer )
%
% Requirements
%  PsychToolbox
%
% 2016/12/22 Start coding
% 2017/05/26 Improve instructions and user input
% 2019/04/15 Add trial number to view
%            Dummy first HB to prevent timing error
%            Save image of each trial for post quality assessment
%            Improve timing of beep delivery to
%              prevent beeps from going too late (>200 or >500 msec after
%              spike detection)
%            Allow researcher to repeat a trial mid way if there is an error
%            Add real-time plots of progress
%            Can run task on second monitor (subject views second monitor),
%              researcher views primary monitor

function HeartbeatDetection_expt_Shimmer( experiment_specs )

    %% Parse the input variable
    
    % Number of trials
    Ntrials_HBD = experiment_specs.Number_of_trials_HBD;
    
    % For detecting R spikes
    Minimum_RR_Interval_sec = experiment_specs.Minimum_RR_Interval_sec;
    Minimum_R_Prominence_mV = experiment_specs.Minimum_R_Prominence_mV;
    
    % Read in COM port and format it as a string
    comPort_ECG = sprintf('%d', experiment_specs.COM_ECG);
    
    % Sample rate in [Hz]
    sampling_rate_ECG_Hz = experiment_specs.Sampling_rate_ECG_HBD_Hz;

    PPid                = sprintf('s%03d', experiment_specs.PP_Number);
    Timepoint           = experiment_specs.Timepoint;

    currentDirectory        = experiment_specs.currentDirectory;
    SoftwareVersion         = experiment_specs.SoftwareVersion;    

    HideMousePointer        = experiment_specs.HideMousePointer;
    MinBorderPercent        = experiment_specs.MinBorderPercent;
    FullScreenMode          = experiment_specs.FullScreenMode;    

    %% Other input    
    
    % Use slider scale or use numerical scale?
    RATE_CONFIDENCE_W_SLIDER_SCALE  = false;
    KbName('UnifyKeyNames'); 
    key_conf_guess                  = KbName('`~');
    key_conf_low                    = KbName('1!');
    key_conf_high                   = KbName('2@');
    key_yes                         = KbName('1!');
    key_no                          = KbName('2@');
    key_shift                       = KbName('LeftShift');
    key_esc                         = KbName('ESCAPE');
    KbDevice_Number                 = -1;
    
    % Use Psychtoolox window for instructions, etc.
    USE_PTB_WINDOW = false;
    
    % 2019/06/14 First practice trail is synchronous
    USE_PRACTICE_TRIAL_SYNCHRONOUS = true;
    
    % Show a realtime data plot for each trial with ECG and R spikes
    SHOW_REALTIME_PLOT_PER_TRIAL = true;
    
    % Helps screen out bad trials afterwards
    WRITE_FIGURE_PER_TRIAL = true;
       
    % Fontsize for writing trial number when in fullscreen
    fontsize_questions = 75;
    
    if( ~FullScreenMode )
        fontsize_questions = 25;
    end
    
    % Extra debugging output for detecting R spikes
    OUTPUT_DEBUG = false;
    
    %----------------------------------------------------------------------
    % HBD Task Options
    
    % Don't wait longer than this time between R spikes. If it's been this
    % long then there is something wrong and the trial will be repeated.
    %  DEFAULT = 3 just so there has to be a major error before aborting
    %  the trial
    Maximum_RR_Interval_sec = 3;
    
    % Delay times for coincident and non-coincident trials
    delay_HB_coincident_sec = 0.2;
    delay_HB_noncoincident_sec = 0.5;
    
    % If time to beep is THIS long then there is clearly a problem
    MAX_TIME_TO_BEEP_SEC = 0.7;

    % Number of tones per trial
    Nbeeps_per_HBD_trial = 10;
    
    % The audio tone played by the computer
    tone_frequency_Hz = 523.2511;
    tone_duration_sec = 0.1;
    tone_sampling_rate_Hz = 48000;

    % Minimum delay between showing HBD result (when providng feedback)
    % and proceeding to next trial
    delay_HB_result_sec     = 2; 

    % Delay for when the PP sees they are 25% 50% or 75% done
    delay_break_sec         = 1;

    filename_progress_25_percent = 'data/HeartbeatDetection/instructions/instructions-25_percent_done.jpg';
    filename_progress_50_percent = 'data/HeartbeatDetection/instructions/instructions-50_percent_done.jpg';
    filename_progress_75_percent = 'data/HeartbeatDetection/instructions/instructions-75_percent_done.jpg';

    %----------------------------------------------------------------------
    % Inputs for detecting heartbeats

    % ECG channel used to detect R spikes
    string_ECG_signal_HBD = 'ECG LA-RA';
    %k_ECG_HBD = k_ECG_LA_RA;
    
    % This is how much time between checking the Shimmer for data (smaller
    % numbers increase temporal precision but risk an error if there are no
    % data available)
    sampling_period_check_Shimmer_sec = 0.050;

    %% Program options
    
    % Training mode or Assessment mode
    TrainingMode            = experiment_specs.TrainingMode;

    if( TrainingMode )
        programName = 'HBD-Train';
    else
        programName = 'HBD';
    end

    % Key insruction slides
    %filename_title                      = 'data/HeartbeatDetection/instructions/instructions-title.jpg';    
    %filename_waiting                    = 'data/HeartbeatDetection/instructions/instructions-waiting_to_start.jpg';
    %filename_waiting                    = 'data/HeartbeatDetection/instructions/instructions-proceed.jpg';
    filename_trial_underway             = 'data/HeartbeatDetection/instructions/instructions-trial_underway.jpg';
    filename_request_input              = 'data/HeartbeatDetection/instructions/instructions-request_input.jpg';
    filename_correct_coincident         = 'data/HeartbeatDetection/instructions/instructions-correct-coincident.jpg';
    filename_correct_not_coincident     = 'data/HeartbeatDetection/instructions/instructions-correct-not_coincident.jpg';
    filename_incorrect_coincident       = 'data/HeartbeatDetection/instructions/instructions-incorrect-coincident.jpg';
    filename_incorrect_not_coincident   = 'data/HeartbeatDetection/instructions/instructions-incorrect-not_coincident.jpg';
    filename_proceed_to_trial           = 'data/HeartbeatDetection/instructions/instructions-proceed.jpg';
    filename_task_is_over               = 'data/HeartbeatDetection/instructions/instructions-task_is_over.jpg';
    filename_rate_confidence            = 'data/HeartbeatDetection/instructions/instructions-rate_confidence.jpg';
    
    filename_researcher_retry           = 'data/HeartbeatDetection/instructions/instructions-researcher_retry.jpg';
    filename_Rspike_or_beep_timeout     = 'data/HeartbeatDetection/instructions/instructions-beep_timeout.jpg';

    % Trial Criterion Mode - TONES or DURATION
    %  NUMBER_OF_TONES - when each trial has a specific number of tones (trial
    %  duration may
    %  DURATION - when each trial has a specific duration (number of tones
    %  may vary)
    trialCriteronString = 'NUMBER_OF_TONES';
    %trialCriteronString = 'DURATION';

    % A very short duration (for showing some slides)    
    duration_short = 0.001;
    
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
    %warning('Need to write specifications header once I start using Dynamic Table thing again');
    fprintf(FILE, '\n\n');

    % Backup code into this directory
    %copyCode( output_dir )        
    zip( sprintf('%s/code', output_dir),{'*.m','*.fig'});

    % Start saving contents of command window to log file
    diary( sprintf('%s/diary.txt', output_dir) )

    % Write header for program output
    OUTPUTSTREAM = 1; % 1 => command window
    %writeSpecificationsHeader( OUTPUTSTREAM, experiment_specs );
    warning('Need to write specifications header once I start using Dynamic Table thing again');

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


    %% Set up the display session
    %  See DisplaySession.m for more information

    if( USE_PTB_WINDOW )

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
    end

    %% Audio / Sound Setup    

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

    pahandle = PsychPortAudio('Open', ...
        deviceid, audio_mode, reqlatencyclass, ...
        tone_sampling_rate_Hz, nrchannels, buffersize, suggestedLatencySecs);

    tone_sound_array = MakeBeep(tone_frequency_Hz, tone_duration_sec, tone_sampling_rate_Hz);
    PsychPortAudio('FillBuffer', pahandle, tone_sound_array);

    % Test and run first audio presentation to initialize the system         
    current_time_GetSecs = GetSecs();
    delay_until_play_sec = 0.2;
    time_beep_actual_sec = PsychPortAudio('Start', pahandle, repetitions, current_time_GetSecs + delay_until_play_sec, 1);

    fprintf('\nTest presentation of delay %f msec was late by %f msec', 1e3*delay_until_play_sec, 1e3 * (time_beep_actual_sec - current_time_GetSecs - delay_until_play_sec))

    
    %----------------------------------------------------------------------
    % For the dummy beep 2019/04/23
    tone_frequency_Hz_dummy = 5;
    tone_duration_sec_dummy = 0.1;
    tone_sampling_rate_Hz_dummy = 48000;
    
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

    pahandle_dummy = PsychPortAudio('Open', ...
        deviceid, audio_mode, reqlatencyclass, ...
        tone_sampling_rate_Hz_dummy, nrchannels, buffersize, suggestedLatencySecs);

    tone_sound_array_dummy = MakeBeep(tone_frequency_Hz_dummy, tone_duration_sec_dummy, tone_sampling_rate_Hz_dummy);
    PsychPortAudio('FillBuffer', pahandle_dummy, tone_sound_array_dummy);

    % Test and run first audio presentation to initialize the system         
    current_time_GetSecs = GetSecs();
    delay_until_play_sec = 0.2;
    time_beep_actual_sec = PsychPortAudio('Start', pahandle_dummy, repetitions, current_time_GetSecs + delay_until_play_sec, 1);



    %% User inputs

    %----------------------------------------------------------------------
    % Allow input of buttons 1, 2, 3, or 4
    if( ismac )
        % You pressed key 30 which is 1!
        % You pressed key 31 which is 2@
        % You pressed key 32 which is 3#
        % You pressed key 33 which is 4$
        % You pressed key 41 which is ESCAPE
        responseKeys = {30 31 41};

        responseKey_commit = 44;
        responseKey_abort = 41;
        %error('Program keys for MAC');

    elseif( IsWin )  
        % You pressed key 49 which is 1!
        % You pressed key 50 which is 2@
        % You pressed key 51 which is 3#
        % You pressed key 52 which is 4$
        % You pressed key 27 which is ESCAPE
        % You pressed key 192 which is `~
        responseKeys = {49 50 27};

        responseKey_commit = 32;
        responseKey_abort = 27;
        %error('Program keys for WINDOWS');

    elseif( IsLinux )
        % You pressed key 11 which is 1!
        % You pressed key 12 which is 2@
        % You pressed key 13 which is 3#
        % You pressed key 14 which is 4$
        % You pressed key 10 which is ESCAPE
        responseKeys = {11 12 10};

        responseKey_commit = 66;
        responseKey_abort = 10;
        %error('Program keys for LINUX');

    else
        error('Invalid OS');
    end

    responseValue_coincident = 1;
    responseValue_not_coincident = 2;

    if( USE_PTB_WINDOW )
        % Key1 = 1, Key2=2, KeyEsc=-1, Key`=0 (for confidence)
        responseValues = [responseValue_coincident responseValue_not_coincident -1];  
        % Commit these keys
        DS.setInputKeys( responseKeys, responseValues );

        % Left on, middle on, right on
        DS.setInputKeys_Mouse( [responseValue_coincident, NaN, responseValue_not_coincident] )
        %----------------------------------------------------------------------  

        %----------------------------------------------------------------------      
        % Prepare a slider scale for rating confidence
        sliderScale_Confidence = SliderScale(DS);
        sliderScale_Confidence.setDisplayImage(true, filename_rate_confidence);            
        sliderScale_Confidence.setDisplayText(false);
        sliderScale_Confidence.setResponseKeys(responseKey_commit, responseKey_abort);
    end

    %% Instructions for task
    if( USE_PTB_WINDOW )
        if( ~experiment_specs.SpeedMode )
            DS.addStep('NAVIGATE', 'STEREO', 'data/HeartbeatDetection/instructions/instructions-HBD-01.jpg', []);
            DS.addStep('NAVIGATE', 'STEREO', 'data/HeartbeatDetection/instructions/instructions-HBD-02.jpg', []);
            DS.addStep('NAVIGATE', 'STEREO', 'data/HeartbeatDetection/instructions/instructions-HBD-03.jpg', []);
            DS.addStep('NAVIGATE', 'STEREO', 'data/HeartbeatDetection/instructions/instructions-HBD-04.jpg', []);
            DS.addStep('NAVIGATE', 'STEREO', 'data/HeartbeatDetection/instructions/instructions-HBD-05.jpg', []);    
            DS.addStep('NAVIGATE', 'STEREO', 'data/HeartbeatDetection/instructions/instructions-HBD-06.jpg', []);

            %DS.addStep('NAVIGATE', 'STEREO', 'data/HeartbeatDetection/instructions/instructions-HBD-07.jpg', []);
            %DS.addNewSliderScale( sliderScale_Confidence );

            DS.addStep('NAVIGATE', 'STEREO', 'data/HeartbeatDetection/instructions/instructions-HBD-08.jpg', []);    
            DS.addStep('NAVIGATE', 'STEREO', 'data/HeartbeatDetection/instructions/instructions-HBD-09.jpg', []);    
            DS.addStep('NAVIGATE', 'STEREO', 'data/HeartbeatDetection/instructions/instructions-HBD-10.jpg', []);    
            DS.addStep('NAVIGATE', 'STEREO', 'data/HeartbeatDetection/instructions/instructions-HBD-11.jpg', []);

            if( USE_PRACTICE_TRIAL_SYNCHRONOUS )
                DS.addStep('NAVIGATE', 'STEREO', 'data/HeartbeatDetection/instructions/instructions-HBD-12.jpg', []);
            end
        end        
        DS.displayStep(1);
    else
        fprintf('\n\n***************************************************************');
        fprintf('\nINSTRUCTIONS: The goal of this task is to assess how well you can detect your heartbeat');
        fprintf('\nINSTRUCTIONS: You will hear a series of 10 beeps, each of which is triggered by your heartbeat');
        fprintf('\nINSTRUCTIONS: Sometimes, the 10 beeps occur DURING your heartbeats');
        fprintf('\nINSTRUCTIONS: Other times, the 10 beeps occur IN BETWEEN your heartbeats');
        fprintf('\nINSTRUCTIONS: For each series of 10 beeps, indicate whether you heart the beeps DURING or BETWEEN your heartbeats');
        fprintf('\nINSTRUCTIONS: You will also be asked how CONFIDENT you are in your response');
        fprintf('\nINSTRUCTIONS: Please do not take your pulse directly using your hands or other objects');
        fprintf('\nINSTRUCTIONS: Instead, focus on your heart and chest');
        fprintf('\nINSTRUCTIONS: Finally, please remain still during the beeping part of each trial');
        fprintf('\nINSTRUCTIONS: You will be able to do a practice trial soon with beeps DURING each heartbeat');
        fprintf('\nINSTRUCTIONS: Do you have any questions?');
    end
    
    %% Start the HBD task

    % None of the progress has been shown
    progress_25_percent_was_shown = false;
    progress_50_percent_was_shown = false;
    progress_75_percent_was_shown = false;

    %----------------------------------------------------------------------
    % Set up trial structure
    trialIsCoincident_array = zeros(1,Ntrials_HBD);
    trialIsCoincident_array(1:round(Ntrials_HBD/2)) = 1;
    trialIsCoincident_array = Shuffle(trialIsCoincident_array);
    trialIsCoincident_array = trialIsCoincident_array == 1;

    delay_HB_sec_array = delay_HB_noncoincident_sec * ones(1, Ntrials_HBD);
    delay_HB_sec_array(trialIsCoincident_array) = delay_HB_coincident_sec;

    %----------------------------------------------------------------------
    % Initialize and start trials

    % Array of user input values as to whether HB is coincident with tone    
    userInput_value_array       = [];
    userInput_RT_array          = [];
    userInput_Coincident_array  = [];
    userIsCorrect_array         = [];

    % Signal detection results
    %
    %  __Stimulus__            __PP_Response__         __SDT_Result__
    %  Coicident                Coicident               Hit
    %  Coicident                Not-Coicident           Miss
    %  Not-Coincident           Coincident              False Alarm
    %  Not-Coincident           Not-Coincident          Correct Rejection
    % 
    Nhits                   = 0;
    Nmisses                 = 0;
    Nfalse_alarms           = 0;
    Ncorrect_rejections     = 0;
    
    REPEAT_THIS_TRIAL = false;

    % Header for output file
    fprintf(FILE, '\n\nTrialNum\tAttempt\tTrialDelay(sec)\tTrialIsCoincident\tUserResponse\tResponseTime(sec)\tUserIsCorrect\tConfidence(0_to_1)\tConfidenceRT(sec)');

    % Needs to be a WHILE loop not a FOR loop because of the way that I
    % programmed the ability to repeat trials
    t = 1;
    DONE_WITH_PRACTICE = false;
    while t <= Ntrials_HBD
        if( ~REPEAT_THIS_TRIAL )
            % This is the first time this trial is being attempted
            times_attempting_this_trial = 1;
        end
        
        % Hopefully don't need to repeat the trial...
        REPEAT_THIS_TRIAL = false;
        
        % Proceed to next trial - input from PP
        if( USE_PTB_WINDOW )
            DS.addStep('NAVIGATE_FORWARD_ANYKEY', 'STEREO', filename_proceed_to_trial, filename_proceed_to_trial);
            DS.nextStepDisplay();
        else
            fprintf('\n\n** Press any key to start the next trial');
            
            % Release all keys
            while KbCheck; end
        
            % Press any key
            keyIsDown = false;
            while( ~keyIsDown )
                % Check the state of the keyboard.
                [ keyIsDown, seconds, keyCode ] = KbCheck(KbDevice_Number);
            end
        end
        
        % Reset sampling bufer for first heartbeat
        sampling_buffer_first_HB_sec = 0.1;
        
        % Have not yet delivered the dummy beep (not a beep). This is to
        % help timing of actual beeps because it seems like the first beep
        % is often incorrect [2019/04/15]
        HAVENT_DONE_DUMMY_BEEP = true;
        
        if( USE_PRACTICE_TRIAL_SYNCHRONOUS )       
            delay_HB_msec       = 1000 * delay_HB_coincident_sec;
            delay_HB_sec        = delay_HB_coincident_sec;
            trialIsCoincident   = true;
            
            fprintf('\n\n---------------------------------------------------');
            fprintf('\n\n**** THIS IS A PRACTICE TRIAL ***');
            fprintf('\nTrial is COINCIDENT');
            fprintf('\n\n** PRESS SHIFT + ESC TO END THE BEEPING');
            
        else        
            % Obtain information about trial
            delay_HB_msec       = 1000 * delay_HB_sec_array(t);
            delay_HB_sec        = delay_HB_sec_array(t);
            trialIsCoincident   = trialIsCoincident_array(t);
            
            % Output to console
            fprintf('\n\n---------------------------------------------------');
            fprintf('\nTrial %d / %d, Attempt %d', t, Ntrials_HBD, times_attempting_this_trial);        
            fprintf('\nTrial HBD delay: %0.0f msec', delay_HB_msec);

            if( trialIsCoincident )
                fprintf('\nTrial is COINCIDENT');
            else
                fprintf('\nTrial is NOT COINCIDENT');
            end
        end
        
        % Wait until all keyboard keys are released
        fprintf('\nWaiting for keyboard keys to be released...');
        while KbCheck; end
        fprintf('done');
        
        if( SHOW_REALTIME_PLOT_PER_TRIAL && ~exist('hfig') )
            hfig = figure;
            ax = hfig.CurrentAxes;
        end

        % ID the trial as either a certain number of tones of a certain
        % duration
        switch trialCriteronString
            case 'NUMBER_OF_TONES'
                %DS.addStep('DELAY', 'STEREO', filename_trial_underway, filename_trial_underway, duration_short);
                %DS.nextStepDisplay();
                
                % Clear screen to black for the trial
                if( USE_PTB_WINDOW )
                    Screen('FillRect', window, [0 0 0], rect);
                end
                
                if( USE_PRACTICE_TRIAL_SYNCHRONOUS )
                    if( USE_PTB_WINDOW )
                        oldTextSize = Screen('TextSize', window, round(fontsize_questions));                    
                        trial_num_string = sprintf('Practice (coincident beeps)');
                        [newX,newY,textHeight]=Screen('DrawText', window, trial_num_string, 0.1*rect(4), 0.3*rect(4), 192*[ 1 1 0 ] );

                        trial_num_string = sprintf('Hold Shift + Esc to stop');
                        [newX,newY,textHeight]=Screen('DrawText', window, trial_num_string, 0.1*rect(4), 0.6*rect(4), 192*[ 1 1 0 ] );                    
                    end
                    
                else
                    if( USE_PTB_WINDOW )
                        % Draw the text with trial number
                        trial_num_string = sprintf('Trial %d    ', t);
                        oldTextSize = Screen('TextSize', window, round(fontsize_questions/2));
                        %DrawFormattedText(window, trial_num_string, 'right', 0.99*rect(4), 128*[1 1 1]);
                        [newX,newY,textHeight]=Screen('DrawText', window, trial_num_string, 0, 0, 128*[ 1 1 1 ] );

                        trial_num_string = sprintf('Researcher: Hold Shift + Esc to abort and re-try trial', t);
                        [newX,newY,textHeight]=Screen('DrawText', window, trial_num_string, 0.25*rect(4), 0.9*rect(4), 128*[ 1 1 1 ] );
                    end
                end
                
                if( USE_PTB_WINDOW )
                    Screen('Flip', window);
                end
                
                fprintf('\nNumber of tones: %d', Nbeeps_per_HBD_trial);

                % Initialize timers and other variables
                Npeaks_counted = 0;
                Npeaks_found_that_do_not_count = 0;
                Nbeeps_delivered = 0;                
                ALREADY_BEEPED_FOR_THIS_RSPIKE = true;
                FIRST_READ_HAS_OCCURRED = false;
                FOUND_PEAKS_ON_THIS_TRIAL = false;

                X_time_total_sec = [];
                Y_ECG_total = [];
                X_time_beep_sec_array = [];
                X_time_beep_late_sec_array = [];
                X_time_detected_peak_sec_array = [];
                Y_ECG_detected_peak_mV_array = [];
                Y_peak_prominance_array = [];

                % Start getting data from Shimmer
                shimmer.start;
                
                % Start timers
                tic_at_last_sample = tic;
                    
                % Keep going until trial is done
                while(  USE_PRACTICE_TRIAL_SYNCHRONOUS || ...
                       (~USE_PRACTICE_TRIAL_SYNCHRONOUS && (Nbeeps_delivered < Nbeeps_per_HBD_trial)) )

                    % Check whether it is time to sample the Shimmer again        
                    if( toc(tic_at_last_sample) >= sampling_period_check_Shimmer_sec + sampling_buffer_first_HB_sec )

                        % Get time when data arrive
                        time_current_packet_arrival_pre_GetSecs = GetSecs();
                        
                        % Read the latest data from shimmer data buffer, signalFormatArray defines the format of the data and signalUnitArray the unit
                        [newData, signalNameArray, signalFormatArray, signalUnitArray] = shimmer.getdata('c');

                        time_current_packet_arrival_post_GetSecs = GetSecs();

                        % Updated 2019/04/15
                        NnewData = size(newData,1);
                                
                        if( NnewData >= 2 )

                            if( ~FIRST_READ_HAS_OCCURRED  )
                                k_ECG_Signal_HBD = find(strcmp(string_ECG_signal_HBD, signalNameArray));

                                k_time = find(strcmp('Time Stamp', signalNameArray));
                                %k_ECG_LL_RA = find(strcmp('ECG LL-RA', signalNameArray));
                                %k_ECG_Signal_HBD = find(strcmp('ECG LA-RA', signalNameArray));

                                FIRST_READ_HAS_OCCURRED = true;
                                
                                % Reset this to zero now that it has been
                                % used
                                sampling_buffer_first_HB_sec = 0;
                            end

                            % Parse the data
                            try
                                X_time_current_packet_sec               = newData(:, k_time) / 1e3;
                                X_time_total_sec(end+1 : end+NnewData)  = newData(:, k_time) / 1e3;                                
                                Y_ECG_total(end+1 : end+NnewData)       = newData(:, k_ECG_Signal_HBD);
                                                                
                                Y_ECG_total_adj = ECG_adjust_baseline_spline( X_time_total_sec, Y_ECG_total, Minimum_RR_Interval_sec, Minimum_R_Prominence_mV );

                            catch err
                                fprintf('\n\n----------------------------\n')
                                warning('Error starting up');
                                err
                                fprintf('\n')
                                newData
                                size(newData);                                
                                NnewData
                                k_time
                                k_ECG_Signal_HBD
                                fprintf('\n\n----------------------------')
                                fprintf('\n----------------------------\n\n')
                                
                                % Skip to next loop iteration
                                
                                fprintf('\n** Error while loading data.');                                            
                                REPEAT_THIS_TRIAL = true;

                                if( USE_PTB_WINDOW )
                                    DS.addStep('INPUT', 'STEREO', filename_Rspike_or_beep_timeout, filename_Rspike_or_beep_timeout);
                                else
                                    fprintf('\n\n** Please stand by (electrodes are not playing nicely)');
                                end
                                
                                continue;
                            end

                            if( OUTPUT_DEBUG )
                                fprintf('\nData from %f to %f\t%f sec', ...
                                    X_time_total_sec(1), X_time_total_sec(end), X_time_total_sec(end) - X_time_total_sec(1));
                            end
                            
                            
                            %----------------------------------------------
                            % Real-time plot of trial [2019/04/30]
                            if( SHOW_REALTIME_PLOT_PER_TRIAL )
                                hold('off');
                                plot(X_time_total_sec, Y_ECG_total_adj, '-k');
                                grid('on');
                                title('Hold Shift + Esc to stop or reset trial');
                                ylabel('Ajusted ECG LA-RA (mV)');
                                xlabel('Time (sec)');
                                
                                % Update the plot
                                drawnow()
                            end
                            
                            
                            %------------------------------------------------------
                            % If keyboard button was pressed repeat the
                            % trial [2019/04/30]
                            
                            [ keyIsDown, seconds, keyCode ] = KbCheck;        
                            if( keyIsDown && keyCode(key_esc) && keyCode(key_shift) )
                                if( USE_PRACTICE_TRIAL_SYNCHRONOUS )
                                    % Done with the practice trial. Do not
                                    % attempt again
                                    USE_PRACTICE_TRIAL_SYNCHRONOUS = false;
                                    
                                    % So this doesn't look like it was a
                                    % repeated trial, as this was just
                                    % practice
                                    times_attempting_this_trial = times_attempting_this_trial - 1;
                                    fprintf('\n** Researcher pushed stop key. Moving to first trial.');
                                else
                                    fprintf('\n** Researcher pushed repeat key. Repeating this trial.');
                                    if( USE_PTB_WINDOW )
                                        DS.addStep('INPUT', 'STEREO', filename_researcher_retry, filename_researcher_retry);
                                    end
                                end

                                REPEAT_THIS_TRIAL = true;
                                
                                % Break out of the while loop
                                break;
                            end

                            %------------------------------------------------------
                            % Detection of all R spikes for this trial
                            % so far (between 0 and Ntones)

                            % Ensure there are enough data for peak detection                            
                            time_of_total_recording_sec = X_time_total_sec(end) - X_time_total_sec(1);
                            
                            % If there is at least one peak already
                            if( length(X_time_detected_peak_sec_array) >= 1 )
                                time_since_last_Rspike_sec = X_time_total_sec(end) - X_time_detected_peak_sec_array(end);
                            else
                                % No actual peaks yet so get total duration
                                % of recording
                                time_since_last_Rspike_sec = X_time_total_sec(end) - X_time_total_sec(1);
                            end
                            
                            if( HAVENT_DONE_DUMMY_BEEP )
                                % If the dummy spike hasn't been found yet,
                                % then just use total recording time
                                time_since_dummy_Rspike_sec = X_time_total_sec(end) - X_time_total_sec(1);
                            else
                                % If the dummy beep has been delivered then use
                                % time since that R spike
                                time_since_dummy_Rspike_sec = X_time_total_sec(end) - X_time_detected_peak_dummy_sec;
                            end
                            
                            
                            %fprintf('\nTime since last R spike = %0.3f sec', time_since_last_Rspike_sec);
                            
                            if( time_since_last_Rspike_sec > Maximum_RR_Interval_sec && ...
                                time_since_dummy_Rspike_sec > Maximum_RR_Interval_sec )
                                % There is a problem with peak
                                % detection on this trial, so
                                % stop the trial and start a
                                % new one
                                fprintf('\n** R spike finder wait timed out. Repeating this trial.');                                            
                                REPEAT_THIS_TRIAL = true;

                                if( USE_PTB_WINDOW )
                                    DS.addStep('INPUT', 'STEREO', filename_Rspike_or_beep_timeout, filename_Rspike_or_beep_timeout);
                                else
                                    fprintf('\n\n** Please stand by (electrodes are not playing nicely)');
                                end

                                % Break out of the while loop
                                break;
                            
                            elseif( time_of_total_recording_sec > Minimum_RR_Interval_sec )
                                % Peak detection
                                [peak_Y_array, peak_X_array, peak_width_array, peak_prom_array] = ...
                                    findpeaks(Y_ECG_total_adj, X_time_total_sec, ...
                                    'MinPeakDistance', Minimum_RR_Interval_sec, ...
                                    'MinPeakProminence', Minimum_R_Prominence_mV);
                                
                                % Update plot
                                if( SHOW_REALTIME_PLOT_PER_TRIAL )
                                    hold('all');
                                    plot(peak_X_array, peak_Y_array, 'or', 'LineWidth', 3);
                                
                                    % Update the plot
                                    drawnow()
                                end

                                % Narrow to the valid peaks based on
                                % prominence
                                %k_valid_peaks = peak_prom_array > Minimum_R_Prominence_mV;
                                %peak_Y_array = peak_Y_array(k_valid_peaks);
                                %peak_X_array = peak_X_array(k_valid_peaks);
                                %peak_prom_array = peak_prom_array(k_valid_peaks);

                                Npeaks_found = length(peak_X_array);

                                if( OUTPUT_DEBUG )
                                    fprintf('\n\tFound %d peaks', length(peak_X_array));
                                end

                                % First peak disvorey has to be ONE
                                % peak (not more than that). If there
                                % are 2 peaks on first search, then the
                                % first one doesn't count
                                if( ~FOUND_PEAKS_ON_THIS_TRIAL && Npeaks_found > 0 )                                        
                                    Npeaks_found_that_do_not_count = Npeaks_found - 1;
                                    FOUND_PEAKS_ON_THIS_TRIAL = true;

                                    if( OUTPUT_DEBUG )
                                        fprintf('\n\tInitially found %d peaks so ignoring %d of them', ...
                                            Npeaks_found, Npeaks_found_that_do_not_count);
                                    end
                                end

                                %------------------------------------------------------
                                % If (1) it has been enough time since last
                                % R spike and (2) there is a new peak that
                                % hasn't previously been counted
                                %if( (toc(tic_at_last_Rspike) >= Minimum_RR_Interval_sec) && Npeaks_found > Npeaks_counted )
                                if( Npeaks_found - Npeaks_found_that_do_not_count > Npeaks_counted )
                                    
                                    Npeaks_counted = Npeaks_counted + 1;
                                    %ALREADY_BEEPED_FOR_THIS_RSPIKE = false;

                                    % Set up to present a beep for this
                                    % R spike at the proper time
                                    if( OUTPUT_DEBUG )
                                        fprintf('\n\t\tCounted a new peak (number %d): Peak at %f sec', ...
                                            Npeaks_counted, peak_X_array(end));
                                    end

                                    % Time from packet start to R spike
                                    time_from_packet_start_to_peak_sec = peak_X_array(end) - X_time_current_packet_sec(1);
                                    time_beep_intended_sec = time_current_packet_arrival_pre_GetSecs + time_from_packet_start_to_peak_sec + delay_HB_sec;
                                    
                                    % If it LOOKS like we are waiting TOO
                                    % LONG to beep, then just set it to the
                                    % max waiting duration (the delay
                                    % itself)
                                    time_to_beep_sec = time_beep_intended_sec - GetSecs();
                                    if( time_to_beep_sec > delay_HB_sec )
                                        time_beep_intended_sec = time_beep_intended_sec - (time_to_beep_sec - delay_HB_sec);
                                    end
                                    
                                    if( OUTPUT_DEBUG )
                                        fprintf('\n\nAll peak times (Shimmer time; sec): ');
                                        fprintf('%0.3f\t', peak_X_array);
                                        fprintf('\nTime of current packet arrival (Shimmer time): %0.3f', X_time_current_packet_sec(1));                                    
                                        fprintf('\nPRE Packet arrival (GetSecs time): %0.3f', time_current_packet_arrival_pre_GetSecs);
                                        fprintf('\nPOST Packet arrival (GetSecs time): %0.3f', time_current_packet_arrival_post_GetSecs);
                                        fprintf('\nTime from packet start to most recent peak (sec): %0.3f', time_from_packet_start_to_peak_sec);
                                        fprintf('\nHB delay (sec): %0.3f', delay_HB_sec);
                                        fprintf('\nTime of intended beep (GetSecs): %0.3f', time_beep_intended_sec);
                                        current_time_GetSecs = GetSecs();
                                        fprintf('\nCurrent time (GetSecs): %0.3f', current_time_GetSecs);
                                        fprintf('\nTime to beep (sec): %0.3f', time_beep_intended_sec - current_time_GetSecs)
                                    end
                                    

                                    % Queue up the beep!
                                    %PsychPortAudio('FillBuffer', pahandle, tone_sound_array);
                                    %time_beep_actual_sec = PsychPortAudio('Start', pahandle, repetitions, time_beep_intended_sec, 1);
                                    %waitForEndOfPlayback = 1;
                                    %PsychPortAudio('Stop', pahandle, waitForEndOfPlayback);
                                    %duration_beep_late = time_beep_actual_sec - time_beep_intended_sec;
                                    %fprintf('\n\t\tBeep %d of %d was late by %0.3f msec', Nbeeps_delivered, Nbeeps_per_HBD_trial, 1e3*duration_beep_late);

                                    if( HAVENT_DONE_DUMMY_BEEP )
                                        
                                        fprintf('\n\t\tDUMMY non-beep in %0.3f msec', 1e3 * (time_beep_intended_sec - GetSecs()) );

                                        time_beep_actual_sec = PsychPortAudio('Start', pahandle_dummy, repetitions, time_beep_intended_sec, 1);
                                        duration_beep_late = time_beep_actual_sec - time_beep_intended_sec;
                                        fprintf('\n\t\t\tDUMMY non-beep was late by %0.3f msec', 1e3*duration_beep_late);
                                        
                                        HAVENT_DONE_DUMMY_BEEP = false;
                                        
                                        X_time_detected_peak_dummy_sec = peak_X_array(end);
                                        
                                    else
                                        fprintf('\n\t\tBeeping in %0.3f msec', 1e3 * (time_beep_intended_sec - GetSecs()) );
                                        
                                        time_to_beep_sec = (time_beep_intended_sec - GetSecs()) / 1e3;
                                        
                                        if( time_to_beep_sec > MAX_TIME_TO_BEEP_SEC )
                                            % There is a problem with peak
                                            % detection on this trial, so
                                            % stop the trial and start a
                                            % new one
                                            fprintf('\n** Beep wait timed out. Repeating this trial.');                                            
                                            REPEAT_THIS_TRIAL = true;
                                            
                                            DS.addStep('INPUT', 'STEREO', filename_Rspike_or_beep_timeout, filename_Rspike_or_beep_timeout);
                                            
                                            % Break out of the while loop
                                            break;
                                        end

                                        time_beep_actual_sec = PsychPortAudio('Start', pahandle, repetitions, time_beep_intended_sec, 1);
                                        duration_beep_late = time_beep_actual_sec - time_beep_intended_sec;
                                        fprintf('\n\t\t\tBeep %d of %d was late by %0.3f msec', Nbeeps_delivered, Nbeeps_per_HBD_trial, 1e3*duration_beep_late);

                                        Nbeeps_delivered = Nbeeps_delivered + 1;                                        
                                        
                                        % Save time of beep
                                        X_time_beep_sec_array(end+1)            = time_beep_actual_sec;
                                        X_time_beep_late_sec_array(end+1)       = duration_beep_late;                                        
                                        X_time_detected_peak_sec_array(end+1)   = peak_X_array(end);                                        
                                        Y_ECG_detected_peak_mV_array(end+1)     = peak_Y_array(end);                                        
                                        Y_peak_prominance_array(end+1)          = peak_prom_array(end);
                                        
                                        %fprintf('\nActual time of beep (PsychPort): %f', time_beep_actual_sec);
                                        %fprintf('\nBeep lateness (sec): %f', duration_beep_late);
                                    end
                                end
                            end
                        end
                        tic_at_last_sample = tic;
                    end
                end
                shimmer.stop;
                

                if( REPEAT_THIS_TRIAL )
                    times_attempting_this_trial = times_attempting_this_trial + 1;

                    % Next iteration of for loop
                    continue;
                end      
                

            case 'DURATION'
                fprintf('\nTrial duration: %0.0f sec', duration_HBD_trial);
                
                if( USE_PTB_WINDOW )
                    % Wait some time for duration of trial
                    DS.addStep('DELAY', 'STEREO', filename_trial_underway, filename_trial_underway, duration_HBD_trial);
                else
                    fprintf('\n\nTrial underway for total duration');
                end
        end

        % After trial is done (delay has expired) - get input from PP
        if( USE_PTB_WINDOW )
            DS.addStep('INPUT', 'STEREO', filename_request_input, filename_request_input);
        else
            fprintf('\n\n** Did the beeps occur DURING your heartbeats? (Yes=1, No=2)\n')
            
            tic_answer = tic();
            while( ~keyCode(key_conf_low) && ~keyCode(key_conf_high) )
                % Check the state of the keyboard.
                [ keyIsDown, seconds, keyCode ] = KbCheck(KbDevice_Number);
            end
            RT_answer_sec = toc(tic_answer);
            fprintf('\n');

            if( keyCode(key_yes) )
                userInput_value_array(t) = 1;
                userInput_RT_array(t) = RT_answer_sec;

            elseif( keyCode(key_no) )
                userInput_value_array(t) = 2;
                userInput_RT_array(t) = RT_answer_sec;
            else
                error('invalid key, somehow');
            end
        end
        
        if( USE_PTB_WINDOW && RATE_CONFIDENCE_W_SLIDER_SCALE )
            % Slider scale for confidence
            DS.addNewSliderScale( sliderScale_Confidence );
            
            % Start displaying the steps
            DS.nextStepDisplay();
            
            numQuestions = 2;
        else            
            if( USE_PTB_WINDOW )
                % After trial is done (delay has expired) - get input from PP
                DS.addStep('DELAY', 'STEREO', filename_rate_confidence, filename_rate_confidence, 0.01);
            
                % Start displaying the steps
                DS.nextStepDisplay();
                
                numQuestions = 1;
            end
            
            % Wait until all keys are released
            while KbCheck; end
            [ keyIsDown, seconds, keyCode ] = KbCheck(KbDevice_Number);
            
            fprintf('\n** How would you rate your confidence (`=Guess, 1=Low, 2=High)?\n');
            tic_answer = tic();
            while( ~keyCode(key_conf_guess) && ~keyCode(key_conf_low) && ~keyCode(key_conf_high) )
                % Check the state of the keyboard.
                [ keyIsDown, seconds, keyCode ] = KbCheck(KbDevice_Number);
            end
            RT_answer_sec = toc(tic_answer);
            fprintf('\n');

            if( keyCode(key_conf_guess) )
                confidence_rating = 0;

            elseif( keyCode(key_conf_low) )
                confidence_rating = 0.5;

            elseif( keyCode(key_conf_high) )
                confidence_rating = 1;
            else
                error('invalid key, somehow');
            end
            
            numQuestions = 0;
        end

        
        
        %------------------------------------------------------------------
        % Get input from PP for the most recent questions
        if( USE_PTB_WINDOW )
            for recencyIndex = 1:numQuestions
                % Get response from question number 1, 2, 3, 4, etc.
                [response_value_array(numQuestions-recencyIndex+1), ...
                 response_time_array(numQuestions-recencyIndex+1)] = ...
                 DS.getRecentResponse( recencyIndex );
            end

            userInput_value_array(t) = response_value_array(1);
            userInput_RT_array(t) = response_time_array(1);
        else
            response_value_array = [];
        end

        if( RATE_CONFIDENCE_W_SLIDER_SCALE )
            confidence_array(t) = response_value_array(2);
            confidence_RT_array(t) = response_time_array(2);
        else
            confidence_array(t) = confidence_rating;
            confidence_RT_array(t) = RT_answer_sec; 
        end

        fprintf('\n\tResponse\t%d\tin\t%f sec', userInput_value_array(t), userInput_RT_array(t));
        fprintf('\n\tConfidence\t%f\tin\t%f sec', confidence_array(t), confidence_RT_array(t));

        % User wants to exit
        if( any(response_value_array == -1) || any(userInput_value_array==-1) || any(confidence_array == -1) )
            fprintf('\n\nUser exits');            
            break;
        end

        % Determine if the user responded as "coincident" or "not"
        userInput_Coincident_array(t) = ...
             userInput_value_array(t) == responseValue_coincident;

        % Determine whether PP is correct
        userIsCorrect_array(t) = ...
            userInput_Coincident_array(t) == trialIsCoincident;


        % Provide knowledge of result
        % Determine what filename to use to display the approprite result
        %  The user can be CORRECT or INCORRECT
        %  and the trial can be COINCIDENT or NOT COINCIDENT
        if( userIsCorrect_array(t) )            
            fprintf('\n\tPP is CORRECT');

            if( trialIsCoincident )                
                filename_result = filename_correct_coincident;  
                Nhits = Nhits + 1;
            else
                filename_result = filename_correct_not_coincident; 
                Ncorrect_rejections = Ncorrect_rejections + 1;
            end

        else
            fprintf('\n\tPP is INCORRECT');

            if( trialIsCoincident )
                filename_result = filename_incorrect_coincident;
                Nmisses = Nmisses + 1;
            else
                filename_result = filename_incorrect_not_coincident;
                Nfalse_alarms = Nfalse_alarms + 1;
            end            
        end        

        % Present knowledge of results after the trial?
        if( TrainingMode )
            if( USE_PTB_WINDOW )
                % Display the result (determined above) and wait some delay before proceeding
                DS.addStep('DELAY', 'STEREO', filename_result, filename_result, delay_HB_result_sec);
                DS.nextStepDisplay();
            end
        end
        
        %----------------------------------------------------------
        % Show progress to PP during experiment?
        if( experiment_specs.ShowProgress )

            % Calculate how much the PP has completed
            fraction_completed = double(t) / double(Ntrials_HBD);

            % If PP is more than 75% done and hasn't yet seen the slide
            % that show they are more than 75% done...            
%             if( fraction_completed >= 0.75 && ~progress_75_percent_was_shown )
% 
%                 % Show the slide that the PP is 75% done
%                 DS.addStep('DELAY', 'STEREO', filename_progress_75_percent, [], delay_break_sec);
%                 DS.addStep('NAVIGATE_FORWARD_ANYKEY', 'STEREO', filename_progress_75_percent, []);
%                 progress_75_percent_was_shown = true;
%                 progress_50_percent_was_shown = true;
%                 progress_25_percent_was_shown = true;
%                 DS.nextStepDisplay();
%             end

            if( fraction_completed >= 0.5 && ~progress_50_percent_was_shown )

                % Show the slide that the PP is 50% done
                if( USE_PTB_WINDOW )
                    DS.addStep('DELAY', 'STEREO', filename_progress_50_percent, [], delay_break_sec);
                    DS.addStep('NAVIGATE_FORWARD_ANYKEY', 'STEREO', filename_progress_50_percent, []);
                    DS.nextStepDisplay();
                else
                    fprintf('\n\n** Participant is more than half done');
                end
                progress_50_percent_was_shown = true;
                progress_25_percent_was_shown = true;
                
            end

%             if( fraction_completed >= 0.25 && ~progress_25_percent_was_shown )
% 
%                 % Show the slide that the PP is 25% done
%                 DS.addStep('DELAY', 'STEREO', filename_progress_25_percent, [], delay_break_sec);
%                 DS.addStep('NAVIGATE_FORWARD_ANYKEY', 'STEREO', filename_progress_25_percent, []);
%                 progress_25_percent_was_shown = true;
%                 DS.nextStepDisplay();
%             end
            
            if( USE_PTB_WINDOW )
                % Draw the text with trial number
                oldTextSize = Screen('TextSize', window, round(fontsize_questions/4));
                DrawFormattedText(window, trial_num_string, 'right', 0.99*rect(4), 128*[1 1 1]);
            end
        end

        %------------------------------------------------------------------
        % Write results to file        
        fprintf(FILE, '\n%d\t%d\t%f\t%d\t%d\t%f\t%d\t%f\t%f', ...
            t, ...
            times_attempting_this_trial, ...
            delay_HB_sec_array(t), ...
            trialIsCoincident_array(t), ...
            userInput_Coincident_array(t), ...
            userInput_RT_array(t), ...
            userIsCorrect_array(t), ...
            confidence_array(t), ...
            confidence_RT_array(t) ...
            );
        %------------------------------------------------------------------

        
        
        %----------------------------------------------------------
        % Write a plot showing HB data and times of R spikes and
        % beeps
        if( WRITE_FIGURE_PER_TRIAL )

            if( ~exist('hfig') )
                hfig = figure('visible', 'off');
            else
                figure(hfig);
            end
            hold('off');

            % Plot ECG
            plot(X_time_total_sec, Y_ECG_total_adj, '-k');
            hold('all');                

            % Plot detected peaks
            plot(X_time_detected_peak_sec_array, Y_ECG_detected_peak_mV_array, 'or');

            % Plot beep times
            plot(X_time_detected_peak_sec_array + delay_HB_sec + X_time_beep_late_sec_array, Y_ECG_detected_peak_mV_array, 'sb');

            get(hfig, 'YLim');

            for b = 1:length(X_time_beep_late_sec_array)
                if( b == 1 )
                    extra_text = 'Beep lateness (msec): ';
                    extra_text2 = 'Prominence (mV): ';
                else
                    extra_text = '';
                    extra_text2 = '';
                end

                text(X_time_detected_peak_sec_array(b), min(Y_ECG_total_adj), ...
                    sprintf('%s%0.0f\n%s%0.1f', extra_text, X_time_beep_late_sec_array(b), extra_text2, Y_peak_prominance_array(b)), ...
                    'HorizontalAlignment', 'Right', 'FontSize', 6 );
            end

            title(sprintf('%s, Trial %d, Attempt %d, Delay=%0.1f ms, Worst lateness = %0.2f ms', ...
                PPid, t, times_attempting_this_trial, delay_HB_msec, max(X_time_beep_late_sec_array)));
            xlabel('Time (sec)');
            ylabel('Adjusted ECG (mv)');
            print(hfig, sprintf('%s/out-ECG-Trial_%02d-Attempt_%d.png', output_dir, t, times_attempting_this_trial), '-dpng');                
            %close(hfig);
        end
        
        if( USE_PTB_WINDOW )
            % Flush textures
            DS.flushAllTextures();
        end

%         % If this is NOT the last trial
%         if( t ~= Ntrials_HBD )
% 
%             % Proceed to next trial - input from PP
%             DS.addStep('NAVIGATE_FORWARD_ANYKEY', 'STEREO', filename_proceed_to_trial, filename_proceed_to_trial);
%             DS.nextStepDisplay();
%             
%             [newX,newY,textHeight] = Screen('DrawText', window, sprintf('Finished trial %d    ', t), 0, 0, 128*[ 1 1 1 ] );
%             Screen('Flip', window);
%         end

        t = t + 1;
    end

    if( USE_PTB_WINDOW )
        % The task is over
        DS.addStep('NAVIGATE', 'STEREO', filename_task_is_over, filename_task_is_over);
        DS.nextStepDisplay();    
    else
        fprintf('\n\nCongratulations the task is complete');
    end
    
    
    % Stop timing the task
    timer_taskEnd = toc(timer_taskStart);



    %% Stop hardware communication
    shimmer.disconnect;
    Snd('close');

    %% Write results to file

    fprintf('\n______________________________________________________');
    fprintf('\nResults [%s]', datestr(now));

    % Determine the actual number of trials completed (in case user had
    % exited early)
    Ntrials_HBD_completed = length(userInput_value_array);

    % If the user exited early, then reduce the total number of trials
    % comleted (the last trial was an "abort")
    if( Ntrials_HBD_completed < Ntrials_HBD )
        Ntrials_HBD_completed = Ntrials_HBD_completed - 1;
    end

    % Percent correct
    fractionCorrect = sum(userIsCorrect_array) / Ntrials_HBD_completed;

    % Signal detection stuff    

    % arcsin*sqrt(P(A))

    % The for loop will output to the command window (=1), and then to the FILE
    % stream (=FILE)
    for OUTPUTSTREAM = [1, FILE]
        fprintf(OUTPUTSTREAM, '\n\nBasic Results');
        fprintf(OUTPUTSTREAM, '\nNumTrials\t%d', Ntrials_HBD_completed);
        fprintf(OUTPUTSTREAM, '\nFractionCorrect\t%f', fractionCorrect);        

        fprintf(OUTPUTSTREAM, '\n\nSignal Detection Results');
        fprintf(OUTPUTSTREAM, '\nNumHits\t%d', Nhits);
        fprintf(OUTPUTSTREAM, '\nNumMisses\t%d', Nmisses);
        fprintf(OUTPUTSTREAM, '\nNumFalseAlarms\t%d', Nfalse_alarms);
        fprintf(OUTPUTSTREAM, '\nNumCorrectRejections\t%d', Ncorrect_rejections);

        fprintf(OUTPUTSTREAM, '\n\nHBD - task duration: %0.2f min', timer_taskEnd/60);
    end

    % Close the log file
    fclose(FILE);

    % Close the onscreen window
    Priority(0);
    Screen('CloseAll')
    PsychPortAudio('Close', pahandle);


    %% Stop logging command window output
    fprintf('\nDone!\n');
    diary('off');    

    %fprintf('The percentage of received packets: %d \n', shimmer.getpercentageofpacketsreceived(timeStamp)); % Detect loss packets

end