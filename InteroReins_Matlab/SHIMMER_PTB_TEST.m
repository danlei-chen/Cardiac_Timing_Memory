%% Connect to Shimmer device

pa.comPort_ECG = sprintf('%d', 8);
pa.sampling_rate_ECG_Hz = 1024;
% This is how much time between checking the Shimmer for data (smaller
% numbers increase temporal precision but risk an error if there are no
% data available)
pa.sampling_period_check_Shimmer_sec = 0.050;
% For detecting R spikes
MAX_TIME_TO_BEEP_SEC = 0.7;
Minimum_RR_Interval_sec = 0.5;
Minimum_R_Prominence_mV = 0.6; % run resting state and see 
Npeaks_found_that_do_not_count = 0;
Npeaks_counted = 0;
delay_HB_sec = 0;

% Connect to Shimmer device
if strcmpi(pa.shimmer, 'true')
    % Connect to Shimmer at specified COM port
    fprintf('\n\nAttempting to connect to Shimmer, this may take 10 sec...\n');
    shimmer = ShimmerHandleClass(pa.comPort_ECG);

    % assign user friendly macros for setenabledsensors
    SensorMacros = SetEnabledSensorsMacrosClass;                               

    % Ensure connection
    if (shimmer.connect)

        fprintf('\n\nSuccessfully connected to Shimmer device!\n');

        % Select sampling rate
        shimmer.setsamplingrate(pa.sampling_rate_ECG_Hz); 

        % Select internal expansion board; select 'ECG' to enable both SENSOR_EXG1 and SENSOR_EXG2 
        shimmer.setinternalboard('ECG');                                       

        % Disable other sensors
        shimmer.disableallsensors;                                             

        % Enable SENSOR_EXG1 and SENSOR_EXG2
        shimmer.setenabledsensors(SensorMacros.ECG,1);
        
    else
        error('Could not connect to Shimmer device');
    end
end

%%
clear all;
close all;

%%%%%%% read in the stimuli generator
addpath(genpath('/Applications/Psychtoolbox/'))
% addpath(genpath('~/Documents/MATLAB/Psychtoolbox'))
Screen('Preference', 'SkipSyncTests', 1); 
rng('shuffle')
screenNum = max(Screen('Screens'));
% [wPtr, rect] = Screen('OpenWindow', screenNum);
[wPtr,rect] = Screen('OpenWindow',0,0,[0 0 800 600]);
Screen('FillRect', wPtr, [0 0 0], rect)

%%% Read in pics
load faces_stim.mat
load houses_stim.mat

images=cat(3,face(:,:,1),...
    house(:,:,1),...
    face(:,:,2),...
    face(:,:,3),...
    house(:,:,2),...
    house(:,:,3));
start_time = [600, 7000, 14000, 21000, 28500, 35000];
half_beat=500;
ITI = 1000;
presentation_duration=200;
repetitions=2;

presentation_duration=presentation_duration/1000;
start_time=start_time/1000;
half_beat=half_beat/1000;
ITI = ITI/1000;

block_start_time = GetSecs; 
trial_start=block_start_time + 0.04;

%%
% Start getting data from Shimmer
shimmer.start;
% Start timers
tic_at_last_sample = tic;

sampling_period_check_Shimmer_sec = 0.05;
sampling_buffer_first_HB_sec = 0.1;
FIRST_READ_HAS_OCCURRED = false;
string_ECG_signal_HBD = 'ECG LA-RA';

X_time_total_sec = [];
Y_ECG_total = [];
X_time_beep_sec_array = [];
X_time_beep_late_sec_array = [];
X_time_detected_peak_sec_array = [];
Y_ECG_detected_peak_mV_array = [];
Y_peak_prominance_array = [];

WaitSecs(2);

%%
for n=1:size(images,3)
    fixation = ['+'];
    DrawFormattedText(wPtr, fixation, 'center', 'center');
    Screen('TextSize', wPtr, 80);
    Screen('TextColor', wPtr, [255,255,255]);
    Screen('Flip',wPtr);
    
    for r =1:repetitions
        
        % Check whether it is time to sample the Shimmer again        
        if( toc(tic_at_last_sample) >= sampling_period_check_Shimmer_sec + sampling_buffer_first_HB_sec )

            % Get time when data arrive
            time_current_packet_arrival_pre_GetSecs = GetSecs();

            % Read the latest data from shimmer data buffer, signalFormatArray defines the format of the data and signalUnitArray the unit
            [newData, signalNameArray, signalFormatArray, signalUnitArray] = shimmer.getdata('c');

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
                X_time_current_packet_sec               = newData(:, k_time) / 1e3;
                X_time_total_sec(end+1 : end+NnewData)  = newData(:, k_time) / 1e3;                                
                Y_ECG_total(end+1 : end+NnewData)       = newData(:, k_ECG_Signal_HBD);

                Y_ECG_total_adj = ECG_adjust_baseline_spline( X_time_total_sec, Y_ECG_total, Minimum_RR_Interval_sec, Minimum_R_Prominence_mV );

                % Real-time plot of trial [2019/04/30]
                hold('off');
                plot(X_time_total_sec, Y_ECG_total_adj, '-k');
                grid('on');
                ylabel('Ajusted ECG LA-RA (mV)');
                xlabel('Time (sec)');

                % Update the plot
                drawnow()

                %------------------------------------------------------
                % Detection of all R spikes for this trial
                % so far (between 0 and Ntones)

                % Ensure there are enough data for peak
                % detection                            
                time_of_total_recording_sec = X_time_total_sec(end) - X_time_total_sec(1);

                % If there is at least one peak already
                if( length(X_time_detected_peak_sec_array) >= 1 )
                    time_since_last_Rspike_sec = X_time_total_sec(end) - X_time_detected_peak_sec_array(end);
                else
                    % No actual peaks yet so get total duration
                    % of recording
                    time_since_last_Rspike_sec = X_time_total_sec(end) - X_time_total_sec(1);
                end

                if( time_of_total_recording_sec > Minimum_RR_Interval_sec )
                    % Peak detection
                    [peak_Y_array, peak_X_array, peak_width_array, peak_prom_array] = ...
                        findpeaks(Y_ECG_total_adj, X_time_total_sec, ...
                        'MinPeakDistance', Minimum_RR_Interval_sec, ...
                        'MinPeakProminence', Minimum_R_Prominence_mV);

                    % Update plot
                    hold('all');
                    plot(peak_X_array, peak_Y_array, 'or', 'LineWidth', 3);

                    % Update the plot
                    drawnow()

                    Npeaks_found = length(peak_X_array);
                    %------------------------------------------------------
                    % If (1) it has been enough time since last R spike and (2) there is a new peak that hasn't previously been counted
                    %if( (toc(tic_at_last_Rspike) >= Minimum_RR_Interval_sec) && Npeaks_found > Npeaks_counted )
                    if( Npeaks_found - Npeaks_found_that_do_not_count > Npeaks_counted )

                        Npeaks_counted = Npeaks_counted + 1;
                        %ALREADY_BEEPED_FOR_THIS_RSPIKE = false;

                        % Set up to present a beep for this
                        % R spike at the proper time

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
                    end

                    time_to_beep_sec = (time_beep_intended_sec - GetSecs()) / 1e3;

                    if( time_to_beep_sec > MAX_TIME_TO_BEEP_SEC )
                        % There is a problem with peak detection on this trial, so stop the trial and start a new one
                        fprintf('\n** Beep wait timed out. Repeating this trial.');                                            

                        % Break out of the while loop
                        break;
                        sca;
                    end

        %             time_beep_actual_sec = PsychPortAudio('Start', pahandle, repetitions, time_beep_intended_sec, 1);
        %             duration_beep_late = time_beep_actual_sec - time_beep_intended_sec;
        %             fprintf('\n\t\t\tBeep %d of %d was late by %0.3f msec', Nbeeps_delivered, Nbeeps_per_HBD_trial, 1e3*duration_beep_late);

        %             Nbeeps_delivered = Nbeeps_delivered + 1;                                        

                    % Save time of beep
        %             X_time_beep_sec_array(end+1)            = time_beep_actual_sec;
        %             X_time_beep_late_sec_array(end+1)       = duration_beep_late;                                        
                    X_time_detected_peak_sec_array(end+1)   = peak_X_array(end);                                        
                    Y_ECG_detected_peak_mV_array(end+1)     = peak_Y_array(end);                                        
                    Y_peak_prominance_array(end+1)          = peak_prom_array(end);

                    %fprintf('\nActual time of beep (PsychPort): %f', time_beep_actual_sec);
                    %fprintf('\nBeep lateness (sec): %f', duration_beep_late);

                end
            end
        end
    
    
    
        image_to_show=images(:,:,n);
        image_absolute_onset_time=trial_start + half_beat*(r-1);
        image_absolute_offset_time=image_absolute_onset_time+presentation_duration;

        [imageHeight, imageWidth, colorChannels] = size(image_to_show);

        pointer_to_offscreen_window_with_image_on_it = Screen('MakeTexture',wPtr,image_to_show);
        Screen('DrawTexture',wPtr,pointer_to_offscreen_window_with_image_on_it);
        
        fixation = ['+'];
        DrawFormattedText(wPtr, fixation, 'center', 'center');
        Screen('TextSize', wPtr, 80);
        Screen('TextColor', wPtr, [255,255,255]);

        [~, stimOnset] = Screen('Flip',wPtr,image_absolute_onset_time);
        
        [~, stimOffset] = Screen('Flip',wPtr, image_absolute_offset_time);
        
        fixation = ['+'];
        DrawFormattedText(wPtr, fixation, 'center', 'center');
        Screen('TextSize', wPtr, 80);
        Screen('TextColor', wPtr, [255,255,255]);
        Screen('Flip',wPtr);
    end
    
    trial_start = image_absolute_offset_time + ITI;
    
    
    

end

%% Stop hardware communication
shimmer.disconnect;






    