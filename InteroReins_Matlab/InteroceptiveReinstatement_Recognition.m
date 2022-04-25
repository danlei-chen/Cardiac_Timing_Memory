% Set up dir
if (ismac)
    output_dir = sprintf('%ssubj-%s_InteroRein_recognition_output/',base_dir,pa.subjid);
elseif (IsWin)
    output_dir = sprintf('%ssubj-%s_InteroRein_recognition_output\\',base_dir,pa.subjid);
end
mkdir(output_dir);

Screen('Preference', 'SkipSyncTests', 1); 
screenNum = 1;
[wPtr, rect] = Screen('OpenWindow', screenNum);
% [wPtr,rect] = Screen('OpenWindow',screenNum,0,[0 0 800 600]);

%black background
Screen('FillRect', wPtr, [0 0 0], rect);
% HideCursor; 
Screen('TextFont', wPtr, 'Helvetica'); 
xCenter = rect(3)/2;
yCenter = rect(4)/2;

clc;

%draw fixation
Screen('Flip',wPtr);
Screen('TextSize', wPtr, 60);
waitMessage = ['Please wait... Connecting...'];
DrawFormattedText(wPtr, waitMessage, 'center', 'center',[255,255,255]);
Screen('Flip',wPtr);

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
        sca;
        error('Could not connect to Shimmer device');
    end
 
    % ECG Initialize
    % This is how much time between checking the Shimmer for data (smaller
    % numbers increase temporal precision but risk an error if there are no
    % data available)
    sampling_period_check_Shimmer_sec = 0.050;
 
    % Keep this small to let the program run fast
    Plot_viewable_duration_sec = Inf;
 
    % ECG channel used to detect R spikes
    string_ECG_signal_HBD = 'ECG LA-RA';
    
    % ECG output files
    FILE_OUT_ECG = fopen(sprintf('%srecognition-ECG.csv',output_dir),'w');
    % make header
    cHeader = {'Time(Sec)' 'ECG_Raw(mV)'}; %dummy header
    commaHeader = [cHeader;repmat({','},1,numel(cHeader))]; %insert commaas
    commaHeader = commaHeader(:)';
    textHeader = cell2mat(commaHeader); %cHeader in text with commas
    % write header to file
    fprintf(FILE_OUT_ECG,'%s\n',textHeader);
    
    % ECG output files
    FILE_OUT_ECG_Adj = fopen(sprintf('%srecognition-ECG-Baseline_Adjusted.csv',output_dir),'w');
    % write header to file
    fprintf(FILE_OUT_ECG_Adj,'%s\n',textHeader);
 
end

%% recognition
% start instuctions
keyIsDown = 0;
indexPressed = 0;
rand_key = randi([1 10],1);
while keyIsDown == 0 && indexPressed == 0
    [keyIsDown, pressedSecs, keyCodes] = KbCheck(-1);

    if mod(rand_key,2) == 0 %(M)old (N)new
        instructions = 'In this part of the experiment, you will see some images.\n\n After viewing the image, if you think the image is old \nand you have seen it from the previous part of the experiment, \nplease press the (M) button as fast as possible. \n\nIf you think the image is new, \nplease press the (N) button as fast as possible. \n\n\nPlease let the experimenter know now if you have any questions.\n Press either (N) or (M) button to begin.';
        key_press_label = ["old";"new"];
    else %(M)new (N)old
        instructions = 'In this part of the experiment, you will see some images.\n\n After viewing the image, if you think the image is old \nand you have seen it from the previous part of the experiment, \nplease press the (N) button as fast as possible. \n\nIf you think the image is new, \nplease press the (M) button as fast as possible. \n\n\nPlease let the experimenter know now if you have any questions.\n Press either (N) or (M) button to begin.';
        key_press_label = ["new";"old"];
    end
    key_press_button=['m';'n'];
    
    Screen('TextSize', wPtr, 50);
    DrawFormattedText(wPtr, instructions, 'center', 'center',[255,255,255]);
    Screen('Flip',wPtr);
    
    if keyIsDown == 1
        indexPressed = 1;
        keyIsDown = 0;
    end           
end
Screen('Flip',wPtr);

% Start getting data from Shimmer
shimmer.start;
seq.rec.block_start_time=GetSecs();
disp('*********************************************');
disp(['Start of Recognition Block '])
disp('*********************************************');

seq.rec.trial_start_time = nan(1,pa.total_image_used);
seq.rec.stimOnset = nan(1,pa.total_image_used);
seq.rec.real_trial_onset_relative_to_block_start_time = nan(1,pa.total_image_used);

RestrictKeysForKbCheck(active_keys);
rsp.rec.response = nan(1,pa.total_image_used);
rsp.rec.RT = nan(1,pa.total_image_used);
rsp.rec.keyName = nan(1,pa.total_image_used);
rsp.rec.accuracy = nan(1,pa.total_image_used);

%trial level
for trial = 1:pa.total_image_used
    disp(['trial ' num2str(trial) ]);
    seq.rec.trial_start_time(trial) = GetSecs();
    
    image_to_show = seq.rec.trial(:,:,trial);
    [imageHeight, imageWidth, colorChannels] = size(image_to_show);
    pointer_to_offscreen_window_with_image_on_it = Screen('MakeTexture',wPtr,image_to_show);

    %draw fixation
    fixation = '+';
    Screen('TextSize', wPtr, 80);
    DrawFormattedText(wPtr, fixation, 'center', 'center');
    Screen('Flip',wPtr);
    
    label_this_trial = seq.rec.trial_label(trial);
    in_block_ITI_this_trial = seq.rec.in_block_ITI(trial);
    cardiac_timing_this_trial = seq.rec.cardiac_timing_sequence(trial);
    novelty_marker_this_trial = seq.rec.trial_novelty_marker(trial);
    congruency_marker_this_trial = seq.rec.trial_congruency_marker(trial);
 
    % Initialize timers and other variables

    Npeaks_counted = 0;
    Nbeeps_delivered = 0;
    Npeaks_found_that_do_not_count = 0;
    FIRST_READ_HAS_OCCURRED = false;
    FOUND_PEAKS_ON_THIS_TRIAL = false;

    X_time_total_sec = [];
    Y_ECG_total = [];
    X_time_beep_sec_array = [];
    X_time_beep_late_sec_array = [];   
    X_time_detected_peak_sec_array = [];
    Y_ECG_detected_peak_mV_array = [];
    Y_peak_prominance_array = [];
    all_peaks_X = [];
    all_peaks_Y = [];

    tic_at_last_sample = tic;
    % Reset sampling bufer for first heartbeat
    sampling_buffer_first_HB_sec = 0.1;

%     hfig = figure;
    
%     while Npeaks_counted < pa.rec.image_repetition + in_block_ITI_this_trial
    while Npeaks_counted < pa.rec.image_repetition + 1

        % find R peak to determine picture onset
        if toc(tic_at_last_sample) >= (sampling_period_check_Shimmer_sec + sampling_buffer_first_HB_sec)          
%                     disp('trying to find R spikes now');

            % Get time when data arrive
            time_current_packet_arrival_pre_GetSecs = GetSecs();

            % Read the latest data from shimmer data buffer, signalFormatArray defines the format of the data and signalUnitArray the unit
            [newData, signalNameArray, signalFormatArray, signalUnitArray] = shimmer.getdata('c');

            time_current_packet_arrival_post_GetSecs = GetSecs();

            % Updated 2019/04/15
            NnewData = size(newData,1);

            % Process and display data
            if (NnewData >= 2)

                if( ~FIRST_READ_HAS_OCCURRED  )
%                     disp('first reading!');
                    k_ECG_Signal_HBD = find(strcmp(string_ECG_signal_HBD, signalNameArray));

                    k_time = find(strcmp('Time Stamp', signalNameArray));
                    %k_ECG_LL_RA = find(strcmp('ECG LL-RA', signalNameArray));
                    %k_ECG_Signal_HBD = find(strcmp('ECG LA-RA', signalNameArray));

                    FIRST_READ_HAS_OCCURRED = true;

                    % Reset this to zero now that it has been used
                    sampling_buffer_first_HB_sec = 0;
                end

                X_time_current_packet_sec               = newData(:, k_time) / 1e3;
                X_time_total_sec(end+1 : end+NnewData)  = newData(:, k_time) / 1e3;                                
                Y_ECG_total(end+1 : end+NnewData)         = newData(:, k_ECG_Signal_HBD);

                Y_ECG_total_adj = ECG_adjust_baseline_spline( X_time_total_sec, Y_ECG_total, Minimum_RR_Interval_sec, Minimum_R_Prominence_mV );

%                 % Real-time plot of trial 
%                 hold on;
%                 plot(X_time_total_sec, Y_ECG_total_adj, '-k');
%                 drawnow();

                %------------------------------------------------------
                % Detection of all R spikes for this trial
                % so far (between 0 and Ntones)

                % Ensure there are enough data for peak detection                            
                time_of_total_recording_sec = X_time_total_sec(end) - X_time_total_sec(1);

                % If there is at least one peak already
                if( length(X_time_detected_peak_sec_array) >= 1 )
                    time_since_last_Rspike_sec = X_time_total_sec(end) - X_time_detected_peak_sec_array(end);
                else
                    % No actual peaks yet so get total duration of recording
                    time_since_last_Rspike_sec = X_time_total_sec(end) - X_time_total_sec(1);
                end
%                         fprintf('\nTime since last R spike = %0.4f sec \n', time_since_last_Rspike_sec);

                if( time_of_total_recording_sec > Minimum_RR_Interval_sec )
                    % Peak detection
                    [peak_Y_array, peak_X_array, peak_width_array, peak_prom_array] = ...
                        findpeaks(Y_ECG_total_adj, X_time_total_sec, ...
                        'MinPeakDistance', Minimum_RR_Interval_sec, ...
                        'MinPeakProminence', Minimum_R_Prominence_mV);

%                     % Update plot
%                     hold on;
%                     plot(peak_X_array, peak_Y_array, 'or', 'LineWidth', 3);
%                     drawnow()

                    Npeaks_found = length(peak_X_array);
                    
                    % Save time of beep
                    all_peaks_X(end+1:end+Npeaks_found) = peak_X_array;
                    all_peaks_Y(end+1:end+Npeaks_found) = peak_Y_array;
                    
                    % First peak disvorey has to be ONE peak (not more than that). If there
                    % are 2 peaks on first search, then the first one doesn't count
                    if( ~FOUND_PEAKS_ON_THIS_TRIAL && Npeaks_found > 0 )                                        
                        Npeaks_found_that_do_not_count = Npeaks_found - 1;
                        FOUND_PEAKS_ON_THIS_TRIAL = true;
                    end
                
                    %------------------------------------------------------
                    % If (1) it has been enough time since last R spike and (2) there is a new peak that hasn't previously been counted
                    %if( (toc(tic_at_last_Rspike) >= Minimum_RR_Interval_sec) && Npeaks_found > Npeaks_counted )
                    if( Npeaks_found - Npeaks_found_that_do_not_count > Npeaks_counted )

                        Npeaks_counted = Npeaks_counted + 1;

                        if Nbeeps_delivered < pa.rec.image_repetition
                            
                            % Time from packet start to R spike
                            time_from_packet_start_to_peak_sec = peak_X_array(end) - X_time_current_packet_sec(1);
                            time_beep_intended_sec = time_current_packet_arrival_pre_GetSecs + time_from_packet_start_to_peak_sec;

                            % If it LOOKS like we are waiting TOO LONG to beep, then just set it to the
                            % max waiting duration (the delay itself)
                            time_to_beep_sec = time_beep_intended_sec - GetSecs();
                            if( time_to_beep_sec > cardiac_timing_this_trial )
                                time_beep_intended_sec = GetSecs() + cardiac_timing_this_trial;
                            end

                            % image_onset = R spike + cardiac timing offset
                            image_absolute_onset_time = time_beep_intended_sec + cardiac_timing_this_trial;

                            % image offset = onset time + image duration
                            image_absolute_offset_time = image_absolute_onset_time + pa.rec.image_duration;

%                                  %minor fading in effect
%                                 for thisContrast = [20,50,80]
%                                     Screen('DrawTexture',wPtr,pointer_to_offscreen_window_with_image_on_it, [], [], 0, [], thisContrast);
%                                     DrawFormattedText(wPtr, fixation, 'center', 'center');
%                                     Screen('Flip',wPtr);
%                                 end

                            %draw image and fixation
                            Screen('DrawTexture',wPtr,pointer_to_offscreen_window_with_image_on_it);
                            Screen('TextSize', wPtr, 80);
                            DrawFormattedText(wPtr, fixation, 'center', 'center');

                            %onset
                            [~, stimOnset] = Screen('Flip',wPtr,image_absolute_onset_time);
                            seq.rec.stimOnset(trial) = stimOnset;
                            seq.rec.real_trial_onset_relative_to_block_start_time(trial) = stimOnset-seq.rec.block_start_time;
                            
                            duration_beep_late = stimOnset - time_beep_intended_sec;
                            
                            Nbeeps_delivered = Nbeeps_delivered +1;

                            % Save time of beep
                            X_time_beep_sec_array(end+1)            = stimOnset;
                            X_time_beep_late_sec_array(end+1)       = duration_beep_late;   
                            X_time_detected_peak_sec_array(end+1)   = peak_X_array(end);                                        
                            Y_ECG_detected_peak_mV_array(end+1)     = peak_Y_array(end);                                        
                            Y_peak_prominance_array(end+1)          = peak_prom_array(end);

%                             % Update plot
%                             hold on;
%                             plot(X_time_detected_peak_sec_array + X_time_beep_late_sec_array, peak_Y_array(end), 'sb', 'LineWidth', 3);
%                             % Update the plot
%                             drawnow()

                            %offset
                            [~, stimOffset] = Screen('Flip',wPtr, image_absolute_offset_time);
                            Screen('Flip',wPtr);

                        else
                            %draw + during in block ITI time after a image is repeated 4 times
                            Screen('TextSize', wPtr, 80);
                            DrawFormattedText(wPtr, fixation, 'center', 'center');
                            Screen('Flip',wPtr);

%                             % Update plot
%                             hold on;
%                             plot(peak_X_array(end)+0.05, peak_Y_array(end), 'sc', 'LineWidth', 3);
%                             % Update the plot
%                             drawnow()

                        end

                        %always show fixation after image is off screen
                        Screen('TextSize', wPtr, 80);
                        DrawFormattedText(wPtr, fixation, 'center', 'center');
                        Screen('Flip',wPtr);

                    end   

                    % Write data to file
                    time_ECG_sec = newData(:, k_time) / 1e3;
                    Y_ECG_mV = newData(:, k_ECG_Signal_HBD);
                    out = [time_ECG_sec, Y_ECG_mV];
                    % write data to end of file
                    dlmwrite(sprintf('%sencoding-ECG.csv',output_dir),out,'delimiter',',','-append');

                    out = [X_time_total_sec, Y_ECG_total_adj];
                    % write data to end of file
                    dlmwrite(sprintf('%sencoding-ECG-Baseline_Adjusted.csv',output_dir),out,'delimiter',',','-append');
                    
                end 
                
            end
            tic_at_last_sample = tic;
        end
    end
    
    %wait for response
    keyIsDown = 0;
    indexPressed = 0;
    while keyIsDown == 0 && indexPressed == 0
        [keyIsDown, keyTime, keyCodes] = KbCheck(-1);

        if mod(rand_key,2) == 0 %(M)old (N)new
            instructions = 'If you think the image is old, press the (M) button. \nIf you think the image is new, press the (N) button. ';
        else %(M)new (N)old
            instructions = 'If you think the image is old, press the (N) button. \nIf you think the image is new, press the (M) button. ';
        end
%         instruction = 'Please respond now';

        Screen('TextSize', wPtr, 50);
        DrawFormattedText(wPtr, instructions, 'center', 'center',[255,255,255]);
        Screen('Flip',wPtr);
        if keyIsDown == 1
            indexPressed = 1;
            keyIsDown = 0;
            
            disp(KbName(find(keyCodes)));
            rsp.rec.response(trial) = 1;
            rsp.rec.RT(trial)       = keyTime - seq.rec.trial_start_time(trial);
            key_pressed                = find(keyCodes);
            rsp.rec.keyName(trial)  = KbName(key_pressed(1));
            if find(KbName(key_pressed(1))==key_press_button)==find(novelty_marker_this_trial==key_press_label)
                rsp.rec.accuracy(trial) = 1;
            else
                rsp.rec.accuracy(trial) = 0;
            end 
            
        end           
    end
    Screen('Flip',wPtr);

    fixation = '+';
    Screen('TextSize', wPtr, 80);
    DrawFormattedText(wPtr, fixation, 'center', 'center');
    Screen('Flip',wPtr);

%     close(hfig);
    hfig2 = figure;
    plot(X_time_total_sec, Y_ECG_total_adj, '-k');
    hold on;
    plot(all_peaks_X, all_peaks_Y, 'or', 'LineWidth', 3);
    hold on;
    plot(X_time_detected_peak_sec_array + X_time_beep_late_sec_array, repmat(peak_Y_array(end),1,length(X_time_detected_peak_sec_array)), 'sb', 'LineWidth', 3);   
    title(sprintf('Subj %s, Trial %d, Delay=%0.4f s', pa.subjid, trial, cardiac_timing_this_trial));
    xlabel('Time (sec)');
    ylabel('Adjusted ECG (mv)');
    grid('on');
    print(hfig2, sprintf('%sECG-Trial_%02d.png', output_dir, trial), '-dpng');   
    hgsave(hfig2, sprintf('%sECG-Trial_%02d.fig', output_dir, trial));
    close(hfig2);
    
end
Screen('Flip',wPtr);
 
if strcmpi(pa.shimmer, 'true')
    shimmer.stop;
    shimmer.disconnect;
    fclose(FILE_OUT_ECG);
    fclose(FILE_OUT_ECG_Adj);
end 
 
%save csv
trial_num = [1:pa.rec.num_total_trials]';
cardiac_cycle = seq.rec.cardiac_cycle_sequence';
cardiac_timing = seq.rec.cardiac_timing_sequence';
trial_label = seq.rec.trial_label';
trial_novelty_marker = seq.rec.trial_novelty_marker';
trial_congruency_marker = seq.rec.trial_congruency_marker';
response = rsp.rec.response'; 
RT = rsp.rec.RT'; 
key = rsp.rec.keyName'; 
trial_accuracy = rsp.rec.accuracy';
stim_onset = seq.rec.stimOnset';
trial_start_time = seq.rec.trial_start_time';
real_trial_onset_relative_to_block_start_time = seq.rec.real_trial_onset_relative_to_block_start_time';
T = table(trial_num,cardiac_cycle,cardiac_timing,trial_label,trial_novelty_marker,trial_congruency_marker,RT,trial_accuracy,key,response,stim_onset,trial_start_time,real_trial_onset_relative_to_block_start_time);
writetable(T,sprintf('%ssubj-%s_Recognition_Response_by_Trial.csv', output_dir, pa.subjid));

disp('*********************************************');
disp(['End of the Experiment '])
disp('*********************************************');

% end instuctions
instructions = 'You have completed the experiment.\n Thank you!';  
Screen('TextSize', wPtr, 60);
DrawFormattedText(wPtr, instructions, 'center', 'center',[255,255,255]);
Screen('Flip',wPtr);
WaitSecs(3);
Screen('Flip',wPtr);

%save file
save(sprintf('%ssequence.mat',output_dir), 'seq')
save(sprintf('%sstimuli.mat',output_dir), 'stim')
save(sprintf('%sresponseq.mat',output_dir), 'rsp')

sca;

