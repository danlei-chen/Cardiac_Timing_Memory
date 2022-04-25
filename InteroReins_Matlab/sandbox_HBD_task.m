% Ian Kleckner
% Univeristy of Rochester Medical Center
%
% Testing heartbeat detection task using Shimmer ECG hardware
% Some code copied from Shimmer plotandwriteecgexample.m ( (C) Shimmer )
%
% 2016/12/22 Start coding

%% Setup
clear all;
close all;

%% Input
comPort = '18';

captureDuration = 10;

% sample rate in [Hz]
fs = 512; 

VIEW_DATA_TO_CHECK_THRESHOLDS = false;

CONTINUOUS_SIGNAL_MONITORING = true;

RUN_HBD_TASK = false;

%--------------------------------------------------------------------------
% Continuous signal monitoring

% How often to update the plot
update_plot_period_sec = 2;

% Keep this small to let the program run fast
duration_plot_view_sec = 2;

% Total duration to track signals for
DELAY_PERIOD_EXAMINE_SIGNALS = 30;

MINIMUM_RR_INTERVAL_SEC = 0.5;

MIN_R_PEAK_PROMINENCE_MV = 0.3;

%--------------------------------------------------------------------------
% Note: these constants are only relevant to this examplescript and are not used
% by the ShimmerHandle Class

% Number of samples in the plot 
NO_SAMPLES_IN_PLOT = 5000;                                                 

%DELAY_PERIOD_EXAMINE_SIGNALS = 10;

% Delay (in seconds) between data read operations
DELAY_PERIOD = 50e-3;                                                        

numSamples = 0;

%% Settings

%--------------------------------------------------------------------------
% Filtering settings
% mains frequency [Hz]
fm = 50;                                                                   

% corner frequency highpassfilter [Hz]; Shimmer recommends 0.5Hz for monitoring applications, 0.05Hz for diagnostic settings
fchp = 0.5;                                                                

% number of poles (HPF, LPF)
nPoles = 4;                                                                

% pass band ripple (%)
pbRipple = 0.5;                                                            

% enable (true) or disable (false) highpass filter
HPF = true;                                                                

% enable (true) or disable (false) lowpass filter
LPF = true;                                                                

% enable (true) or disable (false) bandstop filter
BSF = true;                                                                

% highpass filters for ExG channels
if (HPF)
    hpfexg1ch1 = FilterClass(FilterClass.HPF,fs,fchp,nPoles,pbRipple);
    hpfexg1ch2 = FilterClass(FilterClass.HPF,fs,fchp,nPoles,pbRipple);
    hpfexg2ch1 = FilterClass(FilterClass.HPF,fs,fchp,nPoles,pbRipple);
    hpfexg2ch2 = FilterClass(FilterClass.HPF,fs,fchp,nPoles,pbRipple);
end
% lowpass filters for ExG channels
if (LPF)
    lpfexg1ch1 = FilterClass(FilterClass.LPF,fs,fs/2-1,nPoles,pbRipple);
    lpfexg1ch2 = FilterClass(FilterClass.LPF,fs,fs/2-1,nPoles,pbRipple);
    lpfexg2ch1 = FilterClass(FilterClass.LPF,fs,fs/2-1,nPoles,pbRipple);
    lpfexg2ch2 = FilterClass(FilterClass.LPF,fs,fs/2-1,nPoles,pbRipple);
end
% bandstop filters for ExG channels;
% cornerfrequencies at +1Hz and -1Hz from mains frequency
if (BSF)
    bsfexg1ch1 = FilterClass(FilterClass.LPF,fs,[fm-1,fm+1],nPoles,pbRipple);
    bsfexg1ch2 = FilterClass(FilterClass.LPF,fs,[fm-1,fm+1],nPoles,pbRipple);
    bsfexg2ch1 = FilterClass(FilterClass.LPF,fs,[fm-1,fm+1],nPoles,pbRipple);
    bsfexg2ch2 = FilterClass(FilterClass.LPF,fs,[fm-1,fm+1],nPoles,pbRipple);
end


%% Connect to Shimmer device

% Connect to Shimmer at specified COM port
fprintf('\n\nAttempting to connect to Shimmer, this may take 10 sec...\n');
shimmer = ShimmerHandleClass(comPort);

% assign user friendly macros for setenabledsensors
SensorMacros = SetEnabledSensorsMacrosClass;                               

% Ensure connection
if (shimmer.connect)
    
    fprintf('\n\nSuccessfully connected to Shimmer device!\n');
    
    % Select sampling rate
    shimmer.setsamplingrate(fs);                                           
    
    % Select internal expansion board; select 'ECG' to enable both SENSOR_EXG1 and SENSOR_EXG2 
    shimmer.setinternalboard('ECG');                                       
    
    % Disable other sensors
    shimmer.disableallsensors;                                             
    
    % Enable SENSOR_EXG1 and SENSOR_EXG2
    shimmer.setenabledsensors(SensorMacros.ECG,1);

    % TRUE if the shimmer starts streaming 
    if( VIEW_DATA_TO_CHECK_THRESHOLDS && shimmer.start )
        
        % Wait desired time
        fprintf('\n\nCollecting data for %0.1f sec...', DELAY_PERIOD_EXAMINE_SIGNALS);
        pause(DELAY_PERIOD_EXAMINE_SIGNALS);
        fprintf('done.');
        fprintf('\nExamine signals in plots to find thresholds');
        
        % Read the latest data from shimmer data buffer, signalFormatArray defines the format of the data and signalUnitArray the unit
        [newData, signalNameArray, signalFormatArray, signalUnitArray] = shimmer.getdata('c');
        
        % Stop data collection
        shimmer.stop;
        
        % Display data
        k_time = find(strcmp('Time Stamp', signalNameArray));
        k_ECG_LL_RA = find(strcmp('ECG LL-RA', signalNameArray));
        k_ECG_LA_RA = find(strcmp('ECG LA-RA', signalNameArray));
        
        time_array_sec = newData(:, k_time) / 1e3;
        sampling_period_sec = time_array_sec(2) - time_array_sec(1);
        
        % Plot ECG data vs. Time
        Y1 = newData(:, k_ECG_LL_RA);
        subplot(3,1,1);
        plot(time_array_sec, Y1, '-k');
        ylabel('LL - RA');
        title(sprintf('ECG at %0.1f Hz', 1/sampling_period_sec));
        
        % Plot other form of ECG
        Y2 = newData(:, k_ECG_LA_RA);
        subplot(3,1,2)        
        plot(time_array_sec, Y2, '-b');
        ylabel('LA - RA');
        
        % Plot time derivative of ECG data vs. Time
        subplot(3,1,3);
        Ydiff = diff(Y2);
        time_array_sec_diff = time_array_sec + sampling_period_sec / 2;
        time_array_sec_diff(end) = [];
        plot(time_array_sec_diff, Ydiff, '-r');
        ylabel('d/dt( LA-RA )');
        
    end
    
    %error('pause');
    
    %% Continuous signal monitoring
    
    
    % TRUE if the shimmer starts streaming 
    if( CONTINUOUS_SIGNAL_MONITORING && shimmer.start )
        
        fprintf('\n\nStarting continuous signal monitoring for %0.1f sec', DELAY_PERIOD_EXAMINE_SIGNALS);
        
        % Wait a short time for first data
        pause(1);
        
        tic_start_time = tic;
        tic_last_plot_update = tic;
        
        time_array_sec = [];
        Y1 = [];
        Y2 = [];
        peak_prom_array_all = [];
                
        while( toc(tic_start_time) < DELAY_PERIOD_EXAMINE_SIGNALS )
            
            if( toc(tic_last_plot_update) > update_plot_period_sec )
                
                % Read the latest data from shimmer data buffer, signalFormatArray defines the format of the data and signalUnitArray the unit
                [newData, signalNameArray, signalFormatArray, signalUnitArray] = shimmer.getdata('c');
                if( ~isempty(newData) )
                    
                    NnewData = length(newData);
                    
                    % Display data
                    k_time = find(strcmp('Time Stamp', signalNameArray));
                    k_ECG_LL_RA = find(strcmp('ECG LL-RA', signalNameArray));
                    k_ECG_LA_RA = find(strcmp('ECG LA-RA', signalNameArray));

                    time_array_sec(end+1 : end+NnewData) = newData(:, k_time) / 1e3;
                    sampling_period_sec = time_array_sec(2) - time_array_sec(1);

                    %------------------------------------------------------
                    % Crop plot to desired duration
                    Nplot_data = duration_plot_view_sec / sampling_period_sec;

                    if( length(time_array_sec) > Nplot_data )
                        % Too much data, crop to most recent portion
                        k_plot_data = length(time_array_sec)-Nplot_data : length(time_array_sec);
                    else
                        % Not enough data yet, so display it all
                        k_plot_data = 1:length(time_array_sec);
                    end

                    time_array_sec = time_array_sec(k_plot_data);

                    %------------------------------------------------------
                    % Plot ECG data vs. Time
                    Y1(end+1 : end+NnewData) = newData(:, k_ECG_LL_RA);
                    Y1 = Y1(k_plot_data);
                    
                    subplot(3,1,1);
                    plot(time_array_sec, Y1, '-k');
                    grid('on');
                    ylabel('LL - RA');
                    title(sprintf('ECG at %0.1f Hz (%0.1f / %0.1f sec)', 1/sampling_period_sec, toc(tic_start_time), DELAY_PERIOD_EXAMINE_SIGNALS));

                    %------------------------------------------------------
                    % Plot other form of ECG
                    Y2(end+1 : end+NnewData) = newData(:, k_ECG_LA_RA);
                    Y2 = Y2(k_plot_data);
                    
                    subplot(3,1,2)        
                    plot(time_array_sec, Y2, '-b');
                    grid('on');
                    ylabel('LA - RA');

                    %------------------------------------------------------
                    % Detection of R spikes
                    [peak_Y_array, peak_X_array, peak_width_array, peak_prom_array] = ...
                        findpeaks(Y2, time_array_sec, ...
                        'MinPeakDistance', MINIMUM_RR_INTERVAL_SEC);
                    
                    k_valid_peaks = peak_prom_array > MIN_R_PEAK_PROMINENCE_MV;
                    
                    peak_Y_array = peak_Y_array(k_valid_peaks);
                    peak_X_array = peak_X_array(k_valid_peaks);
                    peak_prom_array = peak_prom_array(k_valid_peaks);
                    
                    peak_prom_array_all(end+1 : end+length(peak_prom_array)) = peak_prom_array;

                    hold('all');
                    plot(peak_X_array, peak_Y_array, 'or', 'LineWidth', 3);
                    
                    YLIM = get(gca, 'YLim');
                    
                    for p = 1:length(peak_prom_array)
                        text(peak_X_array(p), YLIM(1), ...
                            sprintf('P=%0.2f', peak_prom_array(p)), ...
                            'HorizontalAlignment', 'center', ...
                            'VerticalAlignment', 'bottom');
                    end
                    hold('off');

                    %------------------------------------------------------
                    % Plot time derivative of ECG data vs. Time
                    subplot(3,1,3);
                    Ydiff = diff(Y2);
                    time_array_sec_diff = time_array_sec + sampling_period_sec / 2;
                    time_array_sec_diff(end) = [];
                    plot(time_array_sec_diff, Ydiff, '-r');
                    grid('on');
                    ylabel('d/dt( LA-RA )');
                    
                    %------------------------------------------------------
                    % Update the plot
                    drawnow()

                    % Restart timer
                    tic_last_plot_udpate = tic;
                end
            end
        end  
        
        % Report on peak prominances
        figure;
        plot(sort(peak_prom_array_all), 'ok');
        xlabel('Peak Number');
        ylabel('Peak Prominence (mV) from LA-RA');
        title(sprintf('Prominences of %d peaks range from %0.2f - %0.2f mV', ...
            length(peak_prom_array_all), min(peak_prom_array_all), max(peak_prom_array_all)));
        
    end
    
    shimmer.stop;
    
    
    error('end');
    
    %% Start up streaming again for HBD task
    
    %----------------------------------------------------------------------
    % Inputs
    
    string_ECG_signal_HBD = 'ECG LA-RA';
    %k_ECG_HBD = k_ECG_LA_RA;
    Y_Rspike_ECG_HBD = 2.3;
    
    MINIMUM_RR_INTERVAL_SEC = 0.5;

    time_limit = 5;
    display_sampling_period_sec = 0.025;
    
    delay_Rspike_beep_sec = 0.5;
    
    
    %----------------------------------------------------------------------
    if( RUN_HBD_TASK && shimmer.start )
        fprintf('\n\n**** Running HBD task');
        pause(1);
        
        % Initialize timers and other variables
        tic_at_last_sample = tic;
        tic_at_last_Rspike = tic;
        tic_at_experiment_start = tic;
        ALREADY_BEEPED_FOR_THIS_RSPIKE = true;
        FIRST_READ_HAS_OCCURRED = false;
        
        % Keep going until experiment is done
        while( toc(tic_at_experiment_start) < time_limit )            
            
            % Check whether it is time to sample the Shimmer again        
            if( toc(tic_at_last_sample) >= display_sampling_period_sec )
                
                % Read the latest data from shimmer data buffer, signalFormatArray defines the format of the data and signalUnitArray the unit
                [newData, signalNameArray, signalFormatArray, signalUnitArray] = shimmer.getdata('c');
                
                if( ~isempty(newData) )
                    
                    if( ~FIRST_READ_HAS_OCCURRED  )
                        k_ECG_Signal_HBD = find(strcmp(string_ECG_signal_HBD, signalNameArray));
                        FIRST_READ_HAS_OCCURRED = true;
                    end

                    Y_ECG_current = newData(:, k_ECG_Signal_HBD);

                    % Check for new Rspike based on (1) whether it has been
                    % long enough and (2) whether signal is over threshold
                    if( (toc(tic_at_last_Rspike) >= MINIMUM_RR_INTERVAL_SEC) && any(Y_ECG_current > Y_Rspike_ECG_HBD) )
                        tic_at_last_Rspike = tic;
                        fprintf('\nRspike');                    
                        ALREADY_BEEPED_FOR_THIS_RSPIKE = false;
                    end

                    % Print results from sampling new data
                    %[Nrows, Ncols] = size(newData);
                    %fprintf('\n%0.3f sec\t%d\t%d', time_current, Nrows, Ncols);

                    tic_at_last_sample = tic;
                end
            end
            
            % Present beep at delay from Rspike            
            if( ~ALREADY_BEEPED_FOR_THIS_RSPIKE && toc(tic_at_last_Rspike) > delay_Rspike_beep_sec )
                fprintf('\nBEEP');
                ALREADY_BEEPED_FOR_THIS_RSPIKE = true;
            end            
        end
        
        fprintf('\n\nTime limit is up: %0.1f sec', time_limit);
    end
    
    shimmer.stop;
    
    
    %%
    shimmer.disconnect;
    
    %fprintf('The percentage of received packets: %d \n', shimmer.getpercentageofpacketsreceived(timeStamp)); % Detect loss packets

end

fprintf('\n\nAll done!');