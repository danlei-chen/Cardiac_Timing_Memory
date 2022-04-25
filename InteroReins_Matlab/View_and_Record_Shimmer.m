% Ian Kleckner
% Univeristy of Rochester Medical Center
%
% Testing heartbeat detection task using Shimmer ECG hardware
% Some code copied from Shimmer plotandwriteecgexample.m ( (C) Shimmer )
%
% 2016/12/22 Start coding
% 2017/05/23 Only use one plot for ECG (LA-RA)
% 2019/04/15 Update and fix error on assigning new data

function View_and_Record_Shimmer( experiment_specs )

    %% Input
    
    %----------------------------------------------------------------------
    % ECG inputs
    
    % Whether to use ECG or not
    USE_ECG = experiment_specs.Use_ECG;
    
    % Read in COM port and format it as a string
    comPort_ECG = sprintf('%d', experiment_specs.COM_ECG);
    
    % sample rate in [Hz]
    sampling_rate_ECG_Hz = experiment_specs.Sampling_rate_ECG_Hz; 
    
    %----------------------------------------------------------------------
    % R spike detection settings
    Minimum_RR_Interval_sec = experiment_specs.Minimum_RR_Interval_sec;
    Minimum_R_Prominence_mV = experiment_specs.Minimum_R_Prominence_mV;
    
    
    %----------------------------------------------------------------------
    % Baseline correction [2019/05/03]
    ISOELECTRIC_TIME_BEFORE_RSPIKE_SEC = 66e-3;
    
    %----------------------------------------------------------------------
    % EDA inputs
    
    % Whether to use EDA or not
    USE_EDA = experiment_specs.Use_EDA;
    
    % Read in COM port and format it as a string
    comPort_EDA = sprintf('%d', experiment_specs.COM_EDA);
    
    % sample rate in [Hz]
    sampling_rate_EDA_Hz = experiment_specs.Sampling_rate_EDA_Hz; 
    
    %----------------------------------------------------------------------
    % Display settings
    
    % How much initial data to ignore
    %duration_first_samples_ignore_sec = 0.5;
    
    % How often to update the plot
    Plot_update_period_sec = experiment_specs.Plot_update_period_sec;

    % Keep this small to let the program run fast
    Plot_viewable_duration_sec = experiment_specs.Plot_viewable_duration_sec;

    % Total duration to track signals for
    Test_duration_sec = experiment_specs.Test_duration_sec;
    
    % Write text file with Shimmer data
    WRITE_DATA_TO_FILE = experiment_specs.Write_data_to_file;
    
    % Extra name of this event
    Event_name = experiment_specs.Event_name;
    
    
    OUTPUT_DEBUG = 0;
    
    %% Setup
    
    KbName('UnifyKeyNames')
    key_escape      = KbName('ESCAPE');
    key_left_shift   = KbName('LeftShift');
    key_right_shift  = KbName('RightShift');
    
    %% Connect to Shimmer devices
    
    % assign user friendly macros for setenabledsensors
    SensorMacros = SetEnabledSensorsMacrosClass;

    if( USE_ECG )
        % Connect to Shimmer at specified COM port
        fprintf('\n\n----------------------------------------------------');
        fprintf('\nAttempting to connect to ECG Shimmer, this may take 10 sec...\n');
        shimmer_ECG = ShimmerHandleClass(comPort_ECG);
        
        if (shimmer_ECG.connect)
            fprintf('\nSuccessfully connected to ECG Shimmer!\n');
            
            % Select sampling rate
            shimmer_ECG.setsamplingrate(sampling_rate_ECG_Hz);                                           

            % Select internal expansion board; select 'ECG' to enable both SENSOR_EXG1 and SENSOR_EXG2 
            shimmer_ECG.setinternalboard('ECG');                                       

            % Disable other sensors
            shimmer_ECG.disableallsensors;                                             

            % Enable SENSOR_EXG1 and SENSOR_EXG2
            shimmer_ECG.setenabledsensors(SensorMacros.ECG,1);
            
            shimmer_ECG.setrealtimeclock();
            
        else
            error('Cannot connect to ECG Shimmer');
        end
    end
    
    if( USE_EDA )
        % Connect to Shimmer at specified COM port
        fprintf('\n\n----------------------------------------------------');
        fprintf('\nAttempting to connect to EDA Shimmer, this may take 10 sec...\n');
        shimmer_EDA = ShimmerHandleClass(comPort_EDA);
        
        if (shimmer_EDA.connect)
            fprintf('\nSuccessfully connected to Shimmer device!\n');
            
            % Select sampling rate
            shimmer_EDA.setsamplingrate(sampling_rate_EDA_Hz);                                           

            % Select internal expansion board; select 'ECG' to enable both SENSOR_EXG1 and SENSOR_EXG2 
            shimmer_EDA.setinternalboard('GSR');                                       

            % Disable other sensors
            shimmer_EDA.disableallsensors;                                             

            % Enable SENSOR_EXG1 and SENSOR_EXG2
            shimmer_EDA.setenabledsensors(SensorMacros.GSR,1);
            
            shimmer_EDA.setrealtimeclock();
        else
            error('Cannot connect to EDA Shimmer');
        end
    end
    
    %% Prepare to start streaming data
    
    userInput = questdlg({'Click START to start streaming data', '', 'To STOP press Shift + Esc'}, 'Start streaming', ...
        'START', 'Abort Program', 'Start?');

    if( strcmp(userInput, 'Abort Program') )
        error('User aborted program instead of starting program');
    end
    
    %------------------------------------------------------------------
    % Make output folders, etc.
    % Create an output directory to save results to files
    programName = 'V+R';
    %PPid_nospaces = strrep(experiment_specs.PPid, ' ', '_');
    PPid = sprintf('s%03d', experiment_specs.PP_Number);
    output_dir = sprintf('out-%s-%s_%s-%s-T%d',datestr(now, 'yyyy.mm.dd-HHMM'), programName, Event_name, PPid, experiment_specs.Timepoint);
    mkdir(output_dir);
    
    diary(sprintf('%s/diary.txt', output_dir));
    
    OUT_SETTINGS = fopen(sprintf('%s/SETTINGS.txt', output_dir), 'w');
    writeSpecificationsHeader( OUT_SETTINGS, experiment_specs )
    fclose(OUT_SETTINGS);
    
    %------------------------------------------------------------------
    % Set up for plots
    if( USE_ECG & USE_EDA )
        Nsubplots = 2;
        k_subplot_ECG1 = [];
        k_subplot_ECG2 = 1;
        k_subplot_EDA = 2;

    elseif( USE_ECG & ~USE_EDA )
        Nsubplots = 1;
        k_subplot_ECG1 = [];
        k_subplot_ECG2 = 1;
        k_subplot_EDA = [];

    elseif( ~USE_ECG & USE_EDA )
        Nsubplots = 1;
        k_subplot_ECG1 = [];
        k_subplot_ECG2 = [];
        k_subplot_EDA = 1;

    else
        error('Invalid combination of requested Shimmer data');
    end

    %------------------------------------------------------------------
    % TRUE if the shimmers start streaming
    if( USE_ECG )
        if( shimmer_ECG.start )
            ECG_IS_OK = true;
        else
            ECG_IS_OK = false;
        end
    else
        ECG_IS_OK = true;
    end
    
    if( USE_EDA )
        if( shimmer_EDA.start )
            EDA_IS_OK = true;
        else
            EDA_IS_OK = false;
        end   
    else
        EDA_IS_OK = true;
    end
    
    %----------------------------------------------------------------------
    % Write output files
    if( WRITE_DATA_TO_FILE )
        if( USE_ECG )
            FILE_OUT_ECG = fopen(sprintf('%s/out-data-ECG.txt', output_dir), 'w');
            fprintf(FILE_OUT_ECG, 'Time(sec)\tECG_Raw(mV)');
        end
        
        if( USE_EDA )
            FILE_OUT_EDA = fopen(sprintf('%s/out-data-EDA.txt', output_dir), 'w');    
            fprintf(FILE_OUT_EDA, 'Time(sec)\tEDA(uS)V');
        end
    end
    
    if( ECG_IS_OK && EDA_IS_OK )

        fprintf('\n\nStarting continuous signal monitoring for %0.1f sec', Test_duration_sec);
        fprintf('\n\n** To stop recording early hold Shift + Esc for a couple seconds **');
        
        %------------------------------------------------------------------
        % Initialize
        tic_start_time = tic;
        tic_last_plot_update = tic;

        time_ECG_array_sec_plot = [];
        time_EDA_array_sec_plot = [];
        Y1_ECG_plot = [];
        Y2_ECG_plot = [];
        Y1_EDA_plot = [];
        
        time_ECG_array_sec_all = [];
        time_EDA_array_sec_all = [];
        Y1_ECG_all = [];
        Y2_ECG_all = [];
        Y1_EDA_all = [];
        
        peak_prom_array_all = [];
        HAVENT_READ_ECG_YET = true;
        HAVENT_READ_EDA_YET = true;
        
        USER_EXITS_VIA_KEYBOARD = false;
        
        
        %HAVENT_PAUSED_YET = true;
        
        hfig = figure;
        
        % Wait a short time for first batch of data
        pause(3);
        
        %------------------------------------------------------------------
        while( ~USER_EXITS_VIA_KEYBOARD && toc(tic_start_time) < Test_duration_sec )
            
            % Code to test stopping and restarting data collection
            % 2019/04/30
%             if( HAVENT_PAUSED_YET && toc(tic_start_time) > 10 )
%                 
%                 fprintf('\nStopping for a moment...');
%                 shimmer_ECG.stop;
%                 WaitSecs(2)
%                 shimmer_ECG.start;
%                 fprintf('\nResuming');
%             end
            

            %--------------------------------------------------------------
            % Check for keyboard button
            
            % Check the state of the keyboard.
            [ keyIsDown, seconds, keyCode ] = KbCheck;
                        
            if( (keyCode(key_escape) && keyCode(key_left_shift)) || ...
                (keyCode(key_escape) && keyCode(key_right_shift)) )
                fprintf('\n\nUser stopped recording via keyboard.');                
                USER_EXITS_VIA_KEYBOARD = true;                
            end

            if( toc(tic_last_plot_update) > Plot_update_period_sec )

                % Read the latest data from shimmer data buffer, signalFormatArray defines the format of the data and signalUnitArray the unit
                if( USE_ECG )                    
                    [newData_ECG, signalNameArray_ECG, signalFormatArray, signalUnitArray] = shimmer_ECG.getdata('c');
                end
                
                if( USE_EDA )
                    [newData_EDA, signalNameArray_EDA, signalFormatArray, signalUnitArray] = shimmer_EDA.getdata('c');
                end
                
                % Updated 2019/04/12
                NnewData_ECG = size(newData_ECG,1);
                
                % Process and display data
                if( USE_ECG && NnewData_ECG >= 2 )
                    
                    % Display data
                    k_time = find(strcmp('Time Stamp', signalNameArray_ECG));
                    k_ECG_LL_RA = find(strcmp('ECG LL-RA', signalNameArray_ECG));
                    k_ECG_LA_RA = find(strcmp('ECG LA-RA', signalNameArray_ECG));
                    
                    if( OUTPUT_DEBUG )
                        signalNameArray_ECG
                        k_time
                        k_ECG_LL_RA
                        k_ECG_LA_RA

                        'new data'
                        size(newData_ECG)
                        length(newData_ECG)
                        NnewData_ECG
                        
                        '\n\nexisting data'
                        size(time_ECG_array_sec_plot)

                        newData_ECG
                    end
                    
                    time_ECG_array_sec_all(end+1 : end+NnewData_ECG) = newData_ECG(:, k_time) / 1e3;
                    
                    if( HAVENT_READ_ECG_YET )
                        sampling_period_ECG_sec = time_ECG_array_sec_all(2) - time_ECG_array_sec_all(1);
                        HAVENT_READ_ECG_YET = false;
                    end
                    
                    ISOELECTRIC_SAMPLES_BEFORE_RSPIKE = floor(ISOELECTRIC_TIME_BEFORE_RSPIKE_SEC / sampling_period_ECG_sec);

                    %------------------------------------------------------
                    % Crop plot to desired duration
                    Nplot_data = Plot_viewable_duration_sec / sampling_period_ECG_sec;
                    
                    %Nfirst_samples_ECG_ignore = duration_first_samples_ignore_sec / sampling_period_ECG_sec;

                    if( length(time_ECG_array_sec_all) > Nplot_data )
                        % Too much data, crop to most recent portion
                        %  Fixed bug by adding floor() on 2019/05/02
                        k_plot_data = floor(length(time_ECG_array_sec_all)-Nplot_data) : floor(length(time_ECG_array_sec_all));
                    else
                        % Not too much data yet, so display it all
                        k_plot_data = 1:length(time_ECG_array_sec_all);                        
                    end
                    
                    try
                        time_ECG_array_sec_plot = time_ECG_array_sec_all(k_plot_data);
                    catch err
                        k_plot_data                        
                        err                        
                        error('Busted');
                    end
                        
                    %------------------------------------------------------
                    % Plot ECG data vs. Time
                    %{
                    Y1_ECG(end+1 : end+NnewData_ECG) = newData_ECG(:, k_ECG_LL_RA);
                    Y1_ECG = Y1_ECG(k_plot_data);

                    subplot(Nsubplots, 1, k_subplot_ECG1, 'Parent', hfig);
                    plot(time_ECG_array_sec, Y1_ECG, '-k');
                    grid('on');
                    ylabel('ECG LL - RA');
                    title(sprintf('Prominence Values for ECG at %0.1f Hz (%0.1f / %0.1f sec)', 1/sampling_period_ECG_sec, toc(tic_start_time), Test_duration_sec));
                    %}

                    %------------------------------------------------------
                    % Plot other form of ECG
                    Y2_ECG_all(end+1 : end+NnewData_ECG) = newData_ECG(:, k_ECG_LA_RA);
                    Y2_ECG_plot = Y2_ECG_all(k_plot_data);
                    
                    Y2_ECG_plot_adj = ECG_adjust_baseline_spline( time_ECG_array_sec_plot, Y2_ECG_plot, Minimum_RR_Interval_sec, Minimum_R_Prominence_mV );

                    subplot(Nsubplots, 1, k_subplot_ECG2, 'Parent', hfig);        
                    plot(time_ECG_array_sec_plot, Y2_ECG_plot_adj, '-b');
                    grid('on');
                    ylabel('Adjusted ECG LA-RA (mV)');
                    title(sprintf('ECG at %0.1f Hz (%0.1f / %0.1f sec)', 1/sampling_period_ECG_sec, toc(tic_start_time), Test_duration_sec));

                    %------------------------------------------------------

                    if( Minimum_RR_Interval_sec < ( time_ECG_array_sec_plot(end)-time_ECG_array_sec_plot(1) ) )
                        % Detection of R spikes
                        [peak_Y_array, peak_X_array, peak_width_array, peak_prom_array] = ...
                            findpeaks(Y2_ECG_plot_adj, time_ECG_array_sec_plot, ...
                            'MinPeakDistance', Minimum_RR_Interval_sec, ...
                            'MinPeakProminence', Minimum_R_Prominence_mV);
                        
                        peak_prom_array_all(end+1 : end+length(peak_prom_array)) = peak_prom_array;
                        
                        hold('all');
                        plot(peak_X_array, peak_Y_array, 'or', 'LineWidth', 3);
                        
                        try
                            HR_average_bpm = length(peak_X_array) / range(time_ECG_array_sec_plot) * 60;
                        catch err
                            HR_average_bpm = NaN;
                        end
                        title(sprintf('%s: Avg HR %0.0f bpm, ECG at %0.1f Hz (%0.0f/%0.0f sec). Shift+Esc to stop', ...
                            PPid, HR_average_bpm, 1/sampling_period_ECG_sec, toc(tic_start_time), Test_duration_sec));
                        
                        YLIM = get(gca, 'YLim');

                        for p = 1:length(peak_prom_array)
                            text(peak_X_array(p), YLIM(1), ...
                                sprintf('%0.2f', peak_prom_array(p)), ...
                                'HorizontalAlignment', 'center', ...
                                'VerticalAlignment', 'bottom');
                        end
                        hold('off');


                        %------------------------------------------------------
                        % Plot time derivative of ECG data vs. Time
                        if( 0 )                            
                            subplot(3,1,3, 'Parent', hfig);
                            Ydiff = diff(Y2);
                            time_array_sec_diff = time_array_sec + sampling_period_sec / 2;
                            time_array_sec_diff(end) = [];
                            plot(time_array_sec_diff, Ydiff, '-r');
                            grid('on');
                            ylabel('d/dt( LA-RA )');
                        end                        
                    end
                    
                    
                    %------------------------------------------------------
                    % Write data to file

                    if( WRITE_DATA_TO_FILE )
                        time_ECG_sec = newData_ECG(:, k_time) / 1e3;
                        Y_ECG_mV = newData_ECG(:, k_ECG_LA_RA);
                        out = [time_ECG_sec Y_ECG_mV]';
                        
                        fprintf(FILE_OUT_ECG, '\n%f\t%f', out(:));
                    end
                end
                
                % Process and display data
                if( USE_EDA && ~isempty(newData_EDA) )

                    NnewData_EDA = length(newData_EDA);
                    
                    % Display data
                    k_time = find(strcmp('Time Stamp', signalNameArray_EDA));
                    k_GSR = find(strcmp('GSR', signalNameArray_EDA));

                    time_EDA_array_sec_all(end+1 : end+NnewData_EDA) = newData_EDA(:, k_time) / 1e3;
                    
                    if( HAVENT_READ_EDA_YET )
                        sampling_period_EDA = time_EDA_array_sec_all(2) - time_EDA_array_sec_all(1);
                        HAVENT_READ_EDA_YET = false;
                    end

                    %------------------------------------------------------
                    % Crop plot to desired duration
                    Nplot_data = Plot_viewable_duration_sec / sampling_period_EDA;

                    if( length(time_EDA_array_sec_all) > Nplot_data )
                        % Too much data, crop to most recent portion
                        k_plot_data = length(time_EDA_array_sec_all)-Nplot_data : length(time_EDA_array_sec_all);
                    else
                        % Not enough data yet, so display it all
                        k_plot_data = 1:length(time_EDA_array_sec_all);
                    end

                    time_EDA_array_sec_plot = time_EDA_array_sec_all(k_plot_data);

                    %------------------------------------------------------
                    % Plot EDA data vs. Time (converting GSR in kOhms to EDA in uS)
                    Y1_EDA_all(end+1 : end+NnewData_EDA) = 1e3 ./ newData_EDA(:, k_GSR);
                    Y1_EDA_plot = Y1_EDA_all(k_plot_data);
                    
                    subplot(Nsubplots, 1, k_subplot_EDA, 'Parent', hfig);
                    plot(time_EDA_array_sec_plot, Y1_EDA_plot, '-k');
                    grid('on');
                    ylabel('EDA (uS)');
                    xlabel('Time (sec)');
                    
                    if( ~USE_ECG )
                        title(sprintf('EDA at %0.1f Hz (%0.1f / %0.1f sec)', 1/sampling_period_EDA, toc(tic_start_time), Test_duration_sec));
                    end
                    
                    if( WRITE_DATA_TO_FILE )
                        time_EDA_sec = newData_EDA(:, k_time) / 1e3;
                        Y_EDA_uS = newData_EDA(:, k_GSR);
                        
                        out = [time_EDA_sec Y_EDA_uS]';

                        fprintf(FILE_OUT_EDA, '\n%f\t%f', out(:));
                    end
                end
                
                %------------------------------------------------------
                % Update the plot
                drawnow()
                
                % Restart timer
                tic_last_plot_udpate = tic;
            end
        end
        
        %------------------------------------------------------
        % Save figure to output
        print(hfig, sprintf('%s/out-signals-some.png', output_dir), '-dpng');

        %% Show final signals and prominances
        
        hfig2 = figure;
        
        if( USE_ECG )            
            % Baseline correction 2019/05/03            
            Y2_ECG_all_adj = ECG_adjust_baseline_spline( time_ECG_array_sec_all, Y2_ECG_all, Minimum_RR_Interval_sec, Minimum_R_Prominence_mV );
            
            % Plot corrected signal
            subplot(Nsubplots, 1, 1, 'Parent', hfig2);        
            %plot(time_ECG_array_sec_all, Y2_ECG_all, '-b');
            plot(time_ECG_array_sec_all, Y2_ECG_all_adj, '-b');
            grid('on');
            ylabel('Adjusted ECG LA-RA (mv)');
            
            if( ~USE_EDA )
                xlabel('Time (sec)');
            end
            
            %------------------------------------------------------
            % Detection of R spikes from corrected signal
            [peak_Y_array, peak_X_array, peak_width_array, peak_prom_array] = ...
                findpeaks(Y2_ECG_all_adj, time_ECG_array_sec_all, ...
                'MinPeakDistance', Minimum_RR_Interval_sec, ...
                'MinPeakProminence', Minimum_R_Prominence_mV );

            hold('all');
            plot(peak_X_array, peak_Y_array, 'or', 'LineWidth', 3);
            
            title(sprintf('%s: Avg HR = %0.1f bpm, ECG at %0.1f Hz (%0.1f / %0.1f sec)', ...
                            PPid, HR_average_bpm, 1/sampling_period_ECG_sec, toc(tic_start_time), Test_duration_sec));
            
            
            if( WRITE_DATA_TO_FILE )
                FILE_OUT_ECG_ADJUSTED = fopen(sprintf('%s/out-data-ECG-Baseline_Adjusted.txt', output_dir), 'w');
                fprintf(FILE_OUT_ECG_ADJUSTED, 'Time(sec)\tECG_Adj(mV)');
                
                out = [time_ECG_array_sec_all; Y2_ECG_all_adj];
                fprintf(FILE_OUT_ECG_ADJUSTED, '\n%f\t%f', out(:));
                fclose(FILE_OUT_ECG_ADJUSTED);
                
                FILE_OUT_ECG_RSPIKES = fopen(sprintf('%s/out-data-ECG-Rspikes.txt', output_dir), 'w');
                fprintf(FILE_OUT_ECG_RSPIKES, 'Time(sec)\tRspike(mV)\tProminence(mV)\tIBI(sec)');
                
                IBI_sec_array = diff( [NaN peak_X_array] );
                fprintf(FILE_OUT_ECG_RSPIKES, '\n%f\t%f\t%f\t%f', [peak_X_array; peak_Y_array; peak_prom_array; IBI_sec_array]);
                
                fclose(FILE_OUT_ECG_RSPIKES);
                
            end
        end
        
        if( USE_EDA )
            subplot(Nsubplots, 1, k_subplot_EDA, 'Parent', hfig2);
            plot(time_EDA_array_sec_all, Y1_EDA_all, '-k');
            grid('on');
            ylabel('EDA (uS)');
            xlabel('Time (sec)');
        end
        
        print(hfig2, sprintf('%s/out-signals-all.png', output_dir), '-dpng');
        hgsave(hfig2, sprintf('%s/out-signals-all.fig', output_dir));

        %------------------------------------------------------------------
        % Report on peak prominances
        if( USE_ECG )
            hfig3 = figure;
            plot(sort(peak_prom_array), 'ok');
            xlabel('Peak Number');
            ylabel('Peak Prominence (mV) from LA-RA');
            title(sprintf('Prominences of %d peaks range from %0.2f - %0.2f mV', ...
                length(peak_prom_array), min(peak_prom_array), max(peak_prom_array)));
            
            print(hfig3, sprintf('%s/out-ECG_Prominences.png', output_dir), '-dpng');
        end

        % Stop acquiring data
        if( USE_ECG )
            shimmer_ECG.stop;
        end
        
        if( USE_EDA )
            shimmer_EDA.stop;
        end
    else
        error('Could not start streaming Shimmer. Not sure why');            
    end

    % Disconnect
    if( USE_ECG )
        shimmer_ECG.disconnect;
        if( WRITE_DATA_TO_FILE )
            fclose(FILE_OUT_ECG);
        end
    end
    
    if( USE_EDA )
        shimmer_EDA.disconnect;
        if( WRITE_DATA_TO_FILE )
            fclose(FILE_OUT_EDA);
        end
    end
    
fprintf('\n\nAll done!');
diary('off');
    
end