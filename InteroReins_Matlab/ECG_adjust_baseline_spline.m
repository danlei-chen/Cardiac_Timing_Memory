function [ Y_ECG_adjusted ] = ECG_adjust_baseline_spline( X_sec, Y_ECG, Minimum_RR_Interval_sec, Minimum_R_Prominence_ECG )
%ECG_adjust_baseline_spline Takes ECG data and corrects wandering baselines
%   Detects R spikes, finds isoelectric point by going 66 msec before each
%   R spike, then uses a cubic spline interpolation across all those
%   points. Uses method from Meyer and Keiser 1977, which is apparently
%   very popular according to Romero et al Baseline wander removal methods for
%   ECG signals: A comparative study
%
% INPUT
%   X_sec and Y_ECG are the X, Y pairs for the ECG signal
%       X_sec must be in seconds
%   Minimum_RR_Interval_sec     Minimum time between finding R spikes
%   Minimum_R_Prominence_ECG    Minimum prominence of the R spike
%
%   * Y_ECG and Minimum_R_Prominence_ECG must be in the same units (e.g.,  mV)
%
% CHANGELOG
%   2019/05/03 Start coding
%
% CREDITS
%   Ian Kleckner
%   Ian_Kleckner@URMC.Rochester.edu
%   University of Rochester Medical Center

    %----------------------------------------------------------------------
    % Options
    
    % From C.R. Meyer, H.N. Keiser, Electrocardiogram baseline noise estimation and removal using cubic splines and state-space computation techniques, Comput. Biomed. Res. 10 (1977) 459–470. doi:10.1016/0010-4809(77)90021-0.
    ISOELECTRIC_TIME_BEFORE_RSPIKE_SEC = 66e-3;
    
    
    %----------------------------------------------------------------------
    % Setup
    try
        sampling_period_ECG_sec = X_sec(2) - X_sec(1);
    catch err
        warning('Not enough ECG samples');
        
        % If there are not enough data (duration) then do NOT adjust
        Y_ECG_adjusted = Y_ECG;
        
        return;        
    end
    
    ISOELECTRIC_SAMPLES_BEFORE_RSPIKE = floor(ISOELECTRIC_TIME_BEFORE_RSPIKE_SEC / sampling_period_ECG_sec);


    %----------------------------------------------------------------------
    % Find R spikes from original data
    
%     length(X_sec)
%     X_sec(1)
%     X_sec(end)
%     length(Y_ECG)
%     Minimum_RR_Interval_sec
%     Minimum_R_Prominence_ECG

    duration_data_sec = X_sec(end) - X_sec(1);
    
    if( Minimum_RR_Interval_sec < duration_data_sec )
        
        % Find peaks
        [peak_Y_array, peak_X_array, peak_width_array, peak_prom_array] = ...
                        findpeaks(Y_ECG, X_sec, ...
                        'MinPeakDistance', Minimum_RR_Interval_sec, ...
                        'MinPeakProminence', Minimum_R_Prominence_ECG );

        Npeaks = length(peak_Y_array);
        
        % Only perform adjustment if there are at least 2 R spikes (required
        % for spline fitting)
        if( Npeaks >= 2 )

            % X points are 66 msec before the R spike
            X_for_spline = peak_X_array - ISOELECTRIC_TIME_BEFORE_RSPIKE_SEC;
            Y_for_spline = NaN * X_for_spline;

            % Y points need to be found one by one
            for p = 1:length(X_for_spline)
                k_this_peak     = find( X_sec == peak_X_array(p) );
                
                % If this goes out of range (i.e., R spike happens within
                % the first 66 msec)
                if( k_this_peak - ISOELECTRIC_SAMPLES_BEFORE_RSPIKE <= 0 )
                    X_for_spline(p) = NaN;
                    Y_for_spline(p) = NaN;
                else
                    Y_for_spline(p) =  Y_ECG( k_this_peak - ISOELECTRIC_SAMPLES_BEFORE_RSPIKE );
                end
            end
            
            % Remove invalid entries
            k_invalid = find( isnan(X_for_spline) );
            X_for_spline(k_invalid) = [];
            Y_for_spline(k_invalid) = [];

            try
                % Fit the spline
                Y_spline = spline(X_for_spline, Y_for_spline, X_sec);

                % Perform the adjustment
                Y_ECG_adjusted = Y_ECG - Y_spline;
            catch err
                % Cannot adjust
                Y_ECG_adjusted = Y_ECG;
            end
        else

            % If there are not enough data (R spikes) then do NOT adjust
            Y_ECG_adjusted = Y_ECG;
        end
    else
        % If there are not enough data (duration) then do NOT adjust
        Y_ECG_adjusted = Y_ECG;
    end
end