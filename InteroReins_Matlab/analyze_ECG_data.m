% Ian Kleckner
% ian.kleckner@gmail.com
% Cancer Control Division
% University of Rochester Medical Center
%
% Program to analyze ECG data

%% Input

filename_input = 'C:\Users\ikleckner\Desktop\out-ECG-Trial_01.csv';

Minimum_RR_Interval_sec = 0.25;
Minimum_R_Prominence_mV = 0.1;

%% Load data

data        = csvread(filename_input, 1, 0);
data_time   = data(:,1);
data_ECG_raw    = data(:,2);


%% Process data

data_ECG_adj = ECG_adjust_baseline_spline( data_time, data_ECG_raw, Minimum_RR_Interval_sec, Minimum_R_Prominence_mV );

wt = modwt(data_ECG_adj,5);
wtrec = zeros(size(wt));
wtrec(4:5,:) = wt(4:5,:);
y = imodwt(wtrec,'sym4');
y = abs(y).^2;
[qrspeaks,locs] = findpeaks(y,data_time, 'MinPeakHeight', Minimum_R_Prominence_mV/5, ...
    'MinPeakDistance', Minimum_RR_Interval_sec);

%% Display

figure
plot(data_time, y)
hold on
plot(locs,qrspeaks,'ro')
xlabel('Seconds')
title('R Peaks Localized by Wavelet Transform with Automatic Annotations')


%% Display raw data

figure
plot(data_time, data_ECG_adj)

peak_ECG = NaN * locs;
for p = 1:length(locs)
    peak_ECG(p) = data_ECG_adj(  data_time==locs(p)  );
end

hold('all');

plot(locs, peak_ECG, 'or');




%% Output