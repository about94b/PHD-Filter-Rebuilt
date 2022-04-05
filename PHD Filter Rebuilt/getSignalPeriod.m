function [ period_out ] = getSignalPeriod( signal_list,time_list,search_res,num_bins )
%GETSIGNALPERIOD Finds the main period of an input signal using the FFT and
%periodogram.
%   [ period_out ] = getSignalPeriod( signal_list,time_list,min_res,num_bins )
%   Inputs: signal_list - nx1 array, signal being analyzed
%           time_list - nx1 array, time stamps for values in signal_list
%           search_res - scalar, search resolution for the periodogram 
%               is set to this. If the FFT already gives this resolution or
%               better, the periodogram will not be used.
%           num_bins - scalar, number of bins used in the periodogram.
%               Default value is 25.
%   Output: period_out - scalar, main period of the signal
%
%Author: Alexander Burton
%Created: April 1, 2022

% check the number of bins
if nargin < 4
    num_bins = 25;
end

% initial fft results
[amp_spectrum,period_fft] = fftSingleSided(signal_list,time_list);
[fft_peaks,fft_peak_locs_init] = findpeaks(amp_spectrum(2:end),'SortStr','descend','NPeaks',10);
fft_peak_locs_init = fft_peak_locs_init + 1;

period_fft_peaks_init = period_fft(fft_peak_locs_init);

% combine resonant peaks into a single large peak at highest period
num_peaks_init = length(fft_peak_locs_init);
fft_peaks_combined = fft_peaks;
del_fft_peaks = [];
for i = 1:num_peaks_init
    
    % find integer multiples of a test period
    test_period = period_fft_peaks_init(i);

    % divide through
    period_ratios = period_fft_peaks_init / test_period;
    period_remain = mod(period_ratios,1);

    % find integer multiples
    for j = 1:num_peaks_init

        % do not compare a period to itself or ratios for deletion
        if j == i
            skip_flag = 1;
        elseif ismember(i,del_fft_peaks)
            skip_flag = 1;
        elseif ismember(j,del_fft_peaks)
            skip_flag = 1;
        else
            skip_flag = 0;
        end

        if ~skip_flag

            if period_remain(j) < 1e-8
                
                % add period j to list for deletion
                del_fft_peaks = cat(1,del_fft_peaks,i);

                % combine peak i with peak j
                fft_peaks_combined(j) = fft_peaks_combined(j) + fft_peaks_combined(i);

            end

        end

    end

end

% get highest combined peak
% period_fft_peaks = period_fft_peaks_init;
fft_peak_locs = fft_peak_locs_init;

% period_fft_peaks(del_fft_peaks) = [];
fft_peaks_combined(del_fft_peaks) = [];
fft_peak_locs(del_fft_peaks) = [];

[~,ind_peak_max] = max(fft_peaks_combined);
peak_loc = fft_peak_locs(ind_peak_max);

% check if the periodogram needs to be used
up_res = period_fft(peak_loc - 1) - period_fft(peak_loc);
down_res = period_fft(peak_loc) - period_fft(peak_loc + 1);
if up_res < search_res && down_res < search_res

    % if already good resolution, output period
    period_out = period_fft(peak_loc);

else

    % generate periods to search
    search_start = period_fft(peak_loc + 1);
    search_end = period_fft(peak_loc - 1);
    period_list = (search_start:search_res:search_end)';

    % create periodogram
    chisq_list = getPeriodogram(signal_list,time_list,period_list,num_bins);
    [chisq_peak,chisq_loc] = findpeaks(chisq_list,'MinPeakProminence',20);

    % find the most prominent peak
    period_peaks = period_list(chisq_loc);
    [~,best_peak] = max(chisq_peak);
    
    % output period
    period_out = period_peaks(best_peak);

end

end