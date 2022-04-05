function [ chisq_list ] = getPeriodogram( signal_list,time_list,...
    period_list,num_bins )
%GETPERIODOGRAM Compute the periodogram for a signal (given in signal_list)
%given a set of test periods. Assumes bins are of uniform width.
%   [ chisq_list ] = getPeriodogram( signal_list,time_list,period_list,num_bins )
%   Inputs: signal_list - nx1 array, list of signal values
%           time_list - nx1 array, list of time-values for the signal points
%           period_list - mx1 array, periods to test on the signal
%           num_bins - scalar, number of bins used to compute chi-squared
%   Outputs: chisq_list - mx1 array, chi-squared values for each time peroid
% Author: Alex Burton
% Created: February 3, 2022

%% Initialize

% number of periods to test
num_period = length(period_list);

% statistics of the signal
x_bar = mean(signal_list);
std_tot = std(signal_list);

%% Iterate Through Periods

chisq_list = zeros(num_period,1);
for i = 1:num_period

     % fold the signal
    period_test = period_list(i);
    [signal_fold,time_fold] = shift_points(signal_list,time_list,period_test);
    time_bins = discretize(time_fold,num_bins);
    
    % compute chisq by stepping through bins
    chisq = 0;
    for j = 1:num_bins

        signal_bin = signal_fold(time_bins == j);
        x_bin = mean(signal_bin);
        num_in = length(signal_bin);
        std_bin = std_tot / sqrt(num_in);

        if num_in  > 0
            chisq = chisq + (x_bar - x_bin)^2 / std_bin^2;
        end

    end

    % store value
    chisq_list(i) = chisq;


end


%% Function to Shift Points
% function [signal_shift,time_shift] = shift_points(signal_orig,time_orig,p_test)
% 
%     % shift time points
%     time_shift = time_orig;
%     num_shift = length(time_shift);
%     for q = 1:num_shift
%         
%         while time_shift(q) > p_test
%             time_shift(q) = time_shift(q) - p_test;
%         end
%         
%     end
%     
%     % rearrange
%     [time_shift,ind_shift] = sort(time_shift);
%     signal_shift = signal_orig(ind_shift);
% 
% end

end