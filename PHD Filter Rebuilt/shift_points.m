function [signal_shift,time_shift,loop_track,loop_signals,loop_time] = ...
                            shift_points(signal_orig,time_orig,p_test)
%SHIFT_POINTS Function to "fold" a signal into a length of time equal to
%p_test.
%   [signal_shift,time_shift,loop_track,loop_signals,loop_time] = 
%                           shift_points(signal_orig,time_orig,p_test)
%   Inputs: signal_orig - nx1 array, original (unfolded) signal
%            time_orig - nx1 array, times for each point in signal_orig
%            p_test - scalar, period of time the signal is folded into
%   Outputs: signal_shift - nx1 array, signal points arranged in the order
%               they pass through p_test
%            time_shift - nx1 array, time_orig "looped" through a set of
%               values between 0 and p_test. Time values corresponding to
%               points in signal_shift
%            loop_track - nx1 array, records how many times a point in the
%               signal loops through p_test
%            loop_signals - px1 cell array, each entry is one loop of the
%               signal through p_test
%            loop_time - px1 cell array, each entry is the time values for
%                an entry in loop_signals
%
% Author: Alexander Burton
% Created: February 10, 2022

% shift time points
time_shift = time_orig;
if time_shift(1) ~= 0
    time_shift = time_shif - time_shift(1);
end

num_shift = length(time_shift);
loop_track = zeros(num_shift,1);
for q = 1:num_shift
    
    while time_shift(q) > p_test
        time_shift(q) = time_shift(q) - p_test;
        loop_track(q) = loop_track(q) + 1;
    end
    
end

% rearrange
[time_shift,ind_shift] = sort(time_shift);
signal_shift = signal_orig(ind_shift);
loop_track = loop_track(ind_shift);

% compute loop tracks
num_loops = max(loop_track) + 1;
loop_signals = cell(num_loops,1);
loop_time = cell(num_loops,1);
for i = 1:num_loops

    loop_signals{i} = signal_shift(loop_track == (i - 1));
    loop_time{i} = time_shift(loop_track == (i - 1));

end

end