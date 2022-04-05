function [P1,period_fft,P2] = fftSingleSided( signal_list,time_list )
% FFTSINGLESIDED Compute the single-sided amplitude spectrum
%   [P1,period_fft,P2] = fftSingleSided( signal_list,time_list )

        num_points = length(signal_list);
        Y = fft(signal_list);
        P2 = abs(Y / num_points);
        P1 = P2(1:(floor(num_points/2) + 1));
        P1(2:end-1) = 2 * P1(2:end-1);
        
        Fs = 1 / (time_list(2) - time_list(1));
        freq_hz = Fs * (0:(num_points / 2))' / num_points;
        period_fft = 1 ./ freq_hz;

 end