function [outputArg1,outputArg2] = phdFilterLc( object,vec_list,meas_list,time_list,torque_flag,file_output )
%PHDFILTERLC Uses a PHD filter to find possible attitude time-histories for
%an input light curve and object. Assumes that the object-Sun-observer
%geometry is known.
%   phdFilterLc( object,vec_list,meas_list,time_list,torque_flag,file_output )
%   Input:   object - reflector class object, data about the object being
%               observed
%            vec_list - Nx1 cell array, each cell contains a 7xMi array
%               where each column is a uhat-shat vector pair. Row 1 is the
%               gamma value of the pair; rows 2-4 is the first unit vector;
%               rows 5-7 is the second unit vector.
%            meas_list - Nx1 array, hold the light curve measurements.
%            time_list - Nx1 array, holds time values for each light curve
%               measurement
%            torque_flag - scalar, 0 if assuming no torque
%            file_output - string, saves results to this file


%% Setup

% number of measurements based on time_list
num_meas = length(time_list);

load(file_input,'vecList','nVecList','uTrueList','sTrueList','wTrue','refVecs',...
    'qList','sigM','nRef','nMeas','dt','P0','sigV','sigU',...
    'nMeasLoad','object','obsLoc','sunLoc','timeList','MOI')

%% Generate a Periodogram

%% Generate Initial Guesses

% initial orientation guess taken from vec_list


end