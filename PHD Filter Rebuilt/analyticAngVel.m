function [ w_out,rot_mode,ksq,coef_list,tau_init ] = analyticAngVel( w_init,moi,time_list )
%ANALYTICANGVEL Analytically propogates the angular velocity of an object
%assuming no torque is being applied.
%   [ w_out,rot_mode,ksq,coef_list,tau_init ] = analyticAngVel( w_init,moi,time_list )
%   Input:
%      w_init - 3x1 array, initial angular velocity in the body-fixed frame
%      moi - 3x3 array, body's principal moments of inertia (I3 < I1 < I2)
%      time_list - nx1 array, times to compute angular velocity at
%   Output:
%      w_out - 3xn array, angular velocity components at each time relative
%         to the body-fixed frame
%      rot_mode - string ("Long Axis", "Short Axis", "Edge"), the mode of 
%         rotation the body is undergoing
%      ksq - scalar, k-squared value used to compute angular velocity
%      coef_list - 3x1 array, list of coefficients used to compute angular
%         velocity
%      tau_init - scalar, value of parameter tau at time time_list(1)
%Author: Alexander Burton
%Created: March 22, 2022

%% Compute System Parameters

% angular momentum
H_vec = moi * w_init;   
H = norm(H_vec);

% rotational kinetic energy
T = 0.5 * w_init' * moi * w_init;   % kinetic energy

% dynamic moment of inertia and effective angular velocity
Id = H^2 / (2 * T);     
we = 2 * T / H;         

% get rotational mode
I1 = moi(1,1);

if Id == I1
    rot_mode = "Edge";
elseif Id < I1
    rot_mode = "Long Axis";
else
    rot_mode = "Short Axis";
end

%% Compute angular velocities

if strcmp(rot_mode,"Long Axis")
    [w_out,ksq,coef_list,tau_init] = lamAnalyticVel(we,Id,w_init,moi,time_list);
elseif strcmp(rot_mode,"Short Axis") || strcmp(rot_mode,"Edge")
    [w_out,ksq,coef_list,tau_init] = samAnalyticVel(we,Id,w_init,moi,time_list);
else
    w_out = [];
    ksq = [];
    coef_list = [];
    tau_init = [];
end


end