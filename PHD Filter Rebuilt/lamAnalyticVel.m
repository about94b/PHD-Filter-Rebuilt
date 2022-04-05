function [ w_out,ksq,coef_list,tau_init ] = lamAnalyticVel( we,Id,w_init,moi,time_list )
%LAMANALYTICVEL Analytically propogates the angular velocity of an object
%assuming no torque is being applied for long-axis mode rotations.
%   [ w_out,ksq,coef_list,tau_init ] = lamAnalyticVel( w_init,moi,time_list )
%   Input:
%      we - scalar, effective angular velocity
%      Id - scalar, dynamic moment of inertia
%      w_init - 3x1 array, initial angular velocity in the body-fixed frame
%      moi - 3x3 array, body's principal moments of inertia (I3 < I1 < I2)
%      time_list - nx1 array, times to compute angular velocity at
%   Output:
%      w_out - 3xn array, angular velocity components at each time relative
%         to the body-fixed frame
%      ksq - scalar, k-squared value used to compute angular velocity
%      coef_list - 3x1 array, list of coefficients used to compute angular
%         velocity
%      tau_init - scalar, value of parameter tau at time time_list(1)
%Author: Alexander Burton
%Created: March 22, 2022

%% Compute Coefficients

% extract moments of inertia
I1 = moi(1,1);
I2 = moi(2,2);
I3 = moi(3,3);

% k-squared
ksq = (I2 - I1) * (Id - I3) / ((I1 - I3) * (I2 - Id));

% coefficients
coef_wx = we * sqrt(Id * (Id - I3) / (I1 * (I1 - I3)));
if w_init(3) < 0
    coef_wy = -we * sqrt(Id * (Id - I3) / (I2 * (I2 - I3)));
    coef_wz = -we * sqrt(Id * (I2 - Id) / (I3 * (I2 - I3)));
else
    coef_wy = we * sqrt(Id * (Id - I3) / (I2 * (I2 - I3)));
    coef_wz = we * sqrt(Id * (I2 - Id) / (I3 * (I2 - I3)));
end

coef_list = [coef_wx coef_wy coef_wz]';

%% Get initial tau

% get phi value
sn_init = w_init(1) / coef_wx;
cn_init = w_init(2) / coef_wy;
phi_init = trigInv(cn_init,sn_init,1e-7);

% compute tau_init
tau_init = ellipticF(phi_init,ksq);

% compute initial angular velocity to double-check
% w_test = zeros(3,1);
% w_test(1) = coef_wx * jacobiSN(tau_init,ksq);
% w_test(2) = coef_wy * jacobiCN(tau_init,ksq);
% w_test(3) = coef_wz * jacobiDN(tau_init,ksq);

%% Propagate

% convert time_list to tau
tau_list = tau_init + we * ...
    sqrt(Id * (I2 - Id) * (I1 - I3) / (I1 * I2 * I3)) * (time_list - time_list(1));

% individual components
wx_out = coef_wx * jacobiSN(tau_list,ksq);
wy_out = coef_wy * jacobiCN(tau_list,ksq);
wz_out = coef_wz * jacobiDN(tau_list,ksq);

% construct list
w_out = [wx_out wy_out wz_out]';

end