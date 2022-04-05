function [ dI,coef_list ] = getExpectedDerivative( q,w,object,obs_loc,sun_loc )
%GETEXPECTEDDERIVATIVE Computes the expected derivative of the light curve
%at a specified time given an orientation, angular velocity, and
%object/Sun/observer geometry.
%   [ dI,coef_list ] = getExpectedDerivative( q,w,object,obs_loc,su_loc )
%   Inputs: q - 4x1 array, quaternion representing the rotation from the
%              body-fixed frame to the inertial frame
%           w - 3x1 array, angular velocity of the body-fixed frame
%              relative to the inertial frame
%           object - reflector object, object of interest
%           obs_loc - 3x1 array, position vector of the object
%           sun_loc - 3x1 array, position vector of the Sun
%   Output: dI - scalar, expected derivative of the light curve
%           coef_list - 3x1 array, coefficients Cx, Cy, Cz used to compute
%              dI (dI = coef_list' * w)
%
%Author: Alexander Burton
%Created: March 31, 2022

% get problem geometry values
r = norm(object.pos - obs_loc);
u_hat = (obs_loc - object.pos) / r;
u_hat = quatRotMatrix(q)' * u_hat;
ux = u_hat(1);
uy = u_hat(2);
uz = u_hat(3);

s_hat = (sun_loc - object.pos) / norm(sun_loc - object.pos);
s_hat = quatRotMatrix(q)' * s_hat;
sx = s_hat(1);
sy = s_hat(2);
sz = s_hat(3);

% get normal vectors
norm_vecs = object.alignVecs;%quatRotMatrix(q) * object.alignVecs;

% compute coefficients
num_groups = length(object.alignAreaClamb);
coef_list = zeros(3,1);
for i = 1:num_groups

    % normal vector components
    nx = norm_vecs(1,i);
    ny = norm_vecs(2,i);
    nz = norm_vecs(3,i);

    % dot products
    u_dot = dot(u_hat,norm_vecs(:,i));
    s_dot = dot(s_hat,norm_vecs(:,i));

    % update coefficients
    if u_dot > 0 && s_dot > 0

        coef_list(1) = coef_list(1) + object.alignAreaClamb(i) * ...
            ((ny * uz - nz * uy) * s_dot + (ny * sz - nz * sy) * u_dot);

        coef_list(2) = coef_list(2) + object.alignAreaClamb(i) * ...
            ((nz * ux - nx * uz) * s_dot + (nz * sx - nx * sz) * u_dot);

        coef_list(3) = coef_list(3) + object.alignAreaClamb(i) * ...
            ((nx * uy - ny * ux) * s_dot + (nx * sy - ny * sx) * u_dot);

    end

end

coef_list = coef_list / (pi * r^2);

% compute derivative
dI = coef_list' * w;

end