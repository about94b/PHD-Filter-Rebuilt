function [ axes_list,theta_list,lam_error,sam_error ] = ...
    williamsGeneral( input_lc,object,obs_loc,sun_loc,num_axes,num_theta,num_angles )
%WILLIAMSGENERAL Locates possible angular momentum-angle combinations for
%an input light curve and object. Returns results for both long-axis mode
%(LAM) and short-axis mode (SAM) rotations.
%
%Uses sphere_fibonacci_grid_points by John Burkardt (released under GNU
%LGPL license, https://people.sc.fsu.edu/~jburkardt/m_src/sphere_fibonacci
%_grid/sphere_fibonacci_grid.html)
%   [ axes_list,theta_list,lam_error,sam_error ] = 
%               williamsGeneral( input_lc,object,obs_loc,sun_loc,num_axes,num_theta,num_angles )
%   Inputs: input_lc - nx1 array, intensity light curve to process
%           object - reflector object, object producing the light curve
%           obs_loc - 3x1 array, observer location
%           sun_loc - 3x1 array, Sun location
%           num_axes - scalar, the number of angular momentum axes to test
%           num_theta - scalar, the number of theta angles tested per axis
%           num_angles - scalar, the number of other Euler angles tested
%               per axis
%   Outputs: axes_list, 3 x num_axes array, unit vectors for the angular
%              momentum vectors tested (one vector per column)
%           theta_list - num_theta x 1 array, list of theta values
%           lam_error - num_axes x num_theta array, the error fraction
%               between each axis-angle combination's max/min intensity
%               ratio and the input_lc max/min intensity ratio for
%               LAM rotations.
%           sam_error - num_axes x num_theta array, as lam_error but for
%               SAM rotations.
%
%Author: Alexander Burton
%Created: April 1, 2022

%% Compute the Intensity Ratio for Axis-Theta Pairs

% generate angular momentum axes in inertial space
axes_list = sphere_fibonacci_grid_points(num_axes)';

% the divisions to use for testing the Euler angles
theta_list = linspace(0,pi/2,num_theta)';

phi_list = linspace(0,2 * pi,num_angles + 1)';
phi_list(end) = [];

psi_list = linspace(0,2 * pi,num_angles + 1)';
psi_list(end) = [];

% step through axes
axes_max_long = zeros(num_axes,num_theta);
axes_min_long = inf(num_axes,num_theta);

axes_max_short = zeros(num_axes,num_theta);
axes_min_short = inf(num_axes,num_theta);

fprintf('williamsGeneral: Axis Processing\n')

for i = 1:num_axes

    fprintf('\tAxis %d/%d\n',i,num_axes)

    % set up angular momentum coordinate system
    trial_axis = axes_list(:,i);

    if norm(trial_axis - [0 0 1]') < 1e-6

        Rh_long = eye(3);

    else

        h1 = cross([0 0 1]',trial_axis);
        h1 = h1 / norm(h1);

        h2 = cross(trial_axis,h1);
        h2 = h2 / norm(h2);

        Rh_long = [h1'; h2'; trial_axis'];

    end

    Rh_short = [0 0 1; 1 0 0; 0 1 0] * Rh_long;

    % step through theta angles
    for j = 1:num_theta

        theta = theta_list(j);
    
        % middle rotation for two modes
        R_theta_long = rotMx(theta);
        R_theta_short = rotMy(theta);

        % step through possible object orientations
        for k = 1:num_angles

            phi = phi_list(k);

            % first rotation for the two rotation modes
            R_phi_long = rotMz(phi);
            R_phi_short = rotMx(phi);

            % step through possible orientations about the central axis
            for p = 1:num_angles

                psi = psi_list(p);

                % final rotation for two axis modes
                R_psi_long = rotMz(psi);
                R_psi_short = rotMx(psi);

                % angular momentum to body-frame rotation matrices
                R313_long = R_psi_long * R_theta_long * R_phi_long;
                R121_short = R_psi_short * R_theta_short * R_phi_short;

                % full rotation matrix for both rotation modes
                dcm_long = Rh_long' * R313_long';
                dcm_short = Rh_short' * R121_short';

                % compute reflection
                I_short = object.lambertReflection(obs_loc,sun_loc,dcm_short);
                I_long = object.lambertReflection(obs_loc,sun_loc,dcm_long);

                % check for max and min
                if I_short < axes_min_short(i,j)
                    axes_min_short(i,j) = I_short;
                end

                if I_short > axes_max_short(i,j)
                    axes_max_short(i,j) = I_short;
                end

                if I_long < axes_min_long(i,j)
                    axes_min_long(i,j) = I_long;
                end

                if I_long > axes_max_long(i,j)
                    axes_max_long(i,j) = I_long;
                end

            end


        end

    end

end

%% Compute the Error in each of the Axis-Theta Pairs

fprintf('williamsGeneral: Compute Errors\n')

% Get maximum and minimum points on the light curve
I_max = max(input_lc);
I_min = min(input_lc);
I_ratio_true = I_max / I_min;

% compute the found ratios
I_ratio_test_short = axes_max_short ./ axes_min_short;
I_ratio_test_long = axes_max_long ./ axes_min_long;

% LAM rotation errors
ratio_diff_long = abs(I_ratio_true - I_ratio_test_long);
lam_error = ratio_diff_long ./ I_ratio_true;

% SAM rotation errors
ratio_diff_short = abs(I_ratio_true - I_ratio_test_short);
sam_error = ratio_diff_short ./ I_ratio_true;

end