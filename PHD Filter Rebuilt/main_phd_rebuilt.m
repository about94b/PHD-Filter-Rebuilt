%% Summary
% The script function for running the redone PHD filter.
%
% Author: Alexander Burton
% Created: March 30, 2022

clear
clc
close all
set(0,'DefaultTextInterpreter','Latex')

%% File Names

% load data
% mat_name = 'Input Data\LongCase_NewDappled_1PiBy10_50_2_velX100_5_100_0_401meas_5e3.mat';
mat_name = 'Input Data\LongCase_NewDappled_1PiBy10_50_2_run7.mat';
load(mat_name,'vecList','nVecList','uTrueList','sTrueList','wTrue','refVecs',...
    'qList','sigM','nRef','nMeas','dt','P0','sigV','sigU',...
    'nMeasLoad','object','obsLoc','sunLoc','timeList','MOI')

load 'Reflectors\tet_obj_asymm.mat'
load 'Input Data\test_lc.mat'

fprintf('Old LC\n')
time_list_old = timeList;
num_meas_old = length(time_list_old);
old_lc = zeros(num_meas_old,1);
for i = 1:num_meas_old
    old_lc(i) = object.lambertReflection(obsLoc,sunLoc,quatRotMatrix(qList(:,i)));
end
MOI = diag([1.5 1 0.5]);

% redo light curve
fprintf('Redo\n')
q_init = qList(:,1);
w_init = [0 0.05 1]';%wTrue(:,1);
x_init = [q_init;w_init];
num_meas = 61;%ceil(num_meas_old * 5);
time_list = linspace(0,60,num_meas)';%linspace(time_list_old(1),time_list_old(end),num_meas)';
opt = odeset('RelTol',1e-12,'AbsTol',1e-14);
[t_out,y_out] = ode45(@(t,x) rotOde(t,x,MOI,zeros(3,1)),time_list,x_init,opt);
q_list = y_out(:,1:4)';
w_list = y_out(:,5:7)';
true_lc = zeros(num_meas,1);
% clean_lc = zeros(num_meas,1);
for i = 1:num_meas
%     clean_lc(i) = object.lambertReflection(obsLoc,sunLoc,quatRotMatrix(q_list(:,i)));
    true_lc(i) = object.lambertReflection(obsLoc,sunLoc,quatRotMatrix(q_list(:,i)));% + 8e-14 * rand(1) - 4e-14;
end

%% Check LC Derivative Condition at Each Point

% use FFT and Periodogram to get main period
fprintf('Get Period\n')
main_period = getSignalPeriod(true_lc,time_list,1e-4,25);

% fold the signal back on itself
fprintf('Shift\n')
[shift_lc,time_shift] = shift_points(true_lc,time_list,main_period);

% fit the folded signal with an 8-term Fourier fit
fprintf('Fit\n')
[fit_lc,gof] = fit(time_shift,shift_lc,'Fourier8');
fit_overall = sqrt(mean((true_lc - fit_lc(time_list)).^2 ./ true_lc.^2));

% derivative of modeled light curve
fprintf('Differentiate\n')
fit_deriv = differentiate(fit_lc,time_list);

lc_deriv = zeros(num_meas,1);
for i = 1:num_meas
    
    q = q_list(:,i);
    w = w_list(:,i);
    lc_deriv(i) = getExpectedDerivative(q,w,object,obsLoc,sunLoc);

end

%% Apply the Williams Method

% run the method
num_axes = 10;
num_theta = 10;
num_angles = 5;
[axes_list,theta_list,lam_error,sam_error] = williamsGeneral(true_lc,object,obsLoc,sunLoc,num_axes,num_theta,num_angles);

% compute the angle between each axis and the true angular momentum vector
H_vec = MOI * quatRotMatrix(q_init) * w_init;
axis_true = H_vec / norm(H_vec);
axes_error = zeros(num_axes,1);
for i = 1:num_axes

    axis_dot = dot(axis_true,axes_list(:,i));
    axes_error(i) = acos(axis_dot);

end

%% compute the true averaged theta angle
T = 0.5 * w_init' * MOI * w_init;
H = norm(H_vec);
Id = H^2 / (2 * T);
theta_true_list = zeros(num_meas,1);
axis_true_list = zeros(3,num_meas);
H_vec_list = zeros(3,num_meas);
H_mag_list = zeros(num_meas,1);
main_axis_list = zeros(3,num_meas);
for i = 1:num_meas

    if Id > MOI(2,2)
        main_axis = quatRotMatrix(q_list(:,i)) * [1 0 0]';
    else
        main_axis = quatRotMatrix(q_list(:,i)) * [0 0 1]';
    end

    main_axis_list(:,i) = main_axis;

    H_vec_temp = quatRotMatrix(q_list(:,i)) * MOI * w_list(:,i);
    axis_true_list(:,i) = H_vec_temp / norm(H_vec_temp);
    H_vec_list(:,i) = H_vec_temp;
    H_mag_list(i) = norm(H_vec_temp);

    theta_true_list(i) = acos(dot(main_axis,axis_true_list(:,i)));

end

theta_mean = mean(theta_true_list);



%% Plots

fprintf('Plot\n')

% compare light curves
figure
hold on
plot(time_list_old,old_lc,'b','LineWidth',1)
% plot(time_list,clean_lc,'b','LineWidth',1)
plot(time_list,true_lc,'r','LineWidth',1)
hold off
grid on

% williams plots
figure
surf(theta_list - theta_mean,axes_error,lam_error,'EdgeColor','interp','FaceColor','interp')
ylabel('Axis Error [rad]')
xlabel('$\theta$ [rad]')
zlabel('Perc. Error')
title('LAM Error Fraction (surf)')
colorbar

figure
surf(theta_list - theta_mean,axes_error,sam_error,'EdgeColor','none','FaceColor','interp')
ylabel('Axis Error [rad]')
xlabel('$\theta$ [rad]')
zlabel('Perc. Error')
title('SAM Error Fraction (surf)')
colorbar

figure
subplot(321)
plot(time_list,axis_true_list(1,:),'b','LineWidth',1)
grid on
ylabel('$H_x$')
title('True $H$')

subplot(322)
plot(time_list,main_axis_list(1,:),'b','LineWidth',1)
grid on
ylabel('$A_x$')
title('Main Axis')

subplot(323)
plot(time_list,axis_true_list(2,:),'b','LineWidth',1)
grid on
ylabel('$H_y$')

subplot(324)
plot(time_list,main_axis_list(2,:),'b','LineWidth',1)
grid on
ylabel('$A_y$')

subplot(325)
plot(time_list,axis_true_list(3,:),'b','LineWidth',1)
grid on
ylabel('$H_z$')

subplot(326)
plot(time_list,main_axis_list(3,:),'b','LineWidth',1)
grid on
ylabel('$A_z$')


figure
subplot(311)
hold on
plot(time_list,axis_true_list(1,:),'b','LineWidth',1)
plot(time_list,main_axis_list(1,:),'r--','LineWidth',1)
hold off
grid on
legend('$\hat{H}$','Body','Interpreter','Latex')
ylabel('$\hat{x}$')
title('Main Axis vs. $\hat{H}$')

subplot(312)
hold on
plot(time_list,axis_true_list(2,:),'b','LineWidth',1)
plot(time_list,main_axis_list(2,:),'r--','LineWidth',1)
hold off
grid on
ylabel('$\hat{y}$')

subplot(313)
hold on
plot(time_list,axis_true_list(3,:),'b','LineWidth',1)
plot(time_list,main_axis_list(3,:),'r--','LineWidth',1)
hold off
grid on
ylabel('$\hat{z}$')
xlabel('time [sec]')


figure
subplot(311)
plot(time_list,-axis_true_list(1,:) + main_axis_list(1,:),'b','LineWidth',1)
grid on
ylabel('$\hat{x}$')
title('Main Axis vs. $\hat{H}$ Difference')

subplot(312)
plot(time_list,main_axis_list(2,:) - axis_true_list(2,:),'b','LineWidth',1)
grid on
ylabel('$\hat{y}$')

subplot(313)
plot(time_list,main_axis_list(3,:) - axis_true_list(3,:),'b','LineWidth',1)
grid on
ylabel('$\hat{z}$')
xlabel('time [sec]')

% LC fit
figure
subplot(211)
hold on
% plot(time_list_old,old_lc,'b','LineWidth',1)
plot(time_list,true_lc,'b','LineWidth',1)
plot(time_list,fit_lc(time_list),'r','LineWidth',1)
hold off
grid on
legend('True','Fit')
title('Light Curves')

subplot(212)
hold on
plot(fit_lc,time_shift,shift_lc)%,'b.','MarkerSize',10)
grid on

% LC derivative
figure
hold on
plot(time_list,fit_deriv,'b','LineWidth',1)
plot(time_list,lc_deriv,'r','LineWidth',1)
hold off
grid on
legend('Model','Expected')
title('LC Derivatives')

% angular velocity
figure
subplot(311)
plot(time_list,w_list(1,:),'b','LineWidth',1)
grid on
ylabel('$\omega_x$ [rad/s]')

subplot(312)
plot(time_list,w_list(2,:),'b','LineWidth',1)
grid on
ylabel('$\omega_y$ [rad/s]')

subplot(313)
plot(time_list,w_list(3,:),'b','LineWidth',1)
grid on
ylabel('$\omega_z$ [rad/s]')
xlabel('time [sec]')

% angular momentum
figure
plot(time_list,H_mag_list)
grid on
xlabel('time [sec]')
ylabel('$H$')
title('Angular Momentum Magnitude')


%% Extra Functions
    
% differential equation for the reference quaternion
    function dq = quatAngEstOde(~,q,w)
        
        dq = 0.5 * xi(q) * w;
        
    end

% differential equation for the spinner with rotational kinematics modeled
    function dx = rotOde(t,x,inert,torque)
        
        % angular velocity derivative
        w = x(5:7);
        if inert == 0
            dw = zeros(3,1);
        else
            dw = inert \ (torque + cross(w,inert * w));
        end
        
        % quaternion derivative
        q = x(1:4);
        dq = quatAngEstOde(t,q,w);
        dx = [dq;dw];
    
    end

% cross product matrix
    function matrix = crossMat(vec)
        
        matrix = [0 -vec(3) vec(2); vec(3) 0 -vec(1); -vec(2) vec(1) 0];
        
    end

% xi function (see Crassidis and Junkins p. 612)
    function X = xi(q)
        
        X = zeros(4,3);
        
        X(1:3,:) = q(4) * eye(3) + crossMat(q(1:3));
        X(4,:) = -q(1:3)';
        
    end