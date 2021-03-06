%% Summary
% Finds the possible pseudo-measurements for a given object and light curve
%
% Author: Alex Burton
% Created: March 23, 2021
% Edited: March 31, 2022 - cleaned up for the PHD filter rebuild

clear
clc
close all

set(0,'defaulttextinterpreter','latex')

rng(21)

%% Setup

% object
load '..\Reflectors\tet_obj_asymm.mat'
% load('..\Reflectors\NewDappled_object.mat')

% variables controlling the h5 file name
q = 1;
num_points = 25;
ang_dist = 0.02;

% h5 file
file_name = strcat('Object H5 Files\CheckRun_tet_',int2str(q),...
    'PiBy10_',int2str(num_points),'_',int2str(ang_dist * 100),'_asymm.h5');
% fileName = 'D:\MATLAB\Object H5 Files\CheckRun_NewDappled_1PiBy10_50_2.h5';

% time vector
num_meas = 61;%21;
% dt = 0.5;
tf = 60;%dt * (num_meas - 1);
time_list = linspace(0,tf,num_meas)';%(0:dt:tf)';
dt = time_list(2) - time_list(1);

% reference vectors (inertial frame)
num_ref = 2;
ref_vecs = zeros(3,num_ref);

% objectPos = 4.2e7 * [0 1 0]';
phi = deg2rad(45);
obs_loc = 6378e3 * [cos(phi) 0 sin(phi)]';
obj_obs = obs_loc - object.pos;
r_obs = norm(obj_obs);
ref_vecs(:,1) = obj_obs / r_obs;

delta = q * pi / 10;
ortho_rot = cross(ref_vecs(:,1),[0 1 0]');
ortho_rot = ortho_rot / norm(ortho_rot);
ref_vecs(:,2) = rotMatrix(delta,ortho_rot) * ref_vecs(:,1);
ref_vecs(:,2) = ref_vecs(:,2) / norm(ref_vecs(:,2));

r_sun = 1.5e11;
sun_loc = (r_sun * ref_vecs(:,2)) + object.pos;

% angular velocity
w0 = [0 0.05 1]';%0.25 * [1 2 3]';%[1 1 0]';

% initial orientation
q0 = [0 0 0 1]';

% inertia tensor of the object (Iz < Iy < Ix)
MOI = diag([1.5 1.0 0.5]); %diag([1 2 1]);

% covariance info
sigU = 3.1632e-5;
sigV = sigU;
sigM = 0.1745;
P0 = diag([0.0305 0.0305 0.0305 1 1 1]);


%% Generate True Orientations

% true orientation and angular velocity over time
opt = odeset('RelTol',1e-12,'AbsTol',1e-14);
    
x_start = [q0;w0];
[~,X] = ode45(@(t,x)rotOde(t,x,MOI,zeros(3,1)),time_list,x_start,opt);
Y = X';
q_list = Y(1:4,:);
for i = 1:num_meas
q_list(:,i) = q_list(:,i) / norm(q_list(:,i));
end

w_true = Y(5:7,:);
    

% true body-frame vectors over time
u_true_list = zeros(3,num_meas);
s_true_list = zeros(3,num_meas);
for i = 1:num_meas
   
    % generate vectors by taking ref_vecs from inertial to body frame
    meas_vec = quatRotMatrix(q_list(:,i))' * ref_vecs;
    u_true_list(:,i) = meas_vec(:,1);
    s_true_list(:,i) = meas_vec(:,2);
    
end

%% Create Light Curve

% light curve intensity and relative magnitude
meas_list = zeros(num_meas,1);
mag_list = zeros(num_meas,1);

% generate values
mSun = -26.5;
for i = 1:num_meas

    % object orientation
    rot_matrix = quatRotMatrix(q_list(:,i));

    % light curve intensity/magnitude
    meas_list(i) = object.lambertReflection(obs_loc,sun_loc,rot_matrix);
    mag_list(i) = mSun - 2.5 * log10(meas_list(i));

end
    
plot(time_list,meas_list)

save('..\Input Data\test_lc.mat','meas_list')

%% Find Feasible Orientations Based on LC
    
% storage arrays
num_vec_list = zeros(num_meas,1);
meas_gamma = zeros(num_meas,1);
bounds_list = zeros(num_meas,1);

for i = 1:num_meas

%         tic
    fprintf('Possible Vecs: %d/%d\n',i,num_meas)

    % get gamma value
    gamma = meas_list(i) * pi * r_obs^2;
    meas_gamma(i) = gamma;

    % how much variation allowed in one vector for group with largest Clamb
    big_diff = deg2rad(5e-3);%deg2rad(1e-3); % if spread across both, become sqrt(bigDiff)/vector
    % thetaCrit = acos(sqrt(gamma / max(areaClamb)));
    theta_bound = getEpsilonRange(q * pi/10,meas_list(i) * pi * r_obs^2 / max(object.alignAreaClamb));
    theta_crit = max(theta_bound);
    bound_range = abs(-gamma * tan(theta_crit) * big_diff);
%         boundRange = 1e-3 * gamma;
%         boundRange = 0.85;

    bounds_list(i) = bound_range;

end

% get feasible vectors pairs
vec_list = getPossibleUsVecs_Full(delta,file_name,meas_list,r_obs,bounds_list);

% number of feasible vectors at each time step
for i = 1:num_meas
    num_vec_list(i) = size(vec_list{i},2);
end
    

        

%% Save Results

% save file
% fileSave = strcat('Asymm Tetrahedron\Asymm_tet_',int2str(q),'PiBy10_',...
%     int2str(nPoints),'_',int2str(angDist * 100),'_velX100_',...
%     int2str(floor(abs(w0(1) * 100))),'_',int2str(floor(abs(w0(2) * 100))),'_',...
%     int2str(floor(abs(w0(3) * 100))),'_',int2str(nMeas),'meas_',...
%     int2str(floor(rad2deg(bigDiff) * 1e3)),'e3.mat');

fileSave = strcat('..Input Data\Asymm_tet_',int2str(q),'PiBy10_',...
    int2str(num_points),'_',int2str(ang_dist * 100),'_velX100_',...
    int2str(floor(abs(w0(1) * 100))),'_',int2str(floor(abs(w0(2) * 100))),'_',...
    int2str(floor(abs(w0(3) * 100))),'_',int2str(num_meas),'meas_',...
    int2str(floor(rad2deg(big_diff) * 1e3)),'e3.mat');

save(fileSave)

%{

% assign measurements
uMeasList = zeros(3,nMeas);
sMeasList = zeros(3,nMeas);


zListA = zeros(3 * nRef,nPoints);
zListB = zeros(3 * nRef,nPoints);
for i = 1:nMeas
    
    % true uHat/sHat
    if pmPick == 1
        
        uMeasList(:,i) = uTrueList(:,i);
        sMeasList(:,i) = sTrueList(:,i);
    
    % pseudo-measurements closest to truth
    elseif pmPick == 2
        
        fprintf('Pick Pseuod-Meas:\t%d/%d\n',i,nMeas)
        
        minDist = inf;
        truePair = [uTrueList(:,i) sTrueList(:,i)];
        for j = 1:nVecList(i)
            
            % test vector pair
            testPair = reshape(vecList{i}(2:7,j),3,2);

            % check how close testPair is to truePair
            dist = vecPairDistMeasure(truePair,testPair);
            if dist < minDist
                minDist = dist;
                uMeasList(:,i) = testPair(:,1);
                sMeasList(:,i) = testPair(:,2);
            end

        end
        
    end
    
    % assign measurements
    zListA(:,i) = [uMeasList(:,i);sMeasList(:,i)];
    zListB(:,i) = [sMeasList(:,i);uMeasList(:,i)];
    
end

%% Get Initial Angular Velocity Estimates

measListA = zeros(3,2,nMeas);
for i = 1:nMeas
    measListA(:,1,i) = zListA(1:3,i);
    measListA(:,2,i) = zListA(4:6,i);
end

varList = sqrt(3) * sigM^2 * [1 1]';
wEstArray = pseudo2AngVel(zListA,zListB,timeList,varList,2);

% wHat0 = wEstArray(:,1,1);
% wHat0 = [1.1491 0.6797 0.1319]';
% wHat0 = [0.8412 0.5445 0.0131]';
wHat0 = mean(wEstArray(:,1,:),3);

%% Run Filter

% initial orientation estimates
qHat01 = vecPairs2Quat(reshape(zListA(:,1),3,2),refVecs);
qHat02 = vecPairs2Quat(reshape(zListB(:,1),3,2),refVecs);

% storage for the first track
qHatList1 = zeros(4,nMeas);
qHatList1(:,1) = qHat01;
qHat1 = qHat01;

qPropList1 = zeros(4,nMeas);
qPropList1(:,1) = qHat01;

wHatList1 = zeros(3,nMeas);
wHatList1(:,1) = wHat0;
wHat1 = wHat0;

wPropList1 = zeros(3,nMeas);
wPropList1(:,1) = wHat0;

pkList1 = zeros(6,6,nMeas);
pkList1(:,:,1) = P0;
Pk1 = P0;

% storage for the second track
qHatList2 = zeros(4,nMeas);
qHatList2(:,1) = qHat02;
qHat2 = qHat02;

qPropList2 = zeros(4,nMeas);
qPropList2(:,1) = qHat02;

wHatList2 = zeros(3,nMeas);
wHatList2(:,1) = wHat0;
wHat2 = wHat0;

wPropList2 = zeros(3,nMeas);
wPropList2(:,1) = wHat0;

pkList2 = zeros(6,6,nMeas);
pkList2(:,:,1) = P0;
Pk2 = P0;

% main loop
measPick = zeros(nMeas,2);
measPick(1,:) = [1 2];
for i = 2:nMeas
    
    if i == nMeas || mod(i,100) == 0
        fprintf('%d/%d\n',i,nMeas)
    end
    
    % time span from t(i-1) to t(i)
    tSpan = timeList((i - 1):i);
    
    % propagate
    wReal = wTrue(:,(i - 1));
    [qHatProp1,pkProp1,wHatProp1] = mekfAngEst_Propagate(qHat1,Pk1,wHat1,sigV,sigU,tSpan);
    [qHatProp2,pkProp2,wHatProp2] = mekfAngEst_Propagate(qHat2,Pk2,wHat2,sigV,sigU,tSpan);                                

    % pick measurement closest to the expected
    zHat1 = quatRotMatrix(qHatProp1)' * refVecs;
    zHat2 = quatRotMatrix(qHatProp2)' * refVecs;
    
    zkA = zListA(:,i);
    zkB = zListB(:,i);
    
    [~,order1] = vecPairDistMeasure(zHat1,reshape(zkA,3,2));
    [~,order2] = vecPairDistMeasure(zHat2,reshape(zkA,3,2));
    if order1 == 1
        zk1 = zkA;
        measPick(i,1) = 1;
    else
        zk1 = zkB;
        measPick(i,1) = 2;
    end
    
    if order2 == 1
        zk2 = zkA;
        measPick(i,2) = 1;
    else
        zk2 = zkB;
        measPick(i,2) = 2;
    end
        
    
    % update
%     zk = zList1(:,i);
    [qHat1,wHat1,Pk1] = mekfAngEst_Update(qHatProp1,wHatProp1,pkProp1,zk1,refVecs,Rmeas);
    [qHat2,wHat2,Pk2] = mekfAngEst_Update(qHatProp2,wHatProp2,pkProp2,zk2,refVecs,Rmeas);
%     qHat = qHatProp;
%     Pk = pkProp;

    % store
    qHatList1(:,i) = qHat1;
    qPropList1(:,i) = qHatProp1;
    wPropList1(:,i) = wHatProp1;
    wHatList1(:,i) = wHat1;
    pkList1(:,:,i) = Pk1;
    
    qHatList2(:,i) = qHat2;
    qPropList2(:,i) = qHatProp2;
    wPropList2(:,i) = wHatProp2;
    wHatList2(:,i) = wHat2;
    pkList2(:,:,i) = Pk2;
        
end


%% Compute Orientation Error

% reference axes
vR = refVecs(:,1);%[1 0 0]';

% error in orientation
errRad = zeros(nMeas,2);
errDeg = zeros(nMeas,2);
upRad = zeros(nMeas,2);
upDeg = zeros(nMeas,2);

vArray = zeros(3,5,nMeas);
vDiff = zeros(3,4,nMeas);

wEst = zeros(nMeas,5);
dt = timeList(2) - timeList(1);
for i = 1:nMeas
    
    % compute axes
    vTrue = quatRotMatrix(qList(:,i))' * vR;
    vHat1 = quatRotMatrix(qHatList1(:,i))' * vR;
    vHat2 = quatRotMatrix(qHatList2(:,i))' * vR;
    vProp1 = quatRotMatrix(qPropList1(:,i))' * vR;
    vProp2 = quatRotMatrix(qPropList2(:,i))' * vR;
    
    % store vectors
    vArray(:,:,i) = [vTrue vHat1 vProp1 vHat2 vProp2];
    vDiff(:,:,i) = [(vHat1 - vTrue) (vProp1 - vTrue) (vHat2 - vTrue) (vProp2 - vTrue)];
    
    % estimate angular velocity between updates
    if i ~= 1
        vTp = vArray(:,1,i-1);
        crossTrue = atan2(norm(cross(vTp,vTrue)),dot(vTp,vTrue));
        wEst(i,1) = crossTrue / dt;

        vHp = vArray(:,2,i-1);
        crossHat = atan2(norm(cross(vHp,vHat1)),dot(vHp,vHat1));
        wEst(i,2) = crossHat / dt;
        
        vPp = vArray(:,3,i-1);
        crossProp = atan2(norm(cross(vPp,vProp1)),dot(vPp,vProp1));
        wEst(i,3) = crossProp / dt;
        
        vHp = vArray(:,4,i-1);
        crossHat = atan2(norm(cross(vHp,vHat2)),dot(vHp,vHat2));
        wEst(i,4) = crossHat / dt;
        
        vPp = vArray(:,5,i-1);
        crossProp = atan2(norm(cross(vPp,vProp2)),dot(vPp,vProp2));
        wEst(i,5) = crossProp / dt;

    end
    
    % angle error
    errRad(i,1) = atan2(norm(cross(vTrue,vHat1)),dot(vTrue,vHat1));
    errDeg(i,1) = 360 * errRad(i,1) / (2 * pi);
    
    errRad(i,2) = atan2(norm(cross(vTrue,vHat2)),dot(vTrue,vHat2));
    errDeg(i,2) = 360 * errRad(i,2) / (2 * pi);
    
    % angle change
    upRad(i,1) = atan2(norm(cross(vHat1,vProp1)),dot(vHat1,vProp1));
    upDeg(i,1) = 360 * upRad(i,1) / (2 * pi);
    
    upRad(i,2) = atan2(norm(cross(vHat2,vProp2)),dot(vHat2,vProp2));
    upDeg(i,2) = 360 * upRad(i,1) / (2 * pi);
    
end

%% Predict uHat and sHat

uPred1 = zeros(3,nMeas);
sPred1 = zeros(3,nMeas);

uPred2 = zeros(3,nMeas);
sPred2 = zeros(3,nMeas);
for i = 1:nMeas
    
    % first track
    Rhat1 = quatRotMatrix(qHatList1(:,i));
    
    rotVecs = Rhat1' * refVecs;
    
    uPred1(:,i) = rotVecs(:,1);
    sPred1(:,i) = rotVecs(:,2);
    
    % second track
    Rhat2 = quatRotMatrix(qHatList2(:,i));
    
    rotVecs = Rhat2' * refVecs;
    
    uPred2(:,i) = rotVecs(:,1);
    sPred2(:,i) = rotVecs(:,2);
    
end

%% Covariance Diagonal Values
sigList1 = zeros(6,nMeas);
sigList2 = zeros(6,nMeas);
for i = 1:nMeas
    
    for j = 1:6
        sigList1(j,i) = sqrt(pkList1(j,j,i));
        sigList2(j,i) = sqrt(pkList2(j,j,i));
    end
    
end


%% Reconstruct the ideal wEst based on measPick

wEstIdeal = zeros(3,nMeas,2);
for i = 1:2
    
    % first entry 
    rowInd = 1:3;
    
    % eight possibilities
    measTest = measPick(rowInd,i);
    if isequal(measTest,[1 1 1]')
        wEstIdeal(:,1,i) = wEstArray(:,1,1);

    elseif isequal(measTest,[1 1 2]')
        wEstIdeal(:,1,i) = wEstArray(:,1,2);

    elseif isequal(measTest,[1 2 1]')
        wEstIdeal(:,1,i) = wEstArray(:,1,3);

    elseif isequal(measTest,[1 2 2]')
        wEstIdeal(:,1,i) = wEstArray(:,1,4);

    elseif isequal(measTest,[2 1 1]')
        wEstIdeal(:,1,i) = wEstArray(:,1,5);
        
    elseif isequal(measTest,[2 1 2]')
        wEstIdeal(:,1,i) = wEstArray(:,1,6);

    elseif isequal(measTest,[2 2 1]')
        wEstIdeal(:,1,i) = wEstArray(:,1,7);

    elseif isequal(measTest,[2 2 2]')
        wEstIdeal(:,1,i) = wEstArray(:,1,8);
    end
    
    % final entry
    rowInd = (nMeas - 1):nMeas;
    if isequal(measPick(rowInd,i),[1 1]')
        wEstIdeal(:,end,i) = wEstArray(:,end,1);
    elseif isequal(measPick(rowInd,i),[1 2]')
        wEstIdeal(:,end,i) = wEstArray(:,end,2);
    elseif isequal(measPick(rowInd,i),[2 1]')
        wEstIdeal(:,end,i) = wEstArray(:,end,7);
    elseif isequal(measPick(rowInd,i),[2 2]')
        wEstIdeal(:,end,i) = wEstArray(:,end,8);
    end
    
    % middle entries
    for j = 2:(nMeas - 1)
        
        % row indices
        rowInd = (j - 1):(j + 1);
        
        % eight possibilities
        measTest = measPick(rowInd,i);
        if isequal(measTest,[1 1 1]')
            wEstIdeal(:,j,i) = wEstArray(:,j,1);
            
        elseif isequal(measTest,[1 1 2]')
            wEstIdeal(:,j,i) = wEstArray(:,j,2);
            
        elseif isequal(measTest,[1 2 1]')
            wEstIdeal(:,j,i) = wEstArray(:,j,3);
            
        elseif isequal(measTest,[1 2 2]')
            wEstIdeal(:,j,i) = wEstArray(:,j,4);
            
        elseif isequal(measTest,[2 1 1]')
            wEstIdeal(:,j,i) = wEstArray(:,j,5);
            
        elseif isequal(measTest,[2 1 2]')
            wEstIdeal(:,j,i) = wEstArray(:,j,6);
            
        elseif isequal(measTest,[2 2 1]')
            wEstIdeal(:,j,i) = wEstArray(:,j,7);
            
        elseif isequal(measTest,[2 2 2]')
            wEstIdeal(:,j,i) = wEstArray(:,j,8);
        end
        
    end
    
end

%% Create Plots

% orientation error
figure
subplot(211)
plot(timeList,errDeg(:,1),'b')
grid on
ylabel('Track 1 Error [deg]','FontSize',12)
title('Error in Orientation Estimate','FontSize',14)

subplot(212)
plot(timeList,errDeg(:,2),'b')
grid on
xlabel('time [sec]','FontSize',12)
ylabel('Track 2 Error [deg]','FontSize',12)


% error in quaternion components
figure
subplot(421)
hold on
plot(timeList,qHatList1(1,:) - qList(1,:),'b')
hold off
grid on
ylabel('$\Delta q_1$','FontSize',12)
title('Error in Track 1 $\hat{q}$','FontSize',14)

subplot(422)
hold on
plot(timeList,qHatList2(1,:) - qList(1,:),'r')
hold off
grid on
ylabel('$\Delta q_1$','FontSize',12)
title('Error in Track 2 $\hat{q}$','FontSize',14)

subplot(423)
hold on
plot(timeList,qHatList1(2,:) - qList(2,:),'b')
hold off
grid on
ylabel('$\Delta q_2$','FontSize',12)

subplot(424)
hold on
plot(timeList,qHatList2(2,:) - qList(2,:),'r')
hold off
grid on
ylabel('$\Delta q_2$','FontSize',12)

subplot(425)
hold on
plot(timeList,qHatList1(3,:) - qList(3,:),'b')
hold off
grid on
ylabel('$\Delta q_3$','FontSize',12)

subplot(426)
hold on
plot(timeList,qHatList2(3,:) - qList(3,:),'r')
hold off
grid on
ylabel('$\Delta q_3$','FontSize',12)

subplot(427)
hold on
plot(timeList,qHatList1(4,:) - qList(4,:),'b')
hold off
grid on
xlabel('time [sec]','FontSize',12)
ylabel('$\Delta q_4$','FontSize',12)

subplot(428)
hold on
plot(timeList,qHatList2(4,:) - qList(4,:),'r')
hold off
grid on
xlabel('time [sec]','FontSize',12)
ylabel('$\Delta q_4$','FontSize',12)

% quaternion estimate magnitude
qNorms = zeros(nMeas,1);
qHatNorms = zeros(nMeas,2);
for i = 1:nMeas
    qNorms(i) = norm(qList(:,i));
    qHatNorms(i,1) = norm(qHatList1(:,i));
    qHatNorms(i,2) = norm(qHatList2(:,i));
end

figure 
subplot(211)
plot(timeList,qNorms,'b')
grid on
title('Quaternion Norms','FontSize',14)
ylabel('$q$','FontSize',12)

subplot(212)
hold on
plot(timeList,qHatNorms(:,1),'b')
plot(timeList,qHatNorms(:,2),'r')
hold off
grid on
xlabel('time [sec','FontSize',12)
ylabel('$\hat{q}$','FontSize',12)

% error in angular velocity estimate
figure
subplot(321)
hold on
plot(timeList,wHatList1(1,:) - wTrueProp(1,:),'b')
plot(timeList,3 * sigList1(4,:),'r')
plot(timeList,-3 * sigList1(4,:),'r')
hold off
grid on
title('Error in Track 1 $\hat{\omega}$','FontSize',14)
ylabel('$\hat{\omega}_x$ [rad/s]','FontSize',12)

subplot(322)
hold on
plot(timeList,wHatList2(1,:) - wTrueProp(1,:),'b')
plot(timeList,3 * sigList2(4,:),'r')
plot(timeList,-3 * sigList2(4,:),'r')
hold off
grid on
title('Error in Track 2 $\hat{\omega}$','FontSize',14)
ylabel('$\hat{\omega}_x$ [rad/s]','FontSize',12)

subplot(323)
hold on
plot(timeList,wHatList1(2,:) - wTrueProp(2,:),'b')
plot(timeList,3 * sigList1(5,:),'r')
plot(timeList,-3 * sigList1(5,:),'r')
hold off
grid on
ylabel('$\hat{\omega}_y$ [rad/s]','FontSize',12)

subplot(324)
hold on
plot(timeList,wHatList2(2,:) - wTrueProp(2,:),'b')
plot(timeList,3 * sigList2(5,:),'r')
plot(timeList,-3 * sigList2(5,:),'r')
hold off
grid on
ylabel('$\hat{\omega}_y$ [rad/s]','FontSize',12)

subplot(325)
hold on
plot(timeList,wHatList1(3,:) - wTrueProp(3,:),'b')
plot(timeList,3 * sigList1(6,:),'r')
plot(timeList,-3 * sigList1(6,:),'r')
hold off
grid on
ylabel('$\hat{\omega}_z$ [rad/s]','FontSize',12)
xlabel('time [sec]','FontSize',12)

subplot(326)
hold on
plot(timeList,wHatList2(3,:) - wTrueProp(3,:),'b')
plot(timeList,3 * sigList2(6,:),'r')
plot(timeList,-3 * sigList2(6,:),'r')
hold off
grid on
ylabel('$\hat{\omega}_z$ [rad/s]','FontSize',12)
xlabel('time [sec]','FontSize',12)

% amount of angle change at each update
figure
hold on
plot(timeList,upDeg(:,1),'b')
plot(timeList,upDeg(:,2),'r')
hold off
grid on
xlabel('time [sec]','FontSize',12)
ylabel('$\Delta \phi$ [deg]','FontSize',12)
title('Change in Axis at Update','FontSize',14)

% angular velocity estimate
figure
hold on
plot(timeList(2:end),wEst(2:end,1),'b')
plot(timeList(2:end),wEst(2:end,2),'r')
% plot(timeList(2:end),wEst(2:end,3),'k')
plot(timeList(2:end),wEst(2:end,4),'k')
% plot(timeList(2:end),wEst(2:end,5),'r')
hold off
grid on
xlabel('time [sec]','FontSize',12)
ylabel('$||\hat{\omega}||$ [rad/s]','FontSize',12)
title('Approximated Angular Velocity','FontSize',14)

% plot true reference vector vs estimated in 2D plane
nSet = floor(nMeas / 1);
% indSet = (nPoints - nSet + 1):nPoints;
indSet = 1:nSet;

vTx = reshape(vArray(1,1,indSet),nSet,1);
vTy = reshape(vArray(2,1,indSet),nSet,1);
vTz = reshape(vArray(3,1,indSet),nSet,1);

vHx1 = reshape(vArray(1,2,indSet),nSet,1);
vHy1 = reshape(vArray(2,2,indSet),nSet,1);
vHz1 = reshape(vArray(3,2,indSet),nSet,1);

vPx1 = reshape(vArray(1,3,indSet),nSet,1);
vPy1 = reshape(vArray(2,3,indSet),nSet,1);
vPz1 = reshape(vArray(3,3,indSet),nSet,1);

vHx2 = reshape(vArray(1,4,indSet),nSet,1);
vHy2 = reshape(vArray(2,4,indSet),nSet,1);
vHz2 = reshape(vArray(3,4,indSet),nSet,1);

vPx2 = reshape(vArray(1,5,indSet),nSet,1);
vPy2 = reshape(vArray(2,5,indSet),nSet,1);
vPz2 = reshape(vArray(3,5,indSet),nSet,1);

figure
hold on
plot(vHx1 - vTx,vHy1 - vTy,'b.')
% plot(vHx,vHy,'r')
plot(vHx2 - vTx,vHy2 - vTy,'r.')
hold off
grid on

% plot uHat and sHat vectors
figure
hold on
plot3(uMeasList(1,:),uMeasList(2,:),uMeasList(3,:),'r','Marker','.','MarkerSize',14)
% plot3(uTrueList(1,:),uTrueList(2,:),uTrueList(3,:),'r.')
plot3(uPred1(1,:),uPred1(2,:),uPred1(3,:),'r--','Marker','o')
plot3(uPred2(1,:),uPred2(2,:),uPred1(3,:),'r:','Marker','*')

plot3(sMeasList(1,:),sMeasList(2,:),sMeasList(3,:),'b','Marker','.','MarkerSize',14)
% plot3(sTrueList(1,:),sTrueList(2,:),sTrueList(3,:),'b.')
plot3(sPred1(1,:),sPred1(2,:),sPred1(3,:),'b--','Marker','o')
plot3(sPred2(1,:),sPred2(2,:),sPred2(3,:),'b:','Marker','*')
hold off
grid on

% plot angle estimate standard variations
figure
subplot(311)
hold on
plot(timeList,sigList1(1,:),'b')
plot(timeList,sigList2(1,:),'r')
hold off
grid on
ylabel('$\sigma_{\alpha,x}$ [rad]','FontSize',12)
title('Angle Estimate Std. Dev.','FontSize',14)

subplot(312)
hold on
plot(timeList,sigList1(2,:),'b')
plot(timeList,sigList2(2,:),'r')
hold off
grid on 
ylabel('$\sigma_{\alpha,y}$ [rad]','FontSize',12)

subplot(313)
hold on
plot(timeList,sigList1(3,:),'b')
plot(timeList,sigList2(3,:),'r')
hold off
grid on
ylabel('$\sigma_{\alpha,z}$ [rad]','FontSize',12)
xlabel('time [sec]','FontSize',12)


%% Optional Plots

if lcGen || lcProc
    
    figure 
    yyaxis left
    plot(timeList,measList,'LineWidth',1)
    ylabel('$I_L$','FontSize',12)
    yyaxis right
    plot(timeList,magList,'LineWidth',1)
    set(gca,'ydir','reverse')
    ylabel('Relative Mag.','FontSize',12)
    grid on
    xlabel('time [sec]','FontSize',12)
    % ylabel('$I_L$','FontSize',12)
    title('Tetrahedron Light Curve','FontSize',14)
    
end

%% Save Results
save(fileSave)
%}
%% Extra Functions

% differential equation for the reference quaternion
    function dx = quatGyroOde(~,x,tensor)
        
        % compute the derivative of the orientation quaternion
        quat = x(1:4);
        wVec = x(5:7);
        dq = 0.5 * xi(quat) * -wVec;%wVec;
        
        % compute the derivative of the angular velocity
        Ix = tensor(1,1);
        Iy = tensor(2,2);
        Iz = tensor(3,3);
        
        wx = wVec(1);
        wy = wVec(2);
        wz = wVec(3);
        
        dwx = (Iy - Iz) * (wy * wz) / Ix;
        dwy = (Iz - Ix) * (wx * wz) / Iy;
        dwz = (Ix - Iy) * (wx * wy) / Iz;
        dw = [dwx dwy dwz]';
        
        dx = [dq;dw];
        
    end
    
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
    
% differential equation for the true angular velocity
    function dW = omegaOde(~,~,etaU)
        
        dW = etaU * rand(3,1);
        
    end
    