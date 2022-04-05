function [ q ] = vecPairs2Quat( vecPair1,vecPair2 )
%vecPairs2Quat Gives the quaternion to rotate from vecPair1 to vecPair2.
%The vector pairs are 3x2 arrays.
%   [ q ] = vecPairs2Quat( vecPair1,vecPair2 )

% separate out vectors
u1 = vecPair1(:,1);
s1 = vecPair1(:,2);

u2 = vecPair2(:,1);
s2 = vecPair2(:,2);

% change in vectors
dU = u2 - u1;
dux = dU(1);
duy = dU(2);
duz = dU(3);

dS = s2 - s1;
dsx = dS(1);
dsy = dS(2);
dsz = dS(3);

% check for components approximately zero
dsZeros = (abs(dS) < 1e-8);
duZeros = (abs(dU) < 1e-8);

compZeros = zeros(3,1);
for i = 1:3
    if dsZeros(i) && duZeros(i)
        compZeros(i) = 1;
    end
end

nCompZeros = length(find(compZeros));

% get rotation axis
if nCompZeros == 1
    
    % if both have no change in one direction, then that is the rotation
    % axis
    rotAxis = zeros(3,1);
    rotAxis(compZeros == 1) = 1;
    
elseif nCompZeros == 3
    
    % rotation axis is arbitrary, since rotation angle is zero
    rotAxis = [1 0 0]';
    
elseif nCompZeros == 0
    
    % check if any component of the rotation axis is zero
    rotZero = zeros(3,1);
    
    for i = 1:3
        
        % automatically nonzero if change in other directions is zero
        sTemp = dS;
        sTemp(i) = [];
        
        szTemp = dsZeros;
        szTemp(i) = [];
        
        uTemp = dU;
        uTemp(i) = [];
        
        uzTemp = duZeros;
        uzTemp(i) = [];
        
        % check ratios to find if nonzero or not
        if isempty(find(szTemp,1)) && isempty(find(uzTemp,1))
            
            sRatio = -sTemp(1) / sTemp(2);
            uRatio = -uTemp(1) / uTemp(2);
            
            if abs(sRatio - uRatio) < 1e-8
                rotZero(i) = 1;
            end
            
        end
        
    end
    
    nZeros = length(find(rotZero));
    
    % compute rotation axis
    if nZeros == 1
        
        j = find(rotZero);
        uTemp = dU;
        uTemp(j) = [];
        
        rk = -uTemp(2) / uTemp(1);
        
        if j == 1
            rotAxis = [0 rk 1]';
        elseif j == 2
            rotAxis = [rk 0 1]';
        else
            rotAxis = [rk 1 0]';
        end
        
        rotAxis = rotAxis / norm(rotAxis);
        
    elseif nZeros == 0
        
        % let r3 = 1, then normalize at end
        r3 = 1;
        
        % compute r1
        numer = (duy * dsz / dsy - duz);
        denom = (dux - duy * dsx / dsy);    % this will not be zero unless r3 = 0
        r1 = numer / denom * r3;
        
        % compute r2
        r2 = -(dsx * r1 + dsz * r3) / dsy;
        
        % normalize
        rotAxis = [r1 r2 r3]';
        rotAxis = rotAxis / norm(rotAxis);
        
    % check if dU and dS are equal and opposite
    elseif norm(dU + dS) < 1e-8
        
        rotAxis = (vecPair1(:,1) + vecPair1(:,2)) / 2;
        rotAxis = rotAxis / norm(rotAxis);
        
    else
        
        fprintf('vecPairs2Quat: Wrong number of zeros in rotAxis\n')
        
    end
    
else
    
    fprintf('vecPairs2Quat: Cannot have 2 shared zeros\n')
    
end
        
        
% get angle of rotation
uProj1 = u1 - dot(u1,rotAxis) * rotAxis;
uProj2 = u2 - dot(u2,rotAxis) * rotAxis;

dotProd = dot(uProj1,uProj2) / (norm(uProj1) * norm(uProj2));
crossProd = cross(uProj1,uProj2) / (norm(uProj1) * norm(uProj2));
crossSign = sign(dot(crossProd,rotAxis)) * norm(crossProd);
% angCross = acos(dotProd / (norm(uProj1) * norm(uProj2)));
% angCross = atan2(norm(crossProd),dotProd);
angCross = atan2(crossSign,dotProd);

%{
sProj1 = s1 - dot(s1,rotAxis) * rotAxis;
sProj2 = s2 - dot(s2,rotAxis) * rotAxis;

dotProd2 = dot(sProj1,sProj2) / (norm(sProj1) * norm(sProj2));
crossProd2 = cross(sProj1,sProj2) / (norm(sProj1) * norm(sProj2));
crossSign2 = sign(dot(crossProd2,rotAxis)) * norm(crossProd2);
angCross2 = atan2(crossSign2,dotProd2);

% plot


figure
hold on
plot3([0 uProj1(1)],[0 uProj1(2)],[0 uProj1(3)],'b--')
plot3([0 uProj2(1)],[0 uProj2(2)],[0 uProj2(3)],'b--')

plot3([0 u1(1)],[0 u1(2)],[0 u1(3)],'b')
plot3([0 u2(1)],[0 u2(2)],[0 u2(3)],'b')

plot3([0 dU(1)] + u1(1),[0 dU(2)] + u1(2),[0 dU(3)] + u1(3),'b:')

plot3([0 sProj1(1)],[0 sProj1(2)],[0 sProj1(3)],'r--')
plot3([0 sProj2(1)],[0 sProj2(2)],[0 sProj2(3)],'r--')

plot3([0 s1(1)],[0 s1(2)],[0 s1(3)],'r')
plot3([0 s2(1)],[0 s2(2)],[0 s2(3)],'r')

plot3([0 dS(1)] + s1(1),[0 dS(2)] + s1(2),[0 dS(3)] + s1(3),'r:')

plot3([0 rotAxis(1)],[0 rotAxis(2)],[0 rotAxis(3)],'k')
hold off
grid on
%}
% compute quaternion
sAng = sin(angCross / 2);
cAng = cos(angCross / 2);

q = zeros(4,1);
q(1:3) = sAng * rotAxis;
q(4) = cAng;

q = q / norm(q);

end

