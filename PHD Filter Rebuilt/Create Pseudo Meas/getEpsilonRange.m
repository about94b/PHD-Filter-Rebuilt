function [ epsilonU ] = getEpsilonRange( delta,D )
%GETEPSILONRANGE Uses the method in test_epsilonU_curves to compute the
%range of possible epsilon values between uHat and nHat for a reflecting
%facet.
%   [ epsilonU ] = getEpsilonRange( delta,D )
%   epsilonU - array containing min and max epsilonU
%   delta - angle between sHat and uHat
%   D - ratio I / alpha = I pi robs^2 / ACd

% maximum epsilon
epMax = acos(D);

% intersection of uHatEpsilon and sHatEpsilon
commonEpsilon = acos(sqrt(D));

% test if an intersection occurs on epU + epS
sumTestVal = 2 * D - cos(delta);
sumTest = (abs(sumTestVal) < 1);
if sumTest
    epSumU = zeros(2,1);
    epSumU(1) = 0.5 * (delta - acos(2 * D - cos(delta)));
    epSumU(2) = 0.5 * (delta + acos(2 * D - cos(delta)));
    if epSumU(2) > epMax
        epSumU(2) = [];
    end
    if epSumU(1) < 0
        epSumU(1) = [];
    end
    epsumS = acos(D ./ cos(epSumU));
    nSum = length(epsumS);
    delSum = [];
    for j = 1:nSum
        if abs(epsumS(j) - (delta - epSumU(j))) > 1e-8
            delSum = [delSum;j]; %#ok<AGROW>
        end
    end
    epsumS(delSum) = [];
    epSumU(delSum) = [];
    if isempty(epsumS)
        sumTest = 0;
    end
    
end

% test if an intersection occurs on |epU - epS|

% number of regions
if commonEpsilon > epMax
    nDiff = 1;
else
    nDiff = 2;
end

diffTest = (abs(sumTestVal) < 1);
if diffTest
    epDiffU = zeros(2 * nDiff,1);
    epDiffS = zeros(2 * nDiff,1);
    epDiffU(1) = 0.5 * (-delta - acos(2 * D - cos(delta)));
    epDiffU(2) = 0.5 * (-delta + acos(2 * D - cos(delta)));
    delDiff = [];
    for i = 1:2
        if (epDiffU(i) < 0) || (epDiffU(i) > epMax)
            delDiff = [delDiff;i]; %#ok<AGROW>
        elseif nDiff == 2 && (epDiffU(i) > commonEpsilon)
            delDiff = [delDiff;i]; %#ok<AGROW>
        else
            epDiffS(i) = acos(D / cos(epDiffU(i)));
        end
    end
    if nDiff == 2
        epDiffU(3) = 0.5 * (delta - acos(2 * D - cos(delta)));
        epDiffU(4) = 0.5 * (delta + acos(2 * D - cos(delta)));
        for i = 3:4
            if (epDiffU(i) < commonEpsilon) || (epDiffU(i) > epMax)
                delDiff = [delDiff;i]; %#ok<AGROW>
            else
                epDiffS(i) = acos(D / cos(epDiffU(i)));
                if abs(abs(epDiffS(i) - epDiffU(i)) - delta) > 1e-8
                    delDiff = [delDiff;i]; %#ok<AGROW>
                end
            end
        end
    end

    epDiffU(delDiff) = [];
    epDiffS(delDiff) = [];
    
    if isempty(epDiffS)
        diffTest = 0;
    end
end

%% get the range of possible epsilonU

% the difference test first
if diffTest
    if length(epDiffU) == 1
        if epDiffU < commonEpsilon
            goodAfter = 1;
        else
            goodAfter = 0;
        end
        
        if goodAfter
            epsilonU = [epDiffU epMax]';
        else
            epsilonU = [0 epDiffU]';
        end
    else
        epsilonU = epDiffU;
    end
else
    epsilonU = [0 epMax]';
end

% the summation test
if sumTest
    if length(epSumU) == 1
        if epSumU < commonEpsilon
            goodAfter = 1;
        else
            goodAfter = 0;
        end
        
        if goodAfter
            if epsilonU(2) < epSumU
                epsilonU = [];
            elseif epsilonU(1) < epSumU
                epsilonU(1) = epSumU;
            end
        else
            if epsilonU(1) > epSumU
                epsilonU = [];
            elseif epsilonU(2) > epSumU
                epsilonU(2) = epSumU;
            end
        end
    else
        if epsilonU(1) < epSumU(1)
            epsilonU(1) = epSumU(1);
        end
        if epsilonU(2) > epSumU(2)
            epsilonU(2) = epSumU(2);
        end
    end
end




end

