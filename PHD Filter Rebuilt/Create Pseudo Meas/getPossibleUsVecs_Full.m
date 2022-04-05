function [ vecList ] = getPossibleUsVecs_Full( delta,fileName,measList,dist2obs,boundsList )
%GETPOSSIBLEUSVECS Summary of this function goes here
%   Detailed explanation goes here
%   [ vecList ] = getPossibleUsVecs_Full( delta,fileName,measList,dist2obs,boundsList )
%
%   vecList column format (each contains two vectors)
%       row    1: gamma value of the pair
%       rows 2-4: uHat vec
%       rows 5-7: sHat vec

% database names
ds_sveclist = '/sVecList';
ds_uveclist = '/uVecList';
ds_gammalist = '/gammaList';%'/ds_gammalist';
ds_nummatch = '/uVecNumMatch';

uVecNumMatch = h5read(fileName,ds_nummatch);
uVecList = h5read(fileName,ds_uveclist);
nUvec = length(uVecNumMatch);

cD = cos(delta);


% visibility groups to check
% possibleGroups1 = getPossibleGroups(fileName,meas+bounds,dist2obs);
% possibleGroups2 = getPossibleGroups(fileName,max(meas-bounds,0),dist2obs);
% possibleGroups = unique(cat(1,possibleGroups1,possibleGroups2));

% boundInt = bounds / (dist2obs^2 * pi);

% compute gamma of the measurement
gammaMeasArray = pi * dist2obs.^2 .* measList;

% iterate through uVecs
nMeas = length(measList);
vecList = cell(nMeas,1);
% sVecNum = zeros(nMeas,1);
nPairOld = zeros(nMeas,1);
tic
for i = 1:nUvec
    
    if i == 1 || i == nUvec || mod(i,100) == 0
        toc
        fprintf('\t%d/%d\n',i,nUvec)
        tic
    end
    
    % load the sVec and uVec data
%     uVec = read1D_vecNum(fileName,ds_uveclist,i,3);
    uVec = uVecList((3 * (i - 1) + 1):(3 * i));
    nSmatch = uVecNumMatch(i,1);
    sVecNum = uVecNumMatch(i,2);
    indList = ((sVecNum + 1):(sVecNum + nSmatch))';
    sVecList = read1D_vecNum(fileName,ds_sveclist,indList,3);
    
    % get difference between each sVec's gamma and measured gamma
    gammaList = read1D_vecNum(fileName,ds_gammalist,i,uVecNumMatch);
    
    
    % indices of sVecs where gamma is within bounds of measured
    for k = 1:nMeas
        
        gammaDiff = gammaList - gammaMeasArray(k);
        sIndList = find(abs(gammaDiff) <= boundsList(k));

        % store gammas
        nPairNew = length(sIndList);
        vecList{k} = cat(2,vecList{k},zeros(7,nPairNew));
        columnInds = (nPairOld(k) + 1):(nPairOld(k) + nPairNew);
        vecList{k}(1,columnInds) = gammaList(sIndList)';

        % store sVecs and uVec
        for j = 1:nPairNew

            % get sVec
            indStart = 3 * (sIndList(j) - 1) + 1;
            indEnd = 3 * sIndList(j);
            sVec = sVecList(indStart:indEnd);

            if abs(dot(uVec,sVec) - cD) > 1e-8
                abs(dot(uVec,sVec) - cD)
            end

            % store
            vecList{k}(2:4,j + nPairOld(k)) = uVec;
            vecList{k}(5:7,j + nPairOld(k)) = sVec;

        end

        % update number of pairs
        nPairOld(k) = nPairOld(k) + nPairNew;
        
    end
    
    % update number of sVecs
%     sVecNum = sVecNum + nSmatch;
    
end

% get the best match
% if ~isempty(vecList)
%     bestMatch = zeros(3,2);
%     [~,bestMatchInd] = min(abs(vecList(1,:) - gamma));
%     bestMatch(:,1) = vecList(2:4,bestMatchInd);
%     bestMatch(:,2) = vecList(5:7,bestMatchInd);
%     gammaError = vecList(1,bestMatchInd) - gamma; 
% else
%     bestMatch = [];
%     gammaError = [];
% end

end

