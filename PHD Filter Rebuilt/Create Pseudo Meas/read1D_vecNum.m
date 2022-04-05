function [ output ] = read1D_vecNum(fileName,ds_name,indices,numMatch)
%READ1D_INDEX Get the entries of a 1D array stored in a .h5 file in the
%database ds_name. The vector numbers are indices and the length of each 
%vector in the array is numMatch.
%   [ output ] = read1D_vecNum(fileName,ds_name,indice,numMatch)

% get indices to read out
indList = [];
if size(indices,1) > 1
    indices = indices';
end
for i = indices
    if length(numMatch) == 1
        indStart = numMatch * (i - 1) + 1;
        indEnd = numMatch * i;
        indList = cat(1,indList,(indStart:indEnd)');
    else

        indStart = sum(numMatch(1:(i - 1))) + 1;
        indEnd = sum(numMatch(1:i));
        indList = cat(1,indList,(indStart:indEnd)');

    end
end

output = read1D_index(fileName,ds_name,indList);

end

