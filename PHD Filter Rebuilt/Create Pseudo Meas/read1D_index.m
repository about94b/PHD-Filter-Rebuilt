function [ output ] = read1D_index(fileName,ds_name,indices)
%READ1D_INDEX Get the entries of a 1D vector stored in a .h5 file in the
%database ds_name at locations given by indices.
%   [ output ] = read1D_index(fileName,ds_name,indices)

% initialize output
nEntries = length(indices);
output = zeros(nEntries,1);

% get the series of consecutive numbers in indices
indexStep = diff(indices);
sequenceEnds = find(indexStep > 1)';
if isempty(sequenceEnds)
    sequenceStartInd = 1:nEntries;
    sequenceStart = indices';
    sequenceLength = ones(1,nEntries);
else
    sequenceStartInd = [1 sequenceEnds + 1];
    sequenceStart = indices(sequenceStartInd)';
    sequenceLength = [sequenceEnds(1) diff(sequenceEnds) (nEntries - sequenceEnds(end))];
end

% read data
nSeq = length(sequenceStart);
for i = 1:nSeq
    
    indStart = sequenceStartInd(i);
    indEnd = indStart + sequenceLength(i) - 1;
    indList = (indStart:indEnd)';
    
    output(indList) = h5read(fileName,ds_name,sequenceStart(i),sequenceLength(i));
    
end

end

