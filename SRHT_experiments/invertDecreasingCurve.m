function [indices, lasttarget] = ...
             invertDecreasingCurve(valuematrix, targets, p)
% [indices, lasttarget] = invertDecreasingCurve(valuematrix, targets, p)
%
% Given a matrix valuematrix such that the entries in each row
% decrease with increasing column index, a probability p, and a list of 
% target values in targets, returns a vector of the first column indices 
% into valuematrix so that at most p percentage of the entries in that 
% column are larger than the corresponding target value
%
% It may not be possible to find such a column index for each targetvalue.
% lasttarget is the index of the last target that is achievable

numrows = size(valuematrix, 1);
indices = zeros(1,size(valuematrix,2));
lasttarget = length(targets);

for i=1:length(targets)
    comparisonvec = sum(valuematrix <= targets(i))/numrows;
    if all(comparisonvec <= p)
        lasttarget = i-1;
        return;
    end
    indices(i) = find( comparisonvec > (1-p), 1);
end

end