function [T,densities,stateDefs,states,numTransitions] = findTransitionMatrix(states,eliminateZeros,eliminateSelfTransitions)

if nargin < 3 || isempty(eliminateSelfTransitions)
    eliminateSelfTransitions = true;
end

if nargin < 2 || isempty(eliminateZeros)
    eliminateZeros = true;
end


stateDefs = sort(unique(states));
if eliminateZeros
    stateDefs = setdiff(stateDefs,0);
    states = states(states ~= 0);
end
L = length(stateDefs);


if eliminateSelfTransitions
    a = [true; abs(diff(states)) > 0];
    states = states(a);
end


T = zeros(L);
densities = zeros(L,1);
N = length(states);
numTransitions = 0;
for i=1:L
    idx = setdiff(find(states==stateDefs(i)),N);
    densities(i) = length(idx);
    vals = states(idx + 1);
    
    numTransitions = numTransitions + length(vals);
    for j=1:L
        T(i,j) = sum(vals == stateDefs(j)) ./ length(vals);
    end
    
end

densities = densities ./ sum(densities);



