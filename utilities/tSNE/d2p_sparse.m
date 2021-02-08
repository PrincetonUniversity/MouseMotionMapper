function [P, beta] = d2p_sparse(D, u, tol,maxNeighbors)
%D2P Identifies appropriate sigma's to get kk NNs up to some tolerance 
%
% (C) Laurens van der Maaten, 2008
% Maastricht University
%
% Modifications by Gordon J. Berman, 2014
% Princeton University

    if nargin < 4 || isempty(maxNeighbors)
        maxNeighbors = 150;
    end

    
    if ~exist('u', 'var') || isempty(u)
        u = 15;
    end
    if ~exist('tol', 'var') || isempty(tol)
        tol = 1e-4;
    end
    
    % Initialize some variables
    n = size(D, 1);                     % number of instances
    if maxNeighbors >= n
        maxNeighbors = n - 1;
    end
    
    beta = ones(n, 1);                  % empty precision vector
    logU = log(u);                      % log of perplexity (= entropy)
    
    
    jj = zeros(n,maxNeighbors);
    vals = zeros(size(jj));
    
    % Run over all datapoints
    for i=1:n
        
        if ~rem(i, 500)
            disp(['Computed P-values ' num2str(i) ' of ' num2str(n) ' datapoints...']);
        end
        
        % Set minimum and maximum values for precision
        betamin = -Inf;
        betamax = Inf;
        
        q = D(i,:);
        [sortVals,sortIdx] = sort(q,'ascend');
        sortVals = sortVals(2:(maxNeighbors+1));
        sortIdx = sortIdx(2:(maxNeighbors+1));
        jj(i,:) = sortIdx;
        
        % Compute the Gaussian kernel and entropy for the current precision
        [H, thisP] = Hbeta(sortVals, beta(i));
        
        % Evaluate whether the perplexity is within tolerance
        Hdiff = H - logU;
        tries = 0;
        while abs(Hdiff) > tol && tries < 50
            
            % If not, increase or decrease precision
            if Hdiff > 0
                betamin = beta(i);
                if isinf(betamax)
                    beta(i) = beta(i) * 2;
                else
                    beta(i) = (beta(i) + betamax) / 2;
                end
            else
                betamax = beta(i);
                if isinf(betamin)
                    beta(i) = beta(i) / 2;
                else
                    beta(i) = (beta(i) + betamin) / 2;
                end
            end
            
            % Recompute the values
            [H, thisP] = Hbeta(sortVals, beta(i));
            Hdiff = H - logU;
            tries = tries + 1;
        end
        

        vals(i,:) = thisP;

    end
    
    
    ii = repmat((1:n)',1,maxNeighbors);
    ii = reshape(ii',[n*maxNeighbors 1]);
    jj = reshape(jj',[n*maxNeighbors 1]);
    vals = reshape(vals',[n*maxNeighbors 1]);
    
    
    P = sparse(ii,jj,vals,n,n,n*maxNeighbors);
    
    clear ii jj vals
    
    disp(['Mean value of sigma: ' num2str(mean(sqrt(1 ./ beta)))]);
    disp(['Minimum value of sigma: ' num2str(min(sqrt(1 ./ beta)))]);
    disp(['Maximum value of sigma: ' num2str(max(sqrt(1 ./ beta)))]);
end



% Function that computes the Gaussian kernel values given a vector of
% squared Euclidean distances, and the precision of the Gaussian kernel.
% The function also computes the perplexity of the distribution.
function [H, P] = Hbeta(D, beta)
    P = exp(-D * beta);
    sumP = sum(P);
    H = log(sumP) + beta * sum(D .* P) / sumP;
    P = P / sumP;
end

