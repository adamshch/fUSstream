function [flowEMD, distEMD] = greedyEMD1D(x1, x2, v1, v2)

% [flowEMD, distEMD] = greedyEMD1D(x1, x2, v1, v2)
%
% Greedy implementation of a 1D EMD algorithm to match points on v1 with 
% values x1 to points in v2 with values x2.
%
% Outputs:
%  - flowEMD  - the flow from [v1,x1] to [v2,x2]
%  - distEMD  - the total distance traveled by all mass (EMD cost)
%
% 2020 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Error/input checking

if numel(v1)~=numel(v2)
    error('location vectors must be the same length!')                     % Make sure vectors are the same distance!
end

if abs(sum(x1)-sum(x2)) > 100*eps
    error('data values must be balanced!')                                 % Make sure the problem is balanced (greedy doesn't work for unbalanced)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up computations

[v1, ix1] = sort(v1,'ascend');                                             % Sort data in ascending order
[v2, ix2] = sort(v2,'ascend');                                             % Sort data in descending order

x1 = x1(ix1);                                                              % Reorganize the points in x1 to match v1
x2 = x2(ix2);                                                              % Reorganize the points in x2 to match v2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute EMD in a greedy way

flowEMD = zeros(numel(x1));                                                % Initialize the flow

for ll = 1:numel(x1)                                                       % Loop through elements in [v1,x1]
    kk = 1;                                                                % Initialize a counter
    while (x1(ll)>eps)&&(kk<=numel(x2))                                    % Proceed until x1(ll) runs out or there's nowhere else to put stuff
        if x1(ll) > x2(kk)                                                 % If x1(ll) is too much, put as much into x2(kk) and move on to kk+1
            flowEMD(ll,kk) = x2(kk);
            x2(kk) = 0;
            x1(ll) = x1(ll) - x2(kk);
        else                                                               % Otherwise put all of x1(ll) into x2(kk) and move on to ll+1
            flowEMD(ll,kk) = x1(ll);
            x2(kk) = x2(kk) - x1(ll);
            x1(ll) = 0;
        end
        kk = kk+1;
    end
end

distMat = abs(bsxfun(@minus, v1(:), v2(:)'));                              % Compute the distance matrix between points
distEMD = sum(distMat(:).*flowEMD(:));                                     % Compute the EMD distance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reorganize flow matrix back

TMP          = 1:numel(x1);                                                % Enumerate the indices of x1
flowReOrder1 = TMP(ix1);                                                   % Get reordering of x1 using the original reordering
flowReOrder2 = TMP(ix2);                                                   % Same for x2
flowEMD      = flowEMD(flowReOrder1, flowReOrder2);                        % Reorder the flows to match the original data & return it


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
