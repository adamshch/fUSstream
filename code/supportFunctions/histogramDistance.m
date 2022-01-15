function distVal = histogramDistance(data1, data2, varargin)

% distVal = histogramDistance(data1, data2, ['distSelect', distSelect, 'binNumber', binNumber])
%
% Function to compute the distance between two histograms over 'data1' 
% and 'data2'. Histograms are binned into 'binNumber' number of bins and
% the distance is selected by 'distSelect' (default 'emd')
%
% 2021 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

p = inputParser;                                                           % Set up an object to parse all of the various inputs
p.addParameter('distSelect'    , 'emd'   );                                % Select the metric to use for the distance
p.addParameter('binNumber'     , 50      );                                % Select the bin numbers for the histogram
parse(p,varargin{:});
p = p.Results;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create histograms

maxVal = max(max(data1), max(data2));                                      % Compute the max value over all data
data1  = data1/maxVal;                                                     % Normalize by the max value for sanity & numerical stability
data2  = data2/maxVal;                                                     % Ditto for dataset 2

[counts1,edges1] = histcounts(data1,p.binNumber);                          % Compute the histogram counts for the first dataset
[counts2,edges2] = histcounts(data2,p.binNumber);                          % Compute the histogram counts for the second dataset
counts1          = counts1(:)/sum(counts1);                                % Normalize by number of counts
counts2          = counts2(:)/sum(counts2);                                % Normalize by number of counts
centers1         = 0.5*(edges1(2:end) + edges1(1:(end-1)));                % Compute the centers for the first dataset
centers2         = 0.5*(edges2(2:end) + edges2(1:(end-1)));                % Compute the centers for the second dataset

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute the distance

switch lower(p.distSelect)                                                 % Look at which distance was selected
    case 'emd'                                   
        [~, distVal] = greedyEMD1D(counts1, counts2, centers1, centers2);  % If EMD, use the greedy emd solver
    case 'kl'
        distVal = sum(counts1.*log(counts1./(counts2+eps)));               % KL can be computed in closed form. eps included for stability
    case 'symkl' 
        KL12    = sum(counts1.*log(counts1./(counts2+eps)));               % Symmetric KL computes both directions...
        KL21    = sum(counts1.*log(counts1./(counts2+eps)));               % ...
        distVal = sqrt(KL21.*KL12);                                        % and then returns the geometric average b/w the two
    case 'l2'
        distVal = sum((counts1 - counts2).^2);                             % L2 distance is what it is... don't use this please.
    otherwise
        p.distSelect = 'emd';                                              % Default is EMD
        distVal = histogramDistance(data1, data2, 'distSelect', ...
                                  p.distSelect, 'binNumber', p.binNumber); % Run with default distance by having the function call itself
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
