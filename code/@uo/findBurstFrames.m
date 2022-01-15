function varargout = findBurstFrames(uo, varargin)

% errIDs = findBurstFrames(fData, varargin)
%
% 
% 2019 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

if nargin > 1; thresh = varargin{1};
else         ; thresh = [];
end

if nargin > 2; outlierSel = varargin{2};
else         ; outlierSel = 'direct';
end

if nargin > 3; Nthresh = varargin{3};
else         ; Nthresh = 100;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

if (uo.meta.loadMov)&&(~isempty(uo.movie))
    frameNorms = sqrt(sum(sum(uo.movie.^2 , 1), 2));                       % Get the norms of each frame
else
    frameNorms = sqrt(sum(sum(uo.meta.movMatFile.(uo.meta.movVarName).^2,1),2));% Get the norms of each frame
end
frameNorms = frameNorms(:);                                                % Reorient into a row vector;

if isempty(thresh)                                                         % If no threshold is given...
    thresh = [min(frameNorms), max(frameNorms)];                           %  ... check the entire dynamic range from minimum to maximum values
end

switch outlierSel
    case 'direct'
        errIDs = directOutlierThresholding(frameNorms, thresh, Nthresh);
    case 'robustSTD'
        errIDs = robustSTDOutlierThresholding(frameNorms, thresh, Nthresh);
    otherwise
        error('unknwon outlier selection type. options are direct or robustSTD')
end

uo.errs.burst = errIDs;
if nargout > 0
    varargout{1} = fuObj.errs.burst;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extra functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function errIDs = directOutlierThresholding(frameNorms, thresh, Nthresh)
if isscalar(thresh)                                                        % If only one threshold is given...
    errIDs = frameNorms > thresh;                                          %  ... do a threshold detection.
elseif numel(thresh) >=2                                                   % Otherwise check which threshold gives the best detection
    if numel(thresh)==2                                                    % 
        thresh = linspace(min(thresh), max(thresh), Nthresh);              % If only 2 values are given, check a range of thresholds between them
    end
    minThresh = min(thresh);                                               % Extract the minimum value of the thresholds
    errLevel  = zeros(size(thresh));                                       % Initialize the error summing vector
    for kk = 1:numel(thresh)
        errLevel(kk) = sum(frameNorms > thresh(kk));                       % For each threshold level, check the number of points that fail
    end
    
    thresh     = thresh(errLevel>0);                                       % Remove the zero-error points (no outliers at all) - both the thresholds...
    errLevel   = errLevel(errLevel>0);                                     %  ... and the values
    errLvlRate = diff(errLevel);                                           % Calculate the rate of change in outliers w.r.t. the threshold
    
    if numel(errLevel) <=1                                                 % If the error levels were all zero...
        warning('Thresholds chosen are too high!')                         %  ... Give a warning...
        errIDs = frameNorms > minThresh;                                   %  ... And just use the minimum threshold
    else
        lowHist = find(errLvlRate == min(-errLvlRate), 1, 'first');        % Otherwise, find the first major dip in the rate of change of the outliers w.r.t. the threshold
        errIDs = frameNorms > thresh(lowHist+1);                           % Use that point as the outlier threshold
    end
end
end

function errIDs = robustSTDOutlierThresholding(frameNorms, thresh, Nthresh)

[normSTD, ~, normDist] = robustcov(frameNorms);                            % Calculate the robust standard deviation and distances of each point from the robust mean

if isscalar(thresh)                                                        % If only one threshold is given...
    errIDs = normDist > thresh*normSTD;                                    %  ... do a threshold detection based on being thresh robust STDs away from the robust mean
elseif numel(thresh) >=2   
    if numel(thresh)==2                                                    % 
        thresh = linspace(min(thresh), max(thresh), Nthresh);              % If only 2 values are given, check a range of thresholds between them
    end
    minThresh = min(thresh);                                               % Extract the minimum value of the thresholds
    errLevel  = zeros(size(thresh));                                       % Initialize the error summing vector
    for kk = 1:numel(thresh)
        errLevel(kk) = sum(normDist > thresh(kk)*normSTD);                 % For each threshold level, check the number of points that fail
    end
    
    thresh     = thresh(errLevel>0);                                       % Remove the zero-error points (no outliers at all) - both the thresholds...
    errLevel   = errLevel(errLevel>0);                                     %  ... and the values
    errLvlRate = diff(errLevel);                                           % Calculate the rate of change in outliers w.r.t. the threshold
    
    if numel(errLevel) <=1                                                 % If the error levels were all zero...
        warning('Thresholds chosen are too high!')                         %  ... Give a warning...
        errIDs = frameNorms > minThresh;                                   %  ... And just use the minimum threshold
    else
        lowHist = find(errLvlRate == min(-errLvlRate), 1, 'first');        % Otherwise, find the first major dip in the rate of change of the outliers w.r.t. the threshold
        errIDs = frameNorms > thresh(lowHist+1);                           % Use that point as the outlier threshold
    end
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
