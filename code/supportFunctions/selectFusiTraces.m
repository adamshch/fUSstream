function pixSelect = selectFusiTraces(uo,getTrace,numTraces,pixSelect,varargin)

% pixSelect = selectFusiTraces(uo,getTrace,numTraces,pixSelect, [useMotionCorrSize, maskOnly])
%
% Select which pixel's time traces to plot. The inputs are:
%  - uo        - ultrasound object (see @uo folder)
%  - getTrace  - anonymous function to extract a single time-trace from the
%                   up data
%  - numTraces - scalar number of time-traces to extract (default 10)
%  - pixSelect - how the pixels should ve selected. Options include 
%                    'rand' - completely random selection
%                    'vary' - random slection but with traces that selected
%                             to be less correlated than 'rand'
%
%  Optional inputs: 
%  - useMotionCorrSize - T/F (default false) flag to use the motion corrected
%                           indexing. 
%  - maskOnly          - T/F (default true) flag to select traces only 
%                           within the user-specified mask.
%
% Returns pixSelect: a set of pixels sub-selected from the full data
%
% 2020 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

if isempty(pixSelect);  pixSelect = 'vary';  end
if isempty(numTraces);  numTraces = 10    ;  end

if nargin > 4;     useMotionCorrSize = varargin{1};
else;              useMotionCorrSize = false;
end

if nargin > 5;     maskOnly = varargin{2};
else;              maskOnly = true;
end

if useMotionCorrSize
    frameSize = uo.getMotionCorrectedSize();
else
    frameSize = uo.frameSize;
end

if maskOnly;    uo.ensureMask();   end                                     % Make sure that a mask is available

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Iterate through and select traces

if ~isnumeric(pixSelect)
    switch pixSelect
        case 'vary'
            [a,u] = findSigmoidPars([0.5,0.2], [0.8,0.9]);                 % Set sigmoidal function that will determine the probability of picking a trace
            TMPtraces = zeros(uo.movieLen, numTraces);                     % Initialize temporary storage of time traces
            pixSelect = zeros(numTraces,2);                                % Initialize the pixel selection
            for ll = 1:1:numTraces                                         % Loop to select time traces
                selTrace = false;                                          % Initialize correlation test
                while ~selTrace
                    pixTMP   = randsample(prod(uo.frameSize), 1);          % Sample pixels totally at random
                    [I,J]    = ind2sub(frameSize, pixTMP);                 % Translate the pixel value to a 2D index
                    tSel     = squeeze(getTrace(I,J));                     % Pull out and reshape the time-trace
                    corMat   = corrcoef([tSel(:), TMPtraces(:,1:(ll-1))]); % Get the correlation coefficients for all traces so far
                    corVal   = max(corMat(2:end,1));                       % Extract only the correlation coefficients involving the new trace
                    prSel    = 1./(1+exp(-a*(abs(corVal)-u)));             % Probability that the trace is selected
                    if any(isnan(tSel)); prSel = 0; end                    % In the case where NaNs exist, don't pick that trace (i.e., the field-of-view edges after motion selection
                    if maskOnly&&(uo.mask(I,J)==0); prSel = 0; end         % In the case where NaNs exist, don't pick that trace (i.e., the field-of-view edges after motion selection
                    selTrace = rand(1)<prSel;                              % Select with probability prSel
                end
                TMPtraces(:,ll) = tSel;                                    % Store the selected time trace...
                pixSelect(ll,:) = [I,J];                                   %   ... and the pixel selection
            end
        case 'rand'
            pixSelect = zeros(numTraces,2);                                % Initialize the pixel selection
            for ll = 1:1:numTraces                                         % Loop to select time traces
                selTrace = false;                                          % Initialize correlation test
                while ~selTrace
                    pixTMP   = randsample(prod(uo.frameSize), 1);          % Sample pixels totally at random
                    [I,J]    = ind2sub(frameSize, pixTMP);                 % Translate the pixel value to a 2D index
                    tSel     = squeeze(getTrace(I,J));                     % Pull out and reshape the time-trace
                    selTrace = ~any(isnan(tSel));                          % Select if no NaNs
                end
                pixSelect(ll,:) = [I,J];                                   %   ... and the pixel selection
            end
        otherwise
            warning('Invalid pixel selection: Defaulting to random.\n')
            pixSelect = selectFusiTraces(uo,numTraces,'rand');
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
