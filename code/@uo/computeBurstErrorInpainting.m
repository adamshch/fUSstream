function [frameInPt,frameIDX] = computeBurstErrorInpainting(uo,varargin)

% [frameInPt,frameIDX] = computeBurstErrorInpainting(fuObj,varargin)
%
%
%
% 2020 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

p = inputParser;                                                           % Set up an object to parse all of the various inputs
p.addParameter('fillOpt' , 'spline');                                      % Select how to fill in the missing frames
parse(p,varargin{:});
p = p.Results;

if (uo.meta.loadMov)&&(~isempty(uo.movie))
    getOnePix = @(z1,z2) uo.movie(z1,z2,:);                                % Extract the current trace
else
    tmpMov    = uo.meta.movMatFile.(uo.meta.movVarName);
    getOnePix = @(z1,z2) tmpMov(z1,z2,:);                                  % Extract the current trace
end

nPix      = prod(uo.frameSize);
errFrames = extractBurstErrs(uo);                                          % Get burst error frames
frameIDX  = find(errFrames==1);                                            % Get the index locations of frame errors
frameInPt = zeros(uo.frameSize(1),uo.frameSize(2),numel(frameIDX));        % Initialize the array of inpainted frames
fprintf('Running...')
for ll = 1:nPix
    [iX,iY] = ind2sub(uo.frameSize,ll);                                    % Get the X-Y coordinates of a given trace
    traceNow           = getOnePix(iX,iY);                                 % Extract the current trace
    traceNow(frameIDX) = NaN;                                              % Fill in the burst error times with NaNs
    TMP                = fillmissing(traceNow,p.fillOpt);                  % Use piecewise cubic spline interpolation to fill in the NaNs
    frameInPt(iX,iY,:) = TMP(frameIDX);                                    % Get the locations that correspond to the burst frame
end
fprintf('done.\n')

clear tmpMov

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function to check if burst errors were isolated and return the frame IDs

function errFrames = extractBurstErrs(uo)

if ~(isfield(uo.errs,'burst')&&(numel(uo.errs.burst)==uo.movieLen))
    warning('No burst error detection output found. Running findBurstFrames with default settings...\n')
    uo.findBurstFrames(); 
end
errFrames = uo.errs.burst;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%