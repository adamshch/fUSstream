function varargout = displayExampleTraces(uo, varargin)

% displayExampleTraces(fuObj, varargin)
%
% Function to display example traces from a functional ultrasound movie
%
% 2019 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

p = inputParser;                                                           % Set up an object to parse all of the various inputs
p.addParameter('baseSelect'    , 'median');                                % Select what criteria to calculate the baseline image as
p.addParameter('basePrctile'   , 25      );                                % For the situation that 'prctile' is chosen, which percentile should the baseline image be calculated using
p.addParameter('pixSelect'     , 'vary'  );                                % Select which pixels to select: either an array or 'vary' or 'rand'
p.addParameter('numTraces'     , []      );                                % Select hoe many traces to select
p.addParameter('normTraces'    , false   );                                % Choose if each trace is normalized to its own max or not
p.addParameter('dispBurstErrs' , false   );                                % Select underlay of burst error locations
p.addParameter('dispMotErrs'   , false   );                                % Select underlay of motion correction errors
p.addParameter('smoothLvl'     , 0       );                                % Select denoising level
p.addParameter('motionCorrect' , false   );                                % Select if the motion corrected data should be used
p.addParameter('figNo'         , 105     );                                % Select figure number
p.addParameter('maskOnly'      , true    );                                % Select if only time-traces within the mask should be displayed

parse(p,varargin{:});
p = p.Results;

getTrace    = uo.makeGetTraceFunction(p);                                  % Make the anonymous function that extracts single traces
p.pixSelect = selectFusiTraces(uo,getTrace, p.numTraces,...
                                               p.pixSelect, p.maskOnly);   % Select the pixels to plot
movSel      = extractPixel(uo, getTrace, p.pixSelect);                     % Extract the appropriate pixel time-traces to plot

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate the baseline

switch p.baseSelect
    case 'median'                                                          % If the median is selected...
        baseVals = nanmedian(movSel,1);                                       %  ...Calculate the median of each pixel
    case 'mode'                                                            % If the mode is selected...
        baseVals = halfSampleMode(movSel);                                 %  ...Calculate the half-sample mode for each pixel in the image
    case 'mean'                                                            % If the mean is selected...
        baseVals = nanmean(movSel,1);                                         %  ...Calculate the mean of each pixel in the image
    case 'prctile'                                                         % If a percentile is selected
        baseVals = prctile(movSel,p.basePrctile,1);                        %  ...Calculate the percentile of each pixel's histogram1
    otherwise                                                              % The default is no correction so...
        baseVals = 0;                                                      %  ...Set the baseline image to nothing
end
baseVals(isnan(baseVals)) = 1;                                             % In the chance that a bad pixel is chosen, pick an inoffensive value

movSel = normalizeTraces(movSel, baseVals, p.normTraces);                  % subtract the baseline and normalize the traces

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check for error frames

errFrames = extractBurstErrs(uo);                                          % Get burst error frames
movSel    = denoiseTraces(movSel, p.smoothLvl, errFrames);                 % Optional denoising step
[motDep, motLat, maxMot] = getMotionVectorsWrapper(uo, p);                 % Check for motion errors if needed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate labels and other plotting needs

[tt, ttLab, ylab, pxLab, figTitle] = genLabeling(uo,p.pixSelect,...
                                                            p.normTraces); % Generate labels for the plot
[vertLims, ytk, pxLab] = setVertLims(movSel, p.dispMotErrs, maxMot, pxLab);% Set the vertical limits


if p.dispBurstErrs
    burstLocs              = nan(size(movSel,1),1);
    burstLocs(errFrames,:) = size(movSel,2)+4;
else
    burstLocs = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Display the time series


figure(p.figNo); clf;
figure(p.figNo), hold on;
legendNames = {};
if p.dispBurstErrs
    stem(tt, burstLocs, 'k.');                                              % Plot burst error if requested
    legendNames{1} = 'Burst errors';
end

if p.dispMotErrs
    totalMotErr = sqrt(motDep(:).^2+motLat(:).^2);
    h = fill([tt(:); flipud(tt(:))]',[0*tt(:)+1+size(movSel,2);...
                                flipud(motDep(:)+1+size(movSel,2))]','r');
    set(h,'facealpha',0.3)
    legendNames{2} = 'Depth motion offset';
    h = fill([tt(:); flipud(tt(:))]',[0*tt(:)+1+size(movSel,2);...
                                flipud(motLat(:)+1+size(movSel,2))]','b');
    set(h,'facealpha',0.3)
    legendNames{3} = 'AP motion offset';
    plot(tt,size(movSel,2)+1+totalMotErr,'k')
    legendNames{3} = 'Total motion offset';
end
plot(tt, bsxfun(@plus, movSel, 0:(size(movSel,2)-1)))                      % Plot the actual traces

ylabel(ylab);                                                              % Set the Y-label
xlabel(ttLab);                                                             % Set the X-label
set(gca, 'XLim', [min(tt), max(tt)], 'YTick', ytk, 'YTickLabel', pxLab,...
                                                         'YLim', vertLims);% Set some plotting options
box off; colormap gray
title(figTitle)
figure(p.figNo), hold off
legend(legendNames)
set(gcf, 'color', [1,1,1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output parsing

if nargout > 0
    varargout{1} = p.pixSelect;
end

if nargout > 0
    varargout{2} = movSel;
end 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extra Functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function movSel = normalizeTraces(movSel, baseVals, normTraces)

normMod = 0.8;
movSel  = bsxfun(@plus, movSel, -baseVals);                                % Remove baseline values
if normTraces
    movSel = movSel*diag(1./max(normMod*movSel,[],1));                     % If requested, normalize each trace to its own max
else
    movSel = movSel./(normMod*max(movSel(:)));                             % Otherwise scale everything proportionally
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [tt, ttLab, ylab, pxLab, figTitle] = genLabeling(uo,pixSelect,normTraces)

if isempty(uo.dt)
    tt = 1:uo.movieLen;
    ttLab = 'Frame number';
else 
    tt = 0:uo.dt:uo.dt*(uo.movieLen-1);
    ttLab = 'Time (s)';
end

pxLab = cell(1,size(pixSelect,1));
for ll = 1:size(pixSelect,1)
    pxLab{ll} = sprintf('(%d,%d)', pixSelect(ll,1), pixSelect(ll,2));
end

if normTraces
    figTitle = 'Example normalized traces';
    ylab     = 'CBV (normalized)';
else
    figTitle = 'Example time traces';
    ylab     = 'CBV';
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function to automatically set vertical limits of plot

function [vertLims, ytk, pxLab] = setVertLims(movSel, dispMotErrs, maxMot, pxLab)

minVal = min(movSel(:,1));
if minVal<0; minVal = minVal*1.2;
else       ; minVal = 0.8*minVal;
end


maxVal   = size(movSel,2)-0.5;
vertLims = [minVal, maxVal];

ytk = 0:(size(movSel,2)-1);
if dispMotErrs
    ytk = cat(2,ytk,size(movSel,2)+1+[-0.8,0,0.8]);
    vertLims = vertLims + [0,3];                                           % Extend upper limit by one
    pxLab = cat(2,pxLab,{[num2str(-maxMot,3),'px'], 'Motion', [num2str(maxMot,3),'px']});
else
    vertLims = vertLims + [0,0.5];  
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function to check and extract motion errors

function [motDep, motLat, maxMot] = getMotionVectorsWrapper(uo, p)

if p.dispMotErrs
   [motDep, motLat, maxMot] = getMotionVectors(uo);
    motDep = 0.8*motDep/maxMot;
    motLat = 0.8*motLat/maxMot;
else
    motLat = [];
    motDep = [];
    maxMot = NaN;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function to denoise time traces

function movSel = denoiseTraces(movSel, smoothLvl, errFrames)

iX = sum(isnan(movSel))>0.25*size(movSel,1)*ones(1,size(movSel,2));
movSel(:, iX) = 0;

movSel(errFrames,:) = NaN;                                                 % Designate burst errors as missing data
for ll = 1:size(movSel,2)
    movSel(:,ll) = fillmissing(movSel(:,ll),'spline');                     % Use piecewise cubic spline interpolation to fill in the NaNs
end

if smoothLvl == 0                                                          % No smoothing, nothing to do
elseif smoothLvl > 0
    movSel = wdenoise(double(movSel), smoothLvl, 'DenoisingMethod', ...
                                             'BlockJS', 'Wavelet', 'sym4');% Denoise the traces
else
    warning('Bad smoothLvl value: defaulting to NOT smoothing time traces')
end

movSel(:, iX) = NaN;

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

