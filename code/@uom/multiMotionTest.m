function uom = multiMotionTest(uom, refDataset, varargin)

% uom = multiMotionTest(uom, varargin)
%
% Function to test and align multiple fUSi datasets
%
% 2020 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

p = inputParser;                                                           % Set up an object to parse all of the various inputs
p.addParameter('mcType'      , 'batchMedian');                             % Option of what motion correction type to use ('median', 'batchMedian'). This list will expand
p.addParameter('frameSel'    , []);                                        % Option to select which frames to motion correct
p.addParameter('refImage'    , []);                                        % Option to input a user-defined reference image (useful for multi-dataset computations)
p.addParameter('medBaseLine' , 'firstmedian');                             % Select the reference image
p.addParameter('batchSz'     , 5);                                         % Select batch sizes to reduce computation time
p.addParameter('useDenoised' , false);                                     % Option to motion correct based on the denoised data
p.addParameter('globalOrient', true);                                      % Option to globally orient the data
p.addParameter('searchBlock' , 5);                                         % Set the limit (in pixels) that the motion correction searches for shifts over
p.addParameter('useMask'     , true);                                      % Option to use the mask in order to compute shifts only based on brain motion (ignore out-of-brain motion)                           
% p.addParameter('refDataset'  , 0);                                         % Option to select which reference dataset to align all sessions to
parse(p,varargin{:});
p = p.Results;

if isempty(refDataset); refDataset = 0; end

plotOpt = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find the dataset with the smallest jitter:

for ll = 1:numel(uom.uo)                                                   % Iterate over all the data objects
    fprintf('Getting motion error for dataset %d of %d:\n', ll, numel(uom.uo))
    uom.uo{ll}.findMotionCorrectionError(varargin{:});                     % Get motion errors for llth dataset, passing through parameters
    fprintf('done.\n')
end

for ll = 1:numel(uom.uo)                                                   % Iterate over all the data objects
    uom.uo{ll}.correctResidualMotion('reCalc', true);                      % Run sub-pixel rigid shifts for motion correction
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Use that dataset to align all blocks in other datasets

allMotTots = zeros(numel(uom.uo),3);                                       % Initialize an array to store the total deviations of each dataset
for ll = 1:numel(uom.uo)                                                   % Iterate over all the data objects
    allMotTots(ll,1) = sum(abs(uom.uo{ll}.errs.motion.depth));             % Store the total depth deviations
    allMotTots(ll,2) = sum(abs(uom.uo{ll}.errs.motion.ap));                % Store the total ap deviations
    allMotTots(ll,3) = allMotTots(ll,1)+allMotTots(ll,2);                  % Store the total deviations
end
if isscalar(refDataset)
    if (1<=refDataset)&&(refDataset<=numel(uom.uo))
        refDataset = refDataset; % Keep it
    else
        warning('Bad reference dataset ID given!')
        refDataset = find(allMotTots(:,3) == min(allMotTots(:,3)), 1);     % Get the index for the dataset with the minimal total motion deviation
    end
else
    warning('Bad reference dataset ID given!')
    refDataset = find(allMotTots(:,3) == min(allMotTots(:,3)), 1);         % Get the index for the dataset with the minimal total motion deviation
end
fprintf('selected dataset %d as the global reference...\n', refDataset)
uom.uo{refDataset}.ensureMedian;                                           % Make sure the median image exists


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Store the results in a global vector that can be referenced

if p.globalOrient
    refMot = alignRefs(uom, refDataset, p);
    for ll = 1:numel(uom.uo)                                               % Iterate over all the data objects
        uom.uo{ll}.errs.motion.depth = uom.uo{ll}.errs.motion.depth ...
                                                           + refMot(ll,1); % Update the total depth deviations
        uom.uo{ll}.errs.motion.ap    = uom.uo{ll}.errs.motion.ap ...
                                                           + refMot(ll,2); % Update the total ap deviations
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Recalculate motion corrected videos

for ll = 1:numel(uom.uo)                                                   % Iterate over all the data objects
    uom.uo{ll}.correctResidualMotion('reCalc', true);  
end

if plotOpt 
    plotMotionOutput(uom, p.batchSz);                                      % Optional plotting code
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function refMot = alignRefs(uom, refDataset, p)

p.searchBlock   = min(4*p.searchBlock, 25);                                % Need a larger search radius for inter-session alignment
p.motionCorrect = true;
refImage        = [];                                                      % Initialize an array of reference images
padDiff         = [];
TMPmask         = 0;

for ll = 1:numel(uom.uo)                                                   % Need to get each dataset's reference image
    movieOut = makeGetMovieFunction(uom.uo{ll},p);                         % Temporarily load the data
    refTmp   = median(movieOut,3);                                         % Calculate the median of the data
    [refTmp, padTMP, tmpMask] = resizeFrame(uom, ll, refTmp);
    padDiff  = cat(2, padDiff.', padTMP(:)).';
    refImage = cat(3,refImage,refTmp./max(refTmp(:))); 
    TMPmask  = TMPmask|tmpMask;                                            
end
refImage(isnan(refImage)) = 0;
uoTMP = uo(refImage,'preCalcStats',false, 'drawMask', true);               % Make a temporary uo object that is a movie of reference images
uoTMP.mask = TMPmask;
uoTMP.findMotionCorrectionError('refImage', TMPmask.*refImage(:,:,refDataset), ...
                      'useMask', true, 'searchBlock', p.searchBlock, ...
                                        'mcType', p.mcType, 'batchSz', 1); % Get motion errors for ref images across datasets, passing through parameters

refMot = cat(2,uoTMP.errs.motion.depth(:),uoTMP.errs.motion.ap(:));        % Store the ap and depth deviations

refMot(:,1) = refMot(:,1); % + padDiff(:,1);
refMot(:,2) = refMot(:,2); % + padDiff(:,3);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

function [refIm, padDiff, uoMask] = resizeFrame(uom, ll, refIm)

padSize = zeros(numel(uom.uo), 4);
for kk = 1:numel(uom.uo)
    padSize(kk,:) = uom.uo{kk}.post.motionPad;
end
padGlobal = max(padSize,[],1);

padDiff = padGlobal - padSize(ll,:);
% padMask = padGlobal - padSize(1,:);
uoMask  = padarray(padarray(uom.uo{ll}.mask, padGlobal([1,3]),0,'pre'), ...
                                                  padGlobal([2,4]), 0, 'post');

refIm = padarray(padarray(refIm, ...
             padDiff([1,3]),0,'pre'), padDiff([2,4]), 0, 'post');
% refIm = refIm; %.*uoMask;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

function plotMotionOutput(uom, batchSz)

errsDepth = [];
errsAp    = [];
for ll = 1:numel(uom.uo)
    errsDepth = cat(1, errsDepth, uom.uo{ll}.errs.motion.depth(:));
    errsAp    = cat(1, errsAp,    uom.uo{ll}.errs.motion.ap(:));
end
errsDepth = reshape(ones(batchSz,1)*(errsDepth.'),[],1);
errsAp    = reshape(ones(batchSz,1)*(errsAp.'),[],1);

figure(1)
plot(uom.uo{1}.dt*(0:(size([errsDepth, errsAp],1)-1)), 100*[errsDepth./uom.uo{1}.frameSize(1), errsAp./uom.uo{1}.frameSize(2)])
figure(1), hold on
box off
ylabel('Offset (% field of view)'); xlabel('Frame number'); 
legend('Depth', 'AP')
figure(1), hold off

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
