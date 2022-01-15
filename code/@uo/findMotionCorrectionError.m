function uo = findMotionCorrectionError(uo, varargin)

% function fuObj = findMotionCorrectionError(fuObj, varargin)
% 
% Function that seeks out and quantifies residual sub-pixel motion errors
% in functional ultrasound imaging. 
% 
% 2020 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

p = inputParser;                                                           % Set up an object to parse all of the various inputs
p.addParameter('mcType'      , 'batchMedian');                             % Option of what motion correction type to use ('median', 'batchMedian'). This list will expand
p.addParameter('frameSel'    , []);                                        % Option to select which frames to motion correct
p.addParameter('refImage'    , []);                                        % Option to input a user-defined reference image (useful for multi-dataset computations)
p.addParameter('globalOrient', false);                                      % Option to globally orient the data
p.addParameter('medBaseLine' , 'firstmedian');                             % Select the reference image
p.addParameter('batchSz'     , 5);                                         % Select batch sizes to reduce computation time
p.addParameter('useDenoised' , false);                                     % Option to motion correct based on the denoised data
p.addParameter('searchBlock' , 5);                                         % Set the limit (in pixels) that the motion correction searches for shifts over
p.addParameter('useMask'     , true);                                      % Option to use the mask in order to compute shifts only based on brain motion (ignore out-of-brain motion)                           
parse(p,varargin{:});
p = p.Results;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Some initializations

if p.useDenoised&&(~uo.doesDenoiseExist())                                 % Check if denoised image exists/is requested
    warning('Requested denoised movie for PCA, but no denoised movie available. Running wavelet denoising with basic features...')
    uo.denoiseTracesWavelet();                                             % If requested and does NOT exists, compute it
end

movieOut = makeGetMovieFunction(uo,p); 
uo.ensureMedian;
if isempty(p.refImage); dataMedian = computeComparisonImage(uo,movieOut,p);% If no reference image provided, compute the reference image given the inputs
else;                   dataMedian = p.refImage;                           % Otherwise, use the reference image
end

dataMedian = dataMedian./max(dataMedian(:));                               % Normalize the median for better conditioning

if p.useMask;    uo.ensureMask();   end                                    % Make sure that a mask is available

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start analysis

switch lower(p.mcType)
    case 'fullmedian'
        frameOffsets = fullMedianCorrTest(uo, movieOut, p.frameSel,....
                                                dataMedian, p.searchBlock, p.useMask);% Run one frame at a time
    case 'batchmedian'  
        frameOffsets = batchMedianCorrTest(uo, movieOut, p.frameSel,...
                                     p.batchSz, dataMedian, p.searchBlock, p.useMask);% Run batch median correlation test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% TO IMPLEMENT IN THE FUTURE: %%%%
%     case 'runningcorr'                  
%     case 'normcorre'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    otherwise
        error('Unknown motion error detection scheme')
end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

if nargout > 0
%     varargout{1} = frameOffsets; %fuObj.errs.motion;
end
uo.errs.motion.depth  = frameOffsets(1,:);
uo.errs.motion.ap     = frameOffsets(2,:);
uo.errs.motion.params = p;


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extra functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 

function frameOffsets = fullMedianCorrTest(uo, movieOut, frameSel, dataMedian, searchBlock, useMask)

if isempty(dataMedian); dataMedian     = uo.stats.medImg;     end          % If needed, calculate the median of the data
if isempty(frameSel);   frameSel       = 1:uo.movieLen;       end          % If needed, create a vector of frames to operate on
if isempty(searchBlock); searchBlock   = 10;                  end          % If needed, set a search block to compute correlations over (max motion)

frameOffsets = nan(2,uo.movieLen);                                         % Offsets saved as 'x,y' pairs


if useMask
    dataMedian = uo.mask.*dataMedian;                                      % Update the median image to be masked
    for ll = 1:numel(frameSel)
        tmpFrame = uo.mask.*movieOut(:,:,frameSel(ll));
        frameOffsets(:,frameSel(ll)) = motionFitSingleTest(dataMedian,...
                          tmpFrame./max(tmpFrame(:)), 'cont',...
                                              'searchBlock', searchBlock); % Get the motion vector estimate for the ll^th frame using a mask
    end
else
    for ll = 1:numel(frameSel)
        tmpFrame = movieOut(:,:,frameSel(ll));
        frameOffsets(:,frameSel(ll)) = motionFitSingleTest(dataMedian,...
                    tmpFrame./max(tmpFrame(:)), 'cont', 'searchBlock',...
                                                             searchBlock); % Get the motion vector estimate for the ll^th frame
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Batch medians

function frameOffsets = batchMedianCorrTest(uo, movieOut, frameSel, batchSz, dataMedian, searchBlock, useMask)

if isempty(batchSz);     batchSz    = round(uo.movieLen/100); end          % If needed, calculate the median of the data
if isempty(dataMedian);  dataMedian = uo.stats.medImg;        end          % If needed, calculate the median of the data
if isempty(frameSel);    frameSel   = 1:uo.movieLen;          end          % If needed, create a vector of frames to operate on
if isempty(searchBlock); searchBlock   = 10;                  end          % If needed, set a search block to compute correlations over (max motion)

dataMedian = dataMedian - nanmean(dataMedian(:));

if useMask;  uo.ensureMask();                  end                         % Make sure that a mask is available
if useMask;  dataMedian = uo.mask.*dataMedian; end                         % Update the median image to be masked
    
offCalcType  = 'cont';
numBlks      = ceil(uo.movieLen/batchSz);                                  % Compute the number of blocks
frameOffsets = nan(2,numBlks);                                             % Offsets saved as 'x,y' pairs

fprintf('Computing motion deviations:\n')
for ll = 1:numBlks
    blkIdx = ((ll-1)*batchSz + 1):(ll*batchSz);                            % Get current frame block indexing
    blkIdx = intersect(frameSel, blkIdx);                                  % Find which selected frames are in that batch
    if ~isempty(blkIdx)                                                    % Sweep over batches
        currMed  = nanmedian(movieOut(:,:,blkIdx),3);
        currMed  = currMed - nanmean(currMed(:));
        if useMask;  currMed = uo.mask.*currMed; end                       % Update the current image to be masked
        currMed = currMed./max(currMed(:));                                % Normalize image for better conditioning
        frameOffsets(:,ll) = motionFitSingleTest(dataMedian, currMed, ... 
                                 offCalcType, 'searchBlock', searchBlock); % Estimate the motion correction for the current frame
    end
    if (ll>2)&&(mod(ll,100) == 0);   fprintf('.');  end
    if (ll>2)&&(mod(ll,15000) == 0); fprintf('\n'); end
end
fprintf('done.\n')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function to compute the comparison image to compute offsets W.R.T. 

function dataMedian = computeComparisonImage(uo, movieOut, p)

if p.useDenoised 
    switch lower(p.medBaseLine)
        case 'allmedian'
            dataMedian = median(movieOut,3);                               % Calculate the median of the data
        case 'firstmedian'
            dataMedian = median(movieOut(:,:,1:100), 3);                   % Calculate the median of the data
        case 'lastmedian'
            dataMedian = median(movieOut(:,:,end-99:end), 3);              % Calculate the median of the data
        otherwise
            warning('invalide baseline chosen: using the global median....\n')
            dataMedian = median(movieOut,3);                               % Calculate the median of the data
    end
else
    switch lower(p.medBaseLine)
        case 'allmedian'
            dataMedian = uo.stats.medImg;                                  % Calculate the median of the data
        case 'firstmedian'
            dataMedian = median(movieOut(:,:,1:100), 3);                   % Calculate the median of the data
        case 'lastmedian'
            dataMedian = median(movieOut(:,:,end-99:end), 3);              % Calculate the median of the data
        otherwise
            warning('invalide baseline chosen: using the global median....\n')
            dataMedian = uo.stats.medImg;                                  % Calculate the median of the data
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
