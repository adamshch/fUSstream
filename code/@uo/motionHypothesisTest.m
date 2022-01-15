function motTest = motionHypothesisTest(uo, varargin)

% motTest = motionHypothesisTest(uo, varargin)
%
% This function tests the motion via a test that the pixel histograms are
% similar between the high-motion and low-motion timepoints. In short: 
% if the overall behavior is constant and the same pixel is being imaged
% then we would expect that long batches of time should have similar 
% statistics. Thus the function breaks up the time-courses for each 
% pixel into "minHistNum"-sized batches, and compares using 
% "histogramShuffle". 
%
% Optional parameters are:
%   'motionCorrect' - Select if the motion corrected data should be used:
%                       true/false (default false)                               
%   'denoised'      - Choose to display denoised movie true/false 
%                       (default false) 
%   'plotOpt'       - Optional plotting flag. true/false (default false) 
%   'figNumber'     - Optional figure number to plot to (default 203) 
%   'motionCompare' - Choose which motion metric to correlate to (default
%                       'total')
%   'useMask'       - Choose to use the mask and only consider "inside"
%                       pixels. true/false (default true)
%   'minHistNum'    - Choose to use the mask and only consider "inside" 
%                       pixels (default 1000)
%   'distSelect'    - Choose the distance metric (default 'emd')
%   'shuffle'       - Choose whether to SHUFFLE the data to test the 
%                       metric. true/false (default false)   
%
% Output is a struct with
%   motTest.dists     - The histogram distances per pixel
%   motTest.distRatio - The ratio of high and low motion histogram 
%                       distances per pixel
%
% 2021 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

p = inputParser;                                                           % Set up an object to parse all of the various inputs
p.addParameter('motionCorrect' , false   );                                % Select if the motion corrected data should be used
p.addParameter('denoised'      , false   );                                % Choose to display denoised movie
p.addParameter('plotOpt'       , false   );                                % Optional plotting flag
p.addParameter('figNumber'     , 203     );                                % Optional figure number to plot to
p.addParameter('motionCompare' , 'total' );                                % Choose which motion metric to correlate to
p.addParameter('useMask'       , true    );                                % Choose to use the mask and only consider "inside" pixels
p.addParameter('minHistNum'    , 1000    );                                % Choose to use the mask and only consider "inside" pixels
p.addParameter('distSelect'    , 'emd'   );                                % Choose the distance metric
p.addParameter('shuffle'       , false   );                                % Choose whether to SHUFFLE the data to test the metric. 
parse(p,varargin{:});
p = p.Results;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Some initializations

if p.useMask;    uo.ensureMask();   end                                    % Make sure that a mask is available
movieOut = uo.makeGetMovieFunction(p);                                     % Get the movie given the requested options (denoised/motion corrected)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute motion vector and correlations

motVec    = getMotionVector(uo, p);                                        % Get the motion vector
motTest   = struct;                                                        % Initialize the motion test struct to store variables in
cutoffVal = prctile(abs(motVec), 50);                                      % The median separates the HIGH motion frames from the LOW motion frames (50/50 split for equality)

[dataParts1, dataParts2] = getIdxParts(motVec, cutoffVal, p.minHistNum,...
                                                                p.shuffle);% This function gets the intex parts
padX = p.motionCorrect*uo.post.motionPad(1);                               % Create a padding to compensate for motion correction movement (in x)
padY = p.motionCorrect*uo.post.motionPad(3);                               % Create a padding to compensate for motion correction movement (in y)
distSelect = p.distSelect;                                                 % Save distance selection option separately for use in a parfor loop
for ll = 1:size(movieOut,1)-p.motionCorrect*sum(uo.post.motionPad(1:2))    % LOOP over pixels
    TMPmov = movieOut((ll+padX), :, :);                                    % Pull out one time-trace
    parfor kk = padY+(1:size(movieOut,2)-p.motionCorrect*sum(uo.post.motionPad(3:4))) % LOOP over pixels (again)
        TMPdists(ll,kk) = histogramShuffle(TMPmov(1, kk,:), ...     
                dataParts1,dataParts2,'distSelect',distSelect);            % Get the distances for the {ll, kk}^th pixel
    end
    fprintf('.');
end
fprintf('done\n');
motTest.dists = TMPdists(:, padY+(1:size(movieOut,2)-p.motionCorrect*sum(uo.post.motionPad(3:4))));% 
motTest.distRatio = computeRatioDists(motTest);                            % Compute the image of distance ratios


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Optional plotting

if p.plotOpt
    if isempty(p.figNumber); p.figNumber = 203; end                        % Make sure the figure number is not empty
    %% First a histogram
    figure(p.figNumber)
    subplot(1,2,1), histogram(vec(motTest.distRatio))
    box off; set(gca,'XLim',1.1*[-0.5,max(vec(motTest.distRatio))])        % Remove box and set limits
    title('\rho motion and moving averages')                               % Set the title
   
    subplot(1,2,2), imagesc(motTest.distRatio)
    pbaspect([uo.scanArea(2),uo.scanArea(1),1]);                           % Set the aspect ratio
    axis off;                                                              % Remove axes 
    title(sprintf('Score map (mean)'))                                     % Give a title to the image
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function motVec = getMotionVector(uo, p)

[motDep, motLat, ~] = uo.getMotionVectors();                               % Check for motion errors if needed

motDep  = motDep - motDep(1);                                              % Remove average motion for this session (in-session metric)
motLat  = motLat - motLat(1);                                              % Remove average motion for this session (in-session metric)

switch lower(p.motionCompare)                                              % Choose which motion value to return for comparison
    case 'depth';           motVec = motDep;                               % Only the depth motion is counted
    case {'lateral','lat'}; motVec = motLat;                               % Only the laeral motion is counted
    case 'total';           motVec = sqrt(motDep.^2 + motLat.^2);          % Total distance of the vector (l2 norm)
    otherwise;              error('Bad comparison motion metric provided.')% Return an error
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dataParts1, dataParts2] = partitionIdxs(idx1, idx2, minHistNum)

M1 = numel(idx1);                                                          % Get number of datapoints in first dataset
M2 = numel(idx2);                                                          % Get number of datapoints in second dataset

[blkSizes1, dataBlocks1]  = getBlockSizes(idx1);                           % Get block sizes for first dataset
[blkSizes2, dataBlocks2]  = getBlockSizes(idx2);                           % Get block sizes for second dataset
minBlock = max(min([blkSizes1(:); blkSizes2(:)]), minHistNum);             % Compute the minimal block size

numTotalParts1 = ceil(M1/minHistNum);                                      % Compute the total number of small blocks that need to be combined to make each partition
numTotalParts2 = ceil(M2/minHistNum);                                      % Compute the total number of small blocks that need to be combined to make each partition
dataParts1     = partitionData(dataBlocks1, numTotalParts1, minHistNum,...
                                                                 minBlock);% Partition the data from the first dataset
dataParts2     = partitionData(dataBlocks2, numTotalParts2, minHistNum,...
                                                                 minBlock);% Partition the data from the second dataset
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dataParts1, dataParts2] = getIdxParts(motVec, cutoffVal, minHistNum, varargin)

if nargin > 3;    shuffleOpt = varargin{1};                                % Check if a shuffle test is being requested
else;             shuffleOpt = false;                                      % Default is NO shuffle test
end

idx1 = find(abs(motVec)< cutoffVal);                                       % Get indices of low-motion frames
idx2 = find(abs(motVec)>=cutoffVal);                                       % Get indices of high-motion frames

idx1 = vec(idx1); idx2 = vec(idx2);                                        % Ensure idx1 and idx 2 are COLUMN VECTOS <-- This will mess up indexing down the line!
if max(numel(idx1), numel(idx2)) < minHistNum
    minHistNum = floor(max(numel(idx1), numel(idx2))/2);
end

[dataParts1, dataParts2] = partitionIdxs(idx1, idx2, minHistNum);          % Partition the indices for histogram analysis
dataParts1 = removeEmptyParts(dataParts1);                                 % Remove any empty partitions from the first dataset (low-motion)
dataParts2 = removeEmptyParts(dataParts2);                                 % Remove any empty partitions from the second dataset (high-motion)
                                                          
if shuffleOpt                                                              % If the shuffle option is on, shuffle the blocks between the two conditions
    [dataParts1, dataParts2] = shuffleDataParts(dataParts1, dataParts2);   % Shuffle the data parts
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [blkSizes, dataBlocks] = getBlockSizes(idxIn)

idxJump    = [diff(idxIn(:)'),10] > 1;
blkSizes   = find(idxJump==1);
dataBlocks = cell(numel(blkSizes),1);

blkEnds = blkSizes;                                                        % Initialize the block size
if blkEnds(1)>0;              blkEnds = [0,blkSizes];             end      % Make sure the first block starts at the first index
if blkEnds(end)<numel(idxIn); blkEnds = [blkSizes, numel(idxIn)]; end      % Make sure the final block starts at the final index
    
for ll = 1:size(dataBlocks,1)
    dataBlocks{ll} = idxIn((blkEnds(ll)+1):blkEnds(ll+1));
    blkSizes(ll)   = numel(dataBlocks{ll});
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function to partition the data

function dataParts = partitionData(dataBlocks, numTotalParts, minHistNum, minBlock)%, varargin)

dataParts = cell(numTotalParts, 1);                                        % Initialize the partitions to the elements of an empty cell

% This next for loop looks for blocks of data that are bigger than
% 2*minBlock and breaks them up into blcoks of size minBlock (with the
% residual being placed in the last block). s
for ll = 1:size(dataBlocks,1)                                              % Loop over data blocks
    blkSz = numel(dataBlocks{ll});                                         % Get the size of data blocks
    if blkSz > 2*minBlock                                                  % If there are enough elements in this block to break up into two blocks...
        subBlkNum = floor(blkSz/minBlock);                                 % ... Figure out how many blocks can be made
        TMP       = cell(subBlkNum,1);                                     % ... Temporarily store the number of blocks to break into
        for kk = 1:(size(TMP,1)-1)                                         % ... For each sub-block to be made
            TMP{kk} = dataBlocks{ll}(((kk-1)*minBlock+1):(kk*minBlock));   % ... ... Separate that block away from the rest and save
        end
        TMP{end} = dataBlocks{ll}((kk*minBlock):end);                      % ... The residual should not be *smaller*
        
        dataBlocks = cat(1, dataBlocks(1:(ll-1),:), TMP, ...
                                                dataBlocks((ll+1):end,:)); % ... Replace the block with the diced blocks
    else                                                                   % Otherwise: Do nothing
    end
end

% This next loop partitions out the blocks into sets that have at least
% minHistNum elements.
idxUnused = 1:size(dataBlocks,1);                                          % Initialize a set of contiguous blocks to partition out
for ll = 1:size(dataParts,1)                                               % Loop over the partitions:
    dataParts{ll,1} = [];                                                  % ... initialize the partition indicies (location)
    while (numel(dataParts{ll}) < minHistNum)&&(~isempty(idxUnused))       % ... WHILE the partition is too small AND there are still blocks to partition
        idxBlk        = randsample(idxUnused, 1);                          % ... ... Pick a block at random
        idxUnused     = setdiff(idxUnused, idxBlk);                        % ... ... Add the block to the current partition
        dataParts{ll} = cat(1,vec(dataParts{ll}),dataBlocks{idxBlk,1});    % ... ... Add index set to that partition
    end
end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dataParts] = removeEmptyParts(dataParts)
% function to remove any empty partitions

dataParts = dataParts(~cellfun(@isempty, dataParts));                   % Remove empty partitions from data 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function distRatio = computeRatioDists(motTest) 

distRatio = zeros(size(motTest.dists,1), size(motTest.dists,2));           % Initialize the distance ration matrix

for ll = 1:size(motTest.dists,1)
    for kk = 1:size(motTest.dists,2)
        distRatio(ll,kk) = nanmean([vec(motTest.dists(ll,kk).in + ...
              diag(diag(nan*eye(size(motTest.dists(ll,kk).in))))); ...
                          vec(motTest.dists(ll,kk).out + ...
                    diag(diag(nan*eye(size(motTest.dists(ll,kk).out)))))]);% 
        distRatio(ll,kk) = mean(motTest.dists(ll,kk).between(:))./ ...
                                                          distRatio(ll,kk);% 
    end
    fprintf('.');
end
fprintf('done\n');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

function [dataParts1, dataParts2] = shuffleDataParts(dataParts1, dataParts2)

dataPartsAll = cat(1,dataParts1, dataParts2);
idSamp       = randsample(numel(dataPartsAll), numel(dataParts1));
dataParts1   = dataPartsAll(idSamp);
dataParts2   = dataPartsAll(setdiff(1:numel(dataPartsAll),idSamp));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CODE GRAVEYARD



%         TMPdists(ll,kk) = histogramShuffle(movieOut((ll+padX),...
%                          (kk+padY),:), dataParts1,dataParts2,'distSelect',distSelect);% Get the distances for the {ll, kk}^th pixel
