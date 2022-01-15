function distVal = histogramShuffle(dataStream, idx1, idx2, varargin)

% distVal = histogramShuffle(dataStream, idx1, idx2, ['distSelect', distSelect])
% 
% Compute a suffled distance between two datasets defined by indexing the
% full dataset 'dataStream' with 2 sets of inices 'idx1' and 'idx2'.  The
% distance is selected by 'distSelect' (default 'emd'). 
% 
% 2020 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

p = inputParser;                                                           % Set up an object to parse all of the various inputs
p.addParameter('distSelect'     , 'emd'      );                            % Select which destance metric to use
parse(p,varargin{:});
p = p.Results;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Partition data

if any(reshape(cellfun(@isempty,idx1),[],1))                               % Check for empty partitions in idx1
    fprintf('Empty partitions!')
end
if any(reshape(cellfun(@isempty,idx2),[],1))                               % Check for empty partitions in idx2
    fprintf('Empty partitions!')
end   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loop through and get distnace values

distVal.in      = computeDistValsLoop(dataStream, idx1, idx1, p);        % Interior distances (between blocks of low-motion frames)
distVal.out     = computeDistValsLoop(dataStream, idx2, idx2, p);        % Interior distances (between blocks of high-motion frames)
distVal.between = computeDistValsLoop(dataStream, idx1, idx2, p);        % Exterior distances (between blocks of high- and low-motion frames)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function to partition the data


function distVal = computeDistValsLoop(dataStream, dataParts1, dataParts2, p)

distVal = nan(size(dataParts1,1), size(dataParts2,1));
for ll = 1:size(dataParts1,1)
    for kk = 1:size(dataParts2,1)
        if isempty(dataParts1{ll})||isempty(dataParts2{kk})
            error('Some of the partitions are empty!')
        end
        distVal(ll,kk) = histogramDistance(dataStream(dataParts1{ll}), ...
                             dataStream(dataParts2{kk}), 'distSelect', ...
                                                            p.distSelect); % Compute the distance between all the bins
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% OLD CODE GRAVEYARD::::


% M1 = numel(data1);                                                         % Get number of datapoints in first dataset
% M2 = numel(data2);                                                         % Get number of datapoints in second dataset
% 
% [blkSizes1, dataBlocks1]  = getBlockSizes(data1, idx1);                    % Get block sizes for first dataset
% [blkSizes2, dataBlocks2]  = getBlockSizes(data2, idx2);                    % Get block sizes for second dataset
% minBlock = max(min([blkSizes1(:); blkSizes2(:)]), p.minHistNum);           % Compute the minimal block size
% 
% numTotalParts = ceil(min(M1,M2)/p.minHistNum);                             % Compute the total number of small blocks that need to be combined to make each partition
% dataParts1    = partitionData(dataBlocks1, numTotalParts, p.minHistNum,...
%                                                                 minBlock); % Partition the data from the first dataset
% dataParts2    = partitionData(dataBlocks2, numTotalParts, p.minHistNum,...
%                                                                 minBlock); % Partition the data from the second dataset

% function dataParts = partitionData(dataBlocks, numTotalParts, minHistNum, minBlock)
% 
% dataParts = cell(numTotalParts, 2);
% 
% for ll = 1:size(dataBlocks,1)
%     blkSz = numel(dataBlocks{ll,1});
%     if blkSz > 2*minBlock
%         subBlkNum = floor(blkSz/minBlock);
%         TMP = cell(subBlkNum,2);
%         for kk = 1:(size(TMP,1)-1)
%             TMP{kk,1} = ...
%                       dataBlocks{ll,1}(((kk-1)*minBlock+1):(kk*minBlock)); %
%             TMP{kk,2} = ...
%                       dataBlocks{ll,2}(((kk-1)*minBlock+1):(kk*minBlock)); %
%         end
%         TMP{end,1} = dataBlocks{ll,1}((kk*minBlock):end);                  % The residual should not be *smaller*
%         TMP{end,2} = dataBlocks{ll,2}((kk*minBlock):end);                  % The residual should not be *smaller*
%         
%         dataBlocks = cat(1, dataBlocks(1:(ll-1),:), TMP, ...
%                                                 dataBlocks((ll+1):end,:)); % Replace the block with the diced blocks
%         
%     else                                                                    % Do nothing
%     end
%     
% end
% 
% idxUnused = 1:size(dataBlocks,1);                                          % Initialize a set of contiguous blocks to partition out
% for ll = 1:size(dataParts,1)                                               % Loop over the partitions:
%     dataParts{ll,1} = [];                                                  % ... initialize the partition indicies (location)
%     dataParts{ll,2} = [];                                                  % ... initialize the partition values 
%     while (numel(dataParts{ll,1}) < minHistNum)&&(~isempty(idxUnused))     % ... WHILE the partition is too small AND there are still blocks to partition
%         idxBlk    = randsample(idxUnused, 1);                              % ... ... Pick a block at random
%         idxUnused = setdiff(idxUnused, idxBlk);                            % ... ... Add the block to the current partition
%         
%         dataParts{ll,1} = cat(1,dataParts{ll,1},dataBlocks{idxBlk,1});     % ... ... Add index set to that partition
%         dataParts{ll,2} = cat(1,dataParts{ll,2},dataBlocks{idxBlk,2});     % ... ... Add value set to that partition
%     end
% end
% 
% end


% function [blkSizes, dataBlocks] = getBlockSizes(dataIn, idxIn)
% 
% idxJump    = [diff(idxIn(:)'),10] > 1;
% blkSizes   = find(idxJump==1);
% dataBlocks = cell(numel(blkSizes),2);
% 
% blkEnds = blkSizes;                                                        % Initialize the block size
% if blkEnds(1)>1;              blkEnds = [1,blkSizes];             end      % Make sure the first block starts at the first index
% if blkEnds(end)<numel(idxIn); blkEnds = [blkSizes, numel(idxIn)]; end      % Make sure the final block starts at the final index
%     
% for ll = 1:size(dataBlocks,1)
%     dataBlocks{ll,1} = idxIn(blkEnds(ll):blkEnds(ll+1));
%     dataBlocks{ll,2} = dataIn(blkEnds(ll):blkEnds(ll+1));
% end
% 
% end
