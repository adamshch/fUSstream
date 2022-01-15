function varargout = plotBurstErrors(uo, varargin)

% errIDs = plotBurstErrors(fuobj, varargin)
%
% 
% 2019 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

p = inputParser;                                                           % Set up an object to parse all of the various inputs
p.addParameter('figNo', 101);                                              % Can choose the figure number to plot to
p.addParameter('nPlot', 6);                                                % Choose the number of examples to plot
p.addParameter('nRows', 3);                                                % Choose the number of subfigure rows in the plot
p.addParameter('cLim', 'good');                                            % Can opt to scale the movies in different ways ('all', 'good', 'bad')
parse(p,varargin{:});
p = p.Results;

nCols = ceil(p.nPlot/p.nRows);                                             % Calculate the number of columns

allIdx  = 1:uo.movieLen;                                                   % Make a vector of all index locations
badIdx  = allIdx(uo.errs.burst);                                           % Select all "good" indecies 
goodIdx = allIdx(~uo.errs.burst);                                          % Select all "good" indecies 

goodPlotIdx = randsample(goodIdx,min(p.nPlot, numel(goodIdx)));            % Subselect from the good frames
nPgood      = numel(goodPlotIdx);                                          % Update number of good frames to show
badPlotIdx  = randsample(badIdx,min(p.nPlot, numel(badIdx)));              % Subselect from the bad frames
nPbad       = numel(badPlotIdx);                                           % Update number of bad frames to show


if (uo.meta.loadMov)&&(~isempty(uo.movie))
    getDataBlock = @(z) uo.movie(:,:,z);                                   % Extract a movie block
else
    tmpMov       = uo.meta.movMatFile.(uo.meta.movVarName);
    getDataBlock = @(z) tmpMov(:,:,z);                                     % Extract a movie block
end

clims = chooseColorLimits(getDataBlock, goodPlotIdx, badPlotIdx, p.cLim);  % Compute the color limits

h = figure(p.figNo);
for kk = 1:nCols
    for ll = 1:min(p.nRows,nPgood-((kk-1)*p.nRows))
        idxToPlot = ll+(kk-1)*p.nRows;
        subToPlot = (2*(ll-1)*nCols)+kk;
        subplot(p.nRows,2*nCols,subToPlot), ...
                  imagesc(getDataBlock(goodPlotIdx(idxToPlot)), clims)
        title(sprintf('Good Example %d', idxToPlot))
        axis image; axis off; colormap gray;
    end
    for ll = 1:min(p.nRows,nPbad-((kk-1)*p.nRows))
        idxToPlot = ll+(kk-1)*p.nRows;
        subToPlot = (2*(ll-1)*nCols)+kk+nCols;
        subplot(p.nRows,2*nCols,subToPlot), ...
                   imagesc(getDataBlock(badPlotIdx(idxToPlot)), clims)
        title(sprintf('Bad Example %d', idxToPlot))
        axis image; axis off; colormap gray;
    end
end

if nargout > 0;    varargout{1} = h;    end                                % Output the figure handle if requested

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function clims = chooseColorLimits(getDataBlock, goodPlotIdx, badPlotIdx, cOpt);

if strcmp(cOpt,'good')
    clims = [min(reshape(getDataBlock(goodPlotIdx),1,[])), ...
                            max(reshape(getDataBlock(goodPlotIdx),1,[]))]; % 
elseif strcmp(cOpt,'bad')
    clims = [min(reshape(getDataBlock(badPlotIdx),1,[]))...
                             max(reshape(getDataBlock(badPlotIdx),1,[]))]; % 

else
    cvals_max = max(max(reshape(getDataBlock(goodPlotIdx),1,[])), ...
                             max(reshape(getDataBlock(badPlotIdx),1,[]))); % 
    cvals_min = min(min(reshape(getDataBlock(goodPlotIdx),1,[])), ...
                             min(reshape(getDataBlock(badPlotIdx),1,[]))); % 
    clims = [cvals_min, cvals_max];
end

end