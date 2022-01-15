function trialTraces = alignTrials(fuObj, varargin)

% function fuObj = alignTrials(fuObj, varargin)
% 
% Function to align trials to a behavioral variable.
% 
% 2020 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

p = inputParser;                                                           % Set up an object to parse all of the various inputs
p.addParameter('plotLevel'   , 0);                                         % How often should plots show up? (default 0)
p.addParameter('plotOpt'     , false);                                     % Should plots be displayed at all (default false)
p.addParameter('plotGroup'   , 'byTrial');                                 % Group plotting
p.addParameter('numPlot'     , 1);                                         % Number of plots to display
p.addParameter('toPlot'      , 'PCA');                                     % What to plot (default principal components)
p.addParameter('idxToPlot'   , 1);                                         % Which index of the matrix to plot (default 1)
parse(p,varargin{:});
p = p.Results;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get the time-series of trial starts

EventShape = fuObj.event2timeseries()==0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% extract useful information from the event shape data

evtStarts = find((EventShape(1:(end-1))==0)&(EventShape(2:end)==1)) + 1;   % Get the indecies of all the start frames for each trial
evtStarts = evtStarts(1:(end-1));
evtEnds   = find((EventShape(1:(end-1))==1)&(EventShape(2:end)==0)) + 25;  % Get the indecies of all the end frames for each trial
evtEnds = evtEnds(2:end);
trialLens = evtEnds-evtStarts;                                             % Get the lengths of each trial
[~,ixLen] = sort(trialLens,'descend');                                     % Sort the trial length in descending order

trialTraces = genTrialTraces(fuObj,p,evtStarts,evtEnds);


if p.plotOpt
    switch lower(p.plotGroup)
        case 'bytrial'
            nPlot = min(numel(trialLens),p.numPlot);                                   % Figure out number to plot
            for ll = 1:nPlot
                subplot(nPlot,1,ll), plot(fuObj.dt*(0:(size(trialTraces{ll},2)-1)), trialTraces{ll}')
                box off
            end
            title(sprintf('Trial %d',ll))
            xlabel('Time (s)')
        case 'bypc'
            nPlot = min(numel(p.idxToPlot),size(trialTraces{1},1));        % Figure out number to plot
            nLine = min(numel(trialLens),p.numPlot);                       % Figure out number of lines per plot
            for kk = 1:nPlot
                subplot(nPlot,1,kk), cla
                subplot(nPlot,1,kk), hold on
                for ll = 1:nLine
                    subplot(nPlot,1,kk), plot(fuObj.dt*(0:(size(trialTraces{ll},2)-1)), trialTraces{ll}(kk,:)')
                    box off
                end
                title(sprintf('%s: Component %d',p.toPlot,kk))
                xlabel('Time (s)')
                subplot(nPlot,1,kk), hold off
            end
        otherwise
    end
    linkaxes
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function trialTraces = genTrialTraces(fuObj,p,evtStarts,evtEnds);

nTrials = numel(evtStarts); % Get number of trials
trialTraces = cell(nTrials,1);

 switch lower(p.toPlot)
    case 'pca'
        traceToPlot = fuObj.PCAout.temporal(p.idxToPlot,:);
    case 'area'
        movTMP      = reshape(fuObj.movie,[],fuObj.movieLen);
        traceToPlot = mean(movTMP(p.idxToPlot,:),1);
    case 'pixels'
        movTMP      = reshape(fuObj.movie,[],fuObj.movieLen);
        traceToPlot = movTMP(p.idxToPlot,:);
    otherwise
        p.toPlot = 'pca';
        trialTraces = genTrialTraces(fuObj,p,evtStarts,evtEnds);
end

for ll = 1:nTrials % Loop over trials to extract trace
   trialTraces{ll} = traceToPlot(:,evtStarts(ll):evtEnds(ll));
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
