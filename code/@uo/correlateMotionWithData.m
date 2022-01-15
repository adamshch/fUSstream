function motCorr = correlateMotionWithData(uo, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

p = inputParser;                                                           % Set up an object to parse all of the various inputs
p.addParameter('movingWindow'  , 50     );                                 % Select the size of the moving averages
p.addParameter('motionCorrect' , false   );                                % Select if the motion corrected data should be used
p.addParameter('denoised'      , false   );                                % Choose to display denoised movie
p.addParameter('plotOpt'       , true    );                                % Optional plotting flag
p.addParameter('figNumber'     , 203     );                                % Optional figure number to plot to
p.addParameter('motionCompare' , 'total' );                                % Choose which motion metric to correlate to
p.addParameter('useMask'       , true    );                                % Choose to use the mask and only consider "inside" pixels
parse(p,varargin{:});
p = p.Results;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Some initializations

if p.useMask;    uo.ensureMask();   end                                    % Make sure that a mask is available
movieOut = uo.makeGetMovieFunction(p);                                     % Get the movie given the requested options (denoised/motion corrected)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute motion vector and correlations

motVec  = getMotionVector(uo, p);                                           % Get the motion vector
motCorr = computeMotionCorrsSubFunc(uo, movieOut, motVec, p);              % Compute motion correlations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Optional plotting

if p.plotOpt
    if isempty(p.figNumber); p.figNumber = 203; end                        % Make sure the figure number is not empty
    %% First a histogram
    figure(p.figNumber)
    subplot(2,2,1), histogram(vec(motCorr.mean))
    box off; set(gca,'XLim',[-1,1])                                        % Remove box and set limits
    title('\rho motion and moving averages')                               % Set the title
    subplot(2,2,2), histogram(vec(motCorr.var))
    box off; set(gca,'XLim',[-1,1])                                        % Remove box and set limits
    title('\rho motion and moving variances')                              % Set the title
   
    subplot(2,2,3), imagesc(motCorr.mean,[-1,1])
    pbaspect([uo.scanArea(2),uo.scanArea(1),1]);                           % Set the aspect ratio
    axis off;                                                              % Remove axes 
    title(sprintf('Correlation map (mean)'))                               % Give a title to the image
    subplot(2,2,4), imagesc(motCorr.var,[-1,1])
    pbaspect([uo.scanArea(2),uo.scanArea(1),1]);                           % Set the aspect ratio
    axis off;                                                              % Remove axes 
    title(sprintf('Correlation map (var)'))                                % Give a title to the image
end


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function motVec = getMotionVector(uo, p)

[motDep, motLat, ~] = uo.getMotionVectors();                               % Check for motion errors if needed

switch lower(p.motionCompare)
    case 'depth'
        motVec = motDep;
    case {'lateral','lat'}
        motVec = motLat;
    case 'total'
        motVec = sqrt(motDep.^2 + motLat.^2);
    otherwise
        error('Bad comparison motion metric provided.')
end

motVec = movmean(motVec, p.movingWindow);
motVec = motVec(1:p.movingWindow:end);

end


function motCorr = computeMotionCorrsSubFunc(uo, movieOut, motVec, p)

motCorr.mean = nan*zeros(uo.frameSize);
motCorr.var  = nan*zeros(uo.frameSize);
if p.motionCorrect
    lStart = uo.post.motionPad(1);
    kStart = uo.post.motionPad(2);
else
    lStart = 0;
    kStart = 0;
end

%% Run through all pixels

fprintf('Computing correlations.\n')
for ll = 1:uo.frameSize(1) 
    for kk = 1:uo.frameSize(2)  
        if p.useMask&&(uo.mask(ll,kk)==0)
            % Do nothing
        else
            traceMean = vec(movmean(movieOut(lStart+ll,kStart+kk,:)/max(vec(movieOut(lStart+ll,kStart+kk,:))), p.movingWindow));
            traceMean = traceMean(1:p.movingWindow:end);
            traceVar  = movvar(movieOut(lStart+ll,kStart+kk,:)/max(vec(movieOut(lStart+ll,kStart+kk,:))),  p.movingWindow);
            traceVar  = traceVar(1:p.movingWindow:end);
            tmpMean   = corrcoef(motVec, traceMean);
            tmpVar    = corrcoef(motVec, traceVar);
            motCorr.mean(ll,kk) = tmpMean(1,2);
            motCorr.var(ll,kk)  = tmpVar(1,2);
        end
    end
    fprintf('.')
end
fprintf('done.\n')

end
