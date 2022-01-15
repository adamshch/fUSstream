function varargout = getBasicStats(fuObj, varargin)


p = inputParser;                                                           % Set up an object to parse all of the various inputs
p.addParameter('baseSelect'  , 'median');                                  % Select what criteria to calculate the baseline image as
p.addParameter('basePrctile' , 25);                                        % For the situation that 'prctile' is chosen, which percentile should the baseline image be calculated using
parse(p,varargin{:});
p = p.Results;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Overall statistics

fuObj.stats.max    = max(fuObj.movie(:));                                  % Calculate the maximum value in the data
fuObj.stats.min    = min(fuObj.movie(:));                                  % Calculate the minimum value in the data
fuObj.stats.var    = var(fuObj.movie(:));                                  % Calculate the data valiance 
fuObj.stats.max    = max(fuObj.movie(:));                                  % Calculate the data mean
fuObj.stats.rstd   = robustSTD(fuObj.movie(:));                            % Calculate the robust standard deviation of the data
fuObj.stats.median = median(fuObj.movie(:));                               % Calculate the data mean

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Image-wise statistics

fuObj.stats.medianImg = median(fuObj.movie,3);                             % Calculate the median image
fuObj.stats.modeImg   = halfSampleMode(reshape(fuObj.movie,[],...
                                                       fuObj.movieLen).'); % Calculate the half-sample mode for each pixel in the image
fuObj.stats.modeImg   = reshape(fuObj.stats.modeImg(:),fuObj.frameSize);   % Reshape the vector of modes into an image
fuObj.stats.meanImg   = mean(fuObj.movie,3);                               % Calculate the mean image
fuObj.stats.varImg    = var(fuObj.movie,[],3);                             % Calculate the variance image
fuObj.stats.rstdImg   = robustSTD(fu.movie,3);                             % Calculate the per-pixel robust standard deviation
fuObj.stats.maxImg    = max(fuObj.movie,[],3);                             % Calculate the max image
fuObj.stats.minImg    = min(fuObj.movie,[],3);                             % Calculate the min image
fuObj.stats.prctImg   = prctile(fuObj.movie,p.basePrctile,3);              % Calculate the percentile of each pixel's histogram

if nargout > 0
    varargout{1} = fuObj.stats;
end

end