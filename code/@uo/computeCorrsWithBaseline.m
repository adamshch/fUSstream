function varargout = computeCorrsWithBaseline(fuObj, varargin)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

p = inputParser;                                                           % Set up an object to parse all of the various inputs
p.addParameter('baseSelect'  , 'median');                                  % Select what criteria to calculate the baseline image as
p.addParameter('basePrctile' , 25);                                        % For the situation that 'prctile' is chosen, which percentile should the baseline image be calculated using
parse(p,varargin{:});
p = p.Results;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch p.baseSelect
    case 'mode'                                                            % If the mode is selected...
        baseImg = halfSampleMode(reshape(fuObj.movie,[],fuObj.movieLen).');%  ...Calculate the half-sample mode for each pixel in the image
        baseImg = reshape(baseImg(:),fuObj.frameSize);                     %  ...Reshape the vector of modes into an image
    case 'mean'                                                            % If the mean is selected...
        baseImg = mean(fuObj.movie,3);                                     %  ...Calculate the mean of each pixel in the image
    case 'prctile'                                                         % If a percentile is selected
        baseImg = prctile(fuObj.movie,p.basePrctile,3);                    %  ...Calculate the percentile of each pixel's histogram1
    otherwise                                                              % The default is the median, so if nothing else is selected...
        baseImg = median(fuObj.movie,3);                                   %  ...Calculate the median of each pixel
end

baseMean = mean(baseImg(:));                                               % Get the mean of the base image
baseSTD  = std(baseImg(:));                                                % Get the standard deviation of the base image
imgMeans = sum(sum(fuObj.movie,1),2)/prod(fuObj.frameSize);                % Get the mean for each image
imgVars  = sum(sum(bsxfun(@plus,fuObj.movie,-imgMeans).^2,1),2);           % Get the variance for each image

fuObj.stats.baseCorr = sum(sum(bsxfun(@times, ...
                              bsxfun(@plus,fuObj.movie, -imgMeans), ...
                                              (baseImg - baseMean)),1),2); % Calculate the inner products between the two norms
fuObj.stats.baseCorr = fuObj.stats.baseCorr./(sqrt(imgVars)*baseSTD);      % Calculate the pearson correlation
fuObj.stats.baseCorr = fuObj.stats.baseCorr(:);                            % Reshape into a vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output parsing

if nargout > 0
    varargout{1} = fuObj.stats.baseCorr;
end

end