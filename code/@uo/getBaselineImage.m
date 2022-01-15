function baseImg = getBaselineImage(fuObj, varargin)

p = inputParser;                                                           % Set up an object to parse all of the various inputs
p.addParameter('baseSelect'  , 'median');                                  % Select what criteria to calculate the baseline image as
p.addParameter('basePrctile' , 25);                                        % For the situation that 'prctile' is chosen, which percentile should the baseline image be calculated using
parse(p,varargin{:});
p = p.Results;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate the baseline

switch p.baseSelect
    case 'median'                                                          % If the median is selected...
        baseImg = median(fuObj.movie,3);                                   %  ...Calculate the median of each pixel
    case 'mode'                                                            % If the mode is selected...
        baseImg = halfSampleMode(reshape(fuObj.movie,[],fuObj.movieLen).');%  ...Calculate the half-sample mode for each pixel in the image
        baseImg = reshape(baseImg(:),fuObj.frameSize);                     %  ...Reshape the vector of modes into an image
    case 'mean'                                                            % If the mean is selected...
        baseImg = mean(fuObj.movie,3);                                     %  ...Calculate the mean of each pixel in the image
    case 'prctile'                                                         % If a percentile is selected
        baseImg = prctile(fuObj.movie,p.basePrctile,3);                    %  ...Calculate the percentile of each pixel's histogram1
    otherwise                                                              % The default is no correction so...
        baseImg = 0;                                                       %  ...Set the baseline image to nothing
end

end