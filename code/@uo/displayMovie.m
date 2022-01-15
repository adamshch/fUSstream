function displayMovie(uo, varargin)

% displayMovie(fuObj, varargin)
%
% Function to display a movie. Based off of Sue Ann Koay's amazing
% "MovieSlider" package.
%
% 2019 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

p = inputParser;                                                           % Set up an object to parse all of the various inputs
p.addParameter('baseSelect'    , 'median');                                % Select what criteria to calculate the baseline image as
p.addParameter('basePrctile'   , 25);                                      % For the situation that 'prctile' is chosen, which percentile should the baseline image be calculated using
p.addParameter('normSelect'    , 'base');                                  % Select what criteria to calculate the normalization image as
p.addParameter('normPrctile'   , 25);                                      % For the situation that 'prctile' is chosen, which percentile should the baseline image be calculated using
p.addParameter('denoised'      , false);                                   % Choose to display denoised movie
p.addParameter('removeBursts'  , true);                                    % Choose to remove burst errors
p.addParameter('motionCorrect' , false);                                   % Choose to display motion corrected movie
parse(p,varargin{:});
p = p.Results;

if strcmp(p.normSelect,'base')
    p.normSelect  = p.baseSelect;
    p.normPrctile = p.basePrctile;
end

if p.denoised; uo.ensureDenoisedData(); end                                % Ensure that the denoised data exists

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate the baseline

movieOut           = uo.makeGetMovieFunction(p);
[baseImg, normImg] = createBaseline(movieOut, uo, p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Display the movie

if ~p.denoised
    if isfield(uo.errs,'burst')&&(numel(uo.errs.burst)==size(uo.movie,ndims(uo.movie)))
        uo.findBurstFrames();
        errFrames = p.removeBursts*uo.errs.burst;
    else
        warning('No burst error detection output found. Running findBurstFrames with default settings...\n')
        uo.findBurstFrames();
        errFrames = p.removeBursts*uo.errs.burst;
    end
    mov = MovieSlider(bsxfun(@times, reshape(errFrames==0,[1,1,numel(uo.errs.burst)]), ...
                        bsxfun(@times,bsxfun(@plus, movieOut, -baseImg), 1./normImg)));
else
    mov = MovieSlider(bsxfun(@times,bsxfun(@plus, movieOut, -baseImg), 1./normImg));
end
aspRatio = [uo.frameSize(2)./uo.scanArea(2), ...
                                    uo.frameSize(1)./uo.scanArea(1), 1];   % Set aspect ratio for movie
daspect(mov.axsMovie, aspRatio)                                            % Set the movie's aspect ratio

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function to generate the per-pixel baseline and normalization

function [baseImg, normImg] = createBaseline(mov, uo, p)

switch p.baseSelect
    case 'median'                                                          % If the median is selected...
        baseImg = median(mov,3);                                           %  ...Calculate the median of each pixel
    case 'mode'                                                            % If the mode is selected...
        baseImg = halfSampleMode(reshape(mov, [], uo.movieLen).');         %  ...Calculate the half-sample mode for each pixel in the image
        baseImg = reshape(baseImg(:),uo.frameSize);                        %  ...Reshape the vector of modes into an image
    case 'mean'                                                            % If the mean is selected...
        baseImg = mean(mov,3);                                             %  ...Calculate the mean of each pixel in the image
    case 'prctile'                                                         % If a percentile is selected
        baseImg = prctile(mov,p.basePrctile,3);                            %  ...Calculate the percentile of each pixel's histogram1
    otherwise                                                              % The default is no correction so...
        baseImg = 0;                                                       %  ...Set the baseline image to nothing
end

switch p.normSelect
    case 'median'                                                          % If the median is selected...
        normImg = median(mov,3);                                           %  ...Calculate the median of each pixel
    case 'mode'                                                            % If the mode is selected...
        normImg = halfSampleMode(reshape(mov,[], uo.movieLen).');          %  ...Calculate the half-sample mode for each pixel in the image
        normImg = reshape(normImg(:),uo.frameSize);                        %  ...Reshape the vector of modes into an image
    case 'mean'                                                            % If the mean is selected...
        normImg = mean(mov,3);                                             %  ...Calculate the mean of each pixel in the image
    case 'prctile'                                                         % If a percentile is selected
        normImg = prctile(mov,p.normPrctile,3);                            %  ...Calculate the percentile of each pixel's histogram1
    otherwise                                                              % The default is no correction so...
        normImg = 1;                                                       %  ...Set the baseline image to nothing
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
