function displayExampleFrame(uo, varargin)

% displayMovie(fuObj, varargin)
%
% Function to display a movie
%
% 2019 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

p = inputParser;                                                           % Set up an object to parse all of the various inputs
p.addParameter('baseSelect'  , 'median');                                  % Select what criteria to calculate the baseline image as
p.addParameter('center    '  , true    );                                  % Select whether to center the data (to the temporal per-pixel baseline)
p.addParameter('normalize'   , true    );                                  % Select whether to normalize the data (to the temporal per-pixel baseline)
p.addParameter('basePrctile' , 25      );                                  % For the situation that 'prctile' is chosen, which percentile should the baseline image be calculated using
p.addParameter('frameNo'     , 'active');                                  % Select which frame display
p.addParameter('scaleFactor' , [2,10]  );                                  % Select scaling factor for the image
p.addParameter('figNo'       , 104     );                                  % Select figure number to plot to
p.addParameter('denoised'    , false   );                                  % Choose to display denoised movie

parse(p,varargin{:});
p = p.Results;

if p.denoised&&(~doesDenoiseExist(uo))
    warning('Requested denoised movie for PCA, but no denoised movie available. Running wavelet denoising with basic features...')
    uo.denoiseTracesWavelet();
end

getFrame  = makeGetBlockFunction(uo,p);                                    % Make an anonymous function that loads frames
baseImg   = getBaseImage(getFrame(1:uo.movieLen), uo, p);                  % Get base image
movNorms  = sum(sum(getFrame(1:uo.movieLen).^2,1),2);                      % Get frame norms
p.frameNo = selectFrameNo(uo, movNorms, p.frameNo);                        % Choose which frames to display

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Normalize the data to plot

frameToDisp = getFrame(p.frameNo);                                         % Grab the appropriate frame
if p.center;    frameToDisp = frameToDisp  - baseImg; end                  % Subtract baseline values if required
if p.normalize; frameToDisp = frameToDisp./baseImg;   end                  % Normalize the values if required
frameMode   = halfSampleMode(frameToDisp(:));                              % Get the mode of the current frame
[~,robSTD2] = robustTwoSidedSTD(frameToDisp(:));                           % Compute the robust standard deviation: needed to scale the data for appropriate contrast
colorLims   = p.scaleFactor.*robSTD2;                                      % Select the limits that will guide the contrast (proportional to robust STD)
frameToDisp = softScale(frameToDisp, colorLims, frameMode);                % Softscale the data with a logistic to get rid of large values without severe (ugly) clipping

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Display the frame

figure(p.figNo)                                                            % Open a new figure
imagesc(frameToDisp,[0,1])                                                 % Display the frame
pbaspect([uo.scanArea(2),uo.scanArea(1),1]);                               % Set the correct aspect ratio
axis off; colormap gray;                                                   % Remove axes and set colormap
title(sprintf('Frame number %d', p.frameNo))                               % Give a title to the image

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function frameNo = selectFrameNo(uo, movNorms, frameNo)

if isempty(frameNo)
    frameNo = 'active';
end

if isnumeric(frameNo) % don't need to do anythin
elseif ischar(frameNo)
    switch frameNo
        case 'rand'
            frameNo = randsample(uo.movieLen, 1);
        case 'active'
            movNorms = find(movNorms>prctile(movNorms(:), 90));
            frameNo  = movNorms(randsample(numel(movNorms), 1));
        otherwise
            warning('unknown option given. Choosing an active frame at random\n')
            frameNo = selectFrameNo(uo, 'active');
    end
else
    error('Bad frame number selection given')
end
end

%% 

function baseImg = getBaseImage(mov, uo, p)

switch p.baseSelect
    case 'median'                                                          % If the median is selected...
        if isempty(uo.stats.medImg); uo = uo.calcBasicStats(); end         %  ...Ensure the median image exists
        baseImg = uo.stats.medImg;                                         %  ...Calculate the median of each pixel
    case 'mode'                                                            % If the mode is selected...
        baseImg = halfSampleMode(reshape(mov,[],uo.movieLen).');           %  ...Calculate the half-sample mode for each pixel in the image
        baseImg = reshape(baseImg(:),uo.frameSize);                        %  ...Reshape the vector of modes into an image
    case 'mean'                                                            % If the mean is selected...
        if isempty(uo.stats.medImg); uo = uo.calcBasicStats(); end         %  ...Ensure the mean image exists
        baseImg = uo.stats.meanImg;                                        %  ...Calculate the mean of each pixel in the image
    case 'prctile'                                                         % If a percentile is selected
        baseImg = prctile(mov,p.basePrctile,3);                            %  ...Calculate the percentile of each pixel's histogram1
    otherwise                                                              % The default is no correction so...
        baseImg = 0;                                                       %  ...Set the baseline image to nothing
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%