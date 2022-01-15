function varargout = clusterTraces(uo, varargin)

% uo = clusterTraces(uo, varargin)
% 
% 
% 
% 2020 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

p = inputParser;                                                           % Set up an object to parse all of the various inputs
p.addParameter('numClusters'   , 20);                                      % Number of clusters
p.addParameter('motionCorrect' , false   );                                % Select if the motion corrected data should be used
p.addParameter('denoised'      , false   );                                % Choose to display denoised movie
p.addParameter('plotOpt'       , true    );                                % Choose to display denoised movie
p.addParameter('useMask'       , true    );
p.addParameter('meanSubtract'  , false   );
p.addParameter('normTraces'    , false   );
parse(p,varargin{:});
p = p.Results;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

[mov, frameSize] = loadMaskMovie(uo, p);                                   % Load movie as 3D array

mov       = reshape(mov, prod(frameSize), size(mov, 3));                   % Reshape movie to nPixel x nFrames
mov       = normalizeData(mov);                                            % Normalize the data
pixChoose = sum(abs(mov),2)>0;                                                  % Isolate nonzero pixels

traceSubset = mov(pixChoose,:);

idxClust = kmeans(traceSubset, p.numClusters, 'MaxIter',10000,...
                                        'Display','final','Replicates',1); % Do k-means

clustImg            = zeros(prod(frameSize),1);                            % Initialize the cluster map image
clustImg(pixChoose) = idxClust;                                            % Set the pixel values equal to the cluster numbers
clustImg            = reshape(clustImg, frameSize);                        % Reshape vector to be an image

if p.plotOpt
    figure(412)
    imagesc(clustImg)
    axis off
    pbaspect([uo.scanArea(2),uo.scanArea(1),1]);                           % Set the correct aspect ratio                              
    colormap(cat(1,[0,0,0],distinguishable_colors(p.numClusters,'k')))
end

if nargout > 0; varargout{1} = clustImg; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute correlation map

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extra functions


function [mov, frameSize]  = loadMaskMovie(uo, p)

mov             = uo.makeGetMovieFunction(p);                              % Get the movie given the requested options (denoised/motion corrected)
mov(isnan(mov)) = 0;                                                       % Remove nan values at the padded edges with zeros

if p.useMask;    uo.ensureMask();   end                                    % Make sure that a mask is available
if    p.motionCorrect 
    padVals = uo.post.motionPad;
    frameSize = uo.frameSize + sum(reshape(padVals,2,2),1);                % Get frame size
else
    padVals = [];
    frameSize = uo.frameSize;                                              % Get frame size
end
mov = applyMask(mov, uo.mask, padVals, p);                                 % Apply the pads if motion corrected

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mov = normalizeData(mov)

pixZeros  = sum(abs(mov),2)==0;                                            % Isolate nonzero pixels
mov       = bsxfun(@plus,  mov, -median(mov,1));                             % Remove the median
mov       = bsxfun(@times, mov, 1./sqrt(sum(mov,1)));                      % Normalize the data
% mov       = bsxfun(@plus,  mov, -mean(mov,1));                             % Remove the mean
% mov       = bsxfun(@times, mov, 1./std(mov,[],1));                         % Remove the variance
mov(pixZeros,:) = 0;                                                       % Set zero rows to nans

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%