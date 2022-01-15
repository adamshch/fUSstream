function uo = graphMap(uo, varargin)

% uo = laplaceDiffMap(uo, varargin)
%
% Fuction to create graph laplace maps
%
% 2020 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;                                                           % Set up an object to parse all of the various inputs
p.addParameter('motionCorrect' , false   );                                % Select if the motion corrected data should be used
p.addParameter('denoised'      , false   );                                % Choose to display denoised movie
p.addParameter('plotOpt'       , false   );                                % Optional plotting flag
p.addParameter('figNumber'     , 203     );                                % Optional figure number to plot to
p.addParameter('Knn'           , 15      );                                % Choose number of nearest neighbors
p.addParameter('useMask'       , true    );                                % Choose to use the mask and only consider "inside" pixels
p.addParameter('distSelect'    , 'emd'   );                                % Choose the distance metric
parse(p,varargin{:});
p = p.Results;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Some initializations

if p.useMask;    uo.ensureMask();   end                                    % Make sure that a mask is available
movieOut = uo.makeGetMovieFunction(p);                                     % Get the movie given the requested options (denoised/motion corrected)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute graph embedding

switch lower(p.graph)
    case 'laplace'
        
    case {'diffusion','diff'}
        
    otherwise
        error('Unknown graph creation')
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pixCov = createSparseCov(uo, movieOut, Knn)

covInd  = [];
covVals = [];
for ll = 1:prod(uo.frameSize)
    [Ix,Jx] = ind2sub(uo.frameSize, ll);                                   % 
    corrVec = sum(bsxfun(@times, movieOut, movieOut(Ix,Jx,:)));            % 
    [v,ix]  = find(corrVec > pctile(corrVec, 100*(1-Knn/numel(corrVec)))); % Find the top Knn values
    covInd  = cat(1, covInd, vec(prod(uo.frameSize)*(ll-1) + v));          % 
    covVals = cat(1, covVals, ix(:));                                      % 
end
[Ix,Jx] = ind2sub(uo.frameSize, covInd);
pixCov = sparse(Ix, Jx, covVals, prod(uo.frameSize), prod(uo.frameSize));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

