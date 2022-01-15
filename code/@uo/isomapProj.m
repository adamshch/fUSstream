function varargout = isomapProj(uo, varargin)

% function uo = correlationMap(uo, varargin)
% 
% Function that seeks out and quantifies residual sub-pixel motion errors
% in functional ultrasound imaging. 
% 
% 2020 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

p = inputParser;                                                           % Set up an object to parse all of the various inputs
p.addParameter('neighborFcn' , 'epsilon');                                 % 
p.addParameter('neighborSize', []);                                        % 
p.addParameter('numProj'     , 3);                                         % 
p.addParameter('distance'    , 'euclidean');                               % The distance measure with which to compute the pariwise distances between points with
parse(p,varargin{:});
p = p.Results;

p.neighborFcn = lower(p.neighborFcn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up ISOMAP parameters

movieCopy = makeGetMovieFunction(uo,p);

ISOMAPopts.dims    = p.numProj;
ISOMAPopts.display = 0;
ISOMAPopts.overlay = 0;
ISOMAPopts.verbose = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run ISOMAP

tmp = reshape(movieCopy, [],uo.movieLen);

distMat      = squareform(pdist(tmp,p.distance));                          % Compute pairwise distances between points
clear tmp getMovie                                                         % Clear out the 
neighborSize = selectNeighborSize(distMat,p);                              % Select the neighborhood size (if not provided)
distMat      = double(distMat);
[Y, ~, ~]    = Isomap(distMat, p.neighborFcn, neighborSize, ISOMAPopts);   % Run ISOMAP


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

if ISOMAPopts.dims == 3
    isomapImg = Y.coords;
end

varargout{1} = isomapImg;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function neighborSize = selectNeighborSize(distMat,p)
switch p.neighborFcn
    case 'k'
        if isempty(p.neighborSize)
            neighborSize = 10;
        elseif isscalar(p.neighborSize)
            if isinteger(p.neighborSize)
                neighborSize = p.neighborSize;
            else
                error('k-neighborhoods requires an integer value for k.')
            end
        else
            error('Bad input for neighborhood size')
        end
    case 'euclidean'
        if isempty(p.neighborSize)
            neighborSize = prctile(vec(distMat),5);
        elseif isscalar(p.neighborSize)
            neighborSize = p.neighborSize;
        else
            error('Bad input for neighborhood size')
        end
    otherwise
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
