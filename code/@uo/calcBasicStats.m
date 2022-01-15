function uo = calcBasicStats(uo, varargin)

% uo = calcBasicStats(uo, varargin)
%
% Function to compute some basic statistics for the data
%
% 2020 - Adam Charles


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing 

p = inputParser;                                                           % Set up an object to parse all of the various inputs
p.addParameter('calcMedImg'  , true);                                      % Choose to calculate the mean image
p.addParameter('calcMeanImg' , true);                                      % Choose if to load the entire data to RAM or to just link to the MAT file
p.addParameter('calcMaxImg'  , true);                                      % Need to input either the scan area...
p.addParameter('calcMinImg'  , true);                                      % ...or scan resolution
p.addParameter('calcVarImg'  , true);                                      % ...and/or the times of frames.
parse(p,varargin{:});
p = p.Results;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute the stats

uo = allBasicStats(uo);
uo = computeMedianImage(uo,p);                                             % Compute the median image (if requested)
uo = computeMeanImage(uo,p);                                               % Compute the mean image (if requested)
uo = computeMinImage(uo,p);                                                % Compute the max image (if requested)
uo = computeMaxImage(uo,p);                                                % Compute the min image (if requested)
uo = computeVarImage(uo,p);                                                % Compute the variance image (if requested)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extra functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute all basic stats
function uo = allBasicStats(uo)

if ~uo.isMatFileBased()
    uo.stats.median = median(uo.movie(:));
    uo.stats.mean   = mean(uo.movie(:)  );
    uo.stats.max    = max(uo.movie(:)   );
    uo.stats.min    = min(uo.movie(:)   );
    uo.stats.var    = var(uo.movie(:)   );
else
    uo.stats.median = median(vec(uo.meta.movMatFile.(uo.meta.movVarName)));
    uo.stats.mean   = mean(vec(uo.meta.movMatFile.(uo.meta.movVarName)  ));
    uo.stats.max    = max(vec(uo.meta.movMatFile.(uo.meta.movVarName)   ));
    uo.stats.min    = min(vec(uo.meta.movMatFile.(uo.meta.movVarName)   ));
    uo.stats.var    = var(vec(uo.meta.movMatFile.(uo.meta.movVarName)   ));
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function to compute the median image
function uo = computeMedianImage(uo,p)
if (isempty(uo.stats.medImg))&&(p.calcMedImg)
    if ~uo.isMatFileBased()
        uo.stats.medImg = median(uo.movie,3);
    else
        uo.stats.medImg = median(uo.meta.movMatFile.(uo.meta.movVarName),3);
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function to compute the mean image
function uo = computeMeanImage(uo,p)
if isempty(uo.stats.meanImg)&&(p.calcMeanImg)
    if ~uo.isMatFileBased()
        uo.stats.meanImg = mean(uo.movie,3);
    else
        uo.stats.meanImg = mean(uo.meta.movMatFile.(uo.meta.movVarName),3);
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function to compute the max image
function uo = computeMaxImage(uo,p)
if isempty(uo.stats.maxImg)&&(p.calcMaxImg)
    if ~uo.isMatFileBased()
        uo.stats.maxImg = max(uo.movie,[],3);
    else
        uo.stats.maxImg = max(uo.meta.movMatFile.(uo.meta.movVarName),...
                                                                     [],3);
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function to compute the min image
function uo = computeMinImage(uo,p)
if isempty(uo.stats.minImg)&&(p.calcMinImg)
    if ~uo.isMatFileBased()
        uo.stats.minImg = min(uo.movie,[],3);
    else
        uo.stats.minImg = min(uo.meta.movMatFile.(uo.meta.movVarName),...
                                                                     [],3);
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function to compute the variance image
function uo = computeVarImage(uo,p)
if isempty(uo.stats.varImg)&&(p.calcVarImg)
    if ~uo.isMatFileBased()
        uo.stats.varImg = var(uo.movie,[],3);
    else
        uo.stats.varImg = var(uo.meta.movMatFile.(uo.meta.movVarName),...
                                                                     [],3);
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

