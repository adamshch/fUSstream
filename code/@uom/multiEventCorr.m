function eSTH = multiEventCorr(uom, varargin)

% uom = multiEventCorr(uom, varargin)
%
%
%
% 2020 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

p = inputParser;                                                           % Set up an object to parse all of the various inputs
p.addParameter('kernel'      , 'asymgauss');                               % 
p.addParameter('kernParams'  , {4,[1e-10, 4],[4, 1e-10]});                 % 
p.addParameter('refImage'    , []);                                        % Option to input a user-defined reference image (useful for multi-dataset computations)
p.addParameter('medBaseLine' , 'firstmedian');                             % Select the reference image
p.addParameter('evtName'     , 'sync_cin');                               % Select event to correlate to
p.addParameter('useDenoised' , false);                                     % Option to motion correct based on the denoised data
p.addParameter('useMask'     , true);                                      % Option to use the mask in order to compute shifts only based on brain motion (ignore out-of-brain motion)                           
parse(p,varargin{:});
p = p.Results;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%


evtData = cell(numel(uom.uo),1);
for ll = 1:numel(uom.uo)                                                   % Iterate over all the data objects
    fprintf('Getting event data for session %d of %d:\n', ll, numel(uom.uo))
    TMPfile = load(fullfile(uom.uo{ll}.meta.dataPath,[uom.uo{ll}.meta.fileName,'_meta.mat']));
    try          evtData{ll} = [TMPfile.ratdata.(p.evtName)];
    catch;       error('Unknown meta-data to test against!')
    end
    evtData{ll} = round(evtData{ll}./(uom.uo{ll}.dt));
    fprintf('done.\n')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

if iscell(p.kernParams)
    eSTH.ind = cell(numel(uom.uo), numel(p.kernParams));
    eSTH.all = cell(1, numel(p.kernParams));
    for kk = 1:numel(p.kernParams)
        eSTH.all{kk} = 0;
    end
else
    eSTH.ind = cell(numel(uom.uo),1);
    eSTH.all = 0;
end

for ll = 1:numel(uom.uo)                                                   % Iterate over all the data objects
    fprintf('Getting eSTH %d of %d:\n', ll, numel(uom.uo))
    [uoMask, uoMean] = getMotSizedMaskAndMean(uom.uo{ll});
    if iscell(p.kernParams)
        for kk = 1:numel(p.kernParams)
            if numel(p.kernParams{kk})==1
                kerOpt = 'gauss';
            else
                kerOpt = 'asymgauss';
            end
            
            eSTH.ind{ll,kk}  = uom.uo{ll}.eventSTH('eventVec', evtData{ll},'kernel', kerOpt, 'kernParams', p.kernParams{kk}, 'motionCorrect' , true);
            if p.useMask; eSTH.ind{ll,kk} = uoMask.*(eSTH.ind{ll,kk}-uoMean);
            else;         eSTH.ind{ll,kk} = (eSTH.ind{ll,kk}-uoMean);
            end
%             eSTH.all{kk}     = eSTH.all{kk} + eSTH.ind{ll,kk};
        end
    else
        if numel(p.kernParams)==1
            kerOpt = 'gauss';
        else
            kerOpt = 'asymgauss';
        end
        eSTH.ind{ll}    = uom.uo{ll}.eventSTH('eventVec', evtData{ll},'kernel', kerOpt, 'kernParams', p.kernParams, 'motionCorrect' , true);
        eSTH.ind{ll,kk} = uoMask.*(eSTH.ind{ll}-uoMean);
%         eSTH.all        = eSTH.all + eSTH.ind{ll};
    end
    fprintf('done.\n')
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function [uoMask, uoMean] = getMotSizedMaskAndMean(uo)

padSize = uo.post.motionPad;
uoMask  = padarray(padarray(uo.mask, padSize([1,3]),0,'pre'), ...
                                              padSize([2,4]), 0, 'post');
                                          
% p.motionCorrect = 'true';
% movieOut = uo.makeGetMovieFunction(p);
% uoMean   = mean(movieOut,3);
uoMean  = 0*padarray(padarray(uo.stats.medImg, padSize([1,3]),0,'pre'), ...
                                              padSize([2,4]), 0, 'post');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%