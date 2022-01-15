function fuObj = fuGraFT(fuObj, varargin)

% fu = fuPCA(fuObj)
% 
% Function to apply GraFT within the functional ultrasound struct framework
% 
% 2021 - Ruolan Sun

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing
p = inputParser;                                                           % Set up an object to parse all of the various inputs
p.addParameter('smoothOpt'      , false);                                  % Select whether to use the smoothed traces
p.addParameter('lambda'         , 1.5);                                  % Sparsity parameter
p.addParameter('tau'            , 2);                                      % tau parameter in the weight update tau/(beta + ...)
p.addParameter('n_dict'         , 30);                                     % Number of dictionary elements to learn
p.addParameter('numreps'        , 2);                                      % Default number of repetitions for RWL1 is 2
p.addParameter('beta'           , 0.01);                                   % Beta parameter in the weight update tau/(beta + ...)
p.addParameter('learn_eps'      , 1e-3);                                   % Tolerance parameter in the learning
p.addParameter('max_learn'      , 20);                                     % Maximum number of iterations to run in the learning (weight and dictionary updates)
p.addParameter('grad_type'      , 'full_ls_cor');                          % type of dictionary update
p.addParameter('nonneg'         , true);                                   % Set presence values to be nonnegative
p.addParameter('nneg_dict'      , true);                                   % Set dictionary elements to be nonnegative
p.addParameter('lamForb'        , 0.3);                                    % Parameter to control how much to weigh extra time-traces
p.addParameter('lamCont'        , 0.4);                                    % Parameter to control how much to weigh the previous estimate (continuity term)
p.addParameter('lamCorr'        , 0.3);                                    % Parameter to prevent overly correlated dictionary elements 
p.addParameter('lamContStp'     , 0.9);                                    % Decay rate of the continuation parameter
p.addParameter('plot'           , true);                                   % Set whether to plot intermediary variables
p.addParameter('normalizeSpatial'      , false);                           % Select sparsity constraint for RPCA

corr_kern.w_time   = 0;                                                    % Initialize the correlation kernel struct
corr_kern.corrType = 'convolution';                                        % Set the correlation type to "graph embedding"

parse(p,varargin{:});
p = p.Results;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get the movie to operate on

mov = extractMovie(fuObj,p);                                               % Extract the appropriate version of the movie

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Normalize the data

[r,c,~]=size(mov);
data_min=zeros(r,c); data_max=zeros(r,c);
for i=1:1:r
    for j=1:1:c
        data_min(i,j)=min(mov(i,j,:));
        data_max(i,j)=max(mov(i,j,:));
    end
end
sim_cube=(mov-data_min)./(data_max-data_min);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize the variables needed (embedding kernel & Initialization)
% While GraFT's main function will compute these internally if not
% calculated here, pre-calculating can make running on multiple parameters
% faster, and keep things consistant with using the same (random
% initialization.

ckern    = checkCorrKern(sim_cube, corr_kern, 10);                         % Set up spatial smoothing kernel (odd frames)
dictInit = dictInitialize([], size(sim_cube,3), p.n_dict);                 % Initialized the dictionary

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Actually run GraFT
% This is how to run the main function. The commented line is how to call
% GraFT if you don't want to run the above lines independently

[D, S] = GraFT(sim_cube, dictInit, ckern, p);                              % Learn the dictionary

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up output

GraFTout.spatial = S;
GraFTout.TimeTrace = D;
fuObj.GraFTout = GraFTout;

%% Plot the results

figure;
imagesc(basis2img2(reshape(S,[], size(S,3)), ...
                             size(S(:,:,1)), ceil(sqrt(size(S,3)))*[1,1])) % Plot learned presence maps
axis image; axis off; set(gcf, 'color', [1,1,1])                           % Set up the aesthetics of the plot
title('Learned Presence Maps')

figure; 
plot(D)                                                                    % Plot learned time traces
box off; set(gca, 'TickDir', 'out'); set(gcf, 'color', [1,1,1]);           % Set up the aesthetics of the plot
title('Learned Time Traces')
xlabel('Time (frmaes)')
ylabel('F (AU)')

figure,
imagesc(...
            diag(1./sqrt(sum(D.^2,1)))*(D.'*D)*diag(1./sqrt(sum(D.^2,1)))) %
axis image; axis off; set(gcf, 'color', [1,1,1])
title('Correlation matrix: learned v. learned')
xlabel('Learned component index (sorted)')
ylabel('Learned component index (sorted)')

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
function mov = extractMovie(uo,p)

if ~isMatFileBased(uo)
    if p.smoothOpt
        if ~doesDenoiseExist(uo)% ~(isfield(uo.post,'movDenoised')&&(numel(uo.post.movDenoised)==numel(uo.movie)))
            warning('Requested denoised movie for PCA, but no denoised movie available. Running wavelet denoising with basic features...')
            uo.denoiseTracesWavelet();
        end
        mov = uo.post.movDenoised;                                         % Extract the denoised video
    else
        mov = uo.movie;                                                    % Extract the original movie
        mov = fillBurstErrors(uo, mov);                                    % Correct burst errors if requested
    end
else
    if p.smoothOpt
        if ~doesDenoiseExist(uo)
            warning('Requested denoised movie for PCA, but no denoised movie available. Running wavelet denoising with basic features...')
            uo.denoiseTracesWavelet();
        end
        mov       = uo.post.movDenoised.uoDenoised;
    else
        mov       = uo.meta.movMatFile.(uo.meta.movVarName);
        mov = fillBurstErrors(uo, mov);                                    % Correct burst errors if requested
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%