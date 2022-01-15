function [Sm,Dm,RESULTS,varargout] = patchGraFT(data,K,patches,corr_kern, params)
% RUN_GRAFT_PATCHES - apply GraFT algorithm on overlapping patches in parallel
%
%   [Sm,Dm,RESULTS] = run_GraFT_patches(data,K,patches,corr_kern, params)
%
% Run the GraFT algorithm on a large dataset by operating on
% spatially overlapping patches in parallel and then merging the results.
% Processing in patches also allows the
% identification of weaker neurons without the need of normalization.

% TODO: The inputs is memory mapped, allowing for large datasets to be processed
% with reduced memory requirements. 
% TODO(?): The components are also classified by retaining only the components that
% correlate well with the raw data through classify_comp_corr.m
%
% INPUTS:
% data:    .mat file containing
%            data.Y      (the data matrix in the original dimensions)
%            data.Yr     (the data matrix reshaped in 2d format)
%            data.sizY   (dimensions of the original dataset)
%            data.nY     (minimum value of dataset)
%          OR the original dataset in 3d/4d format in which case the user
%          chooses whether to create a memory mapped file
% K:       number of components to be found in each patch
% patches: cell array containing the start and end points of each patch

%
% OUTPUTS:
% S:       3D array of spatial components
% D:       Matrix of temporal components
% RESULTS: Results of the on individual patches
%
% Adapted from: run_CNMF_patches
% Author: Eftychios A. Pnevmatikakis, Simons Foundation, 2015, 2016
% this version: Gal Mishne, 2019

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing/testing

if nargin < 2 || isempty(K); K = 10;             end                       % Set the default value for the sparsity per patch (K) at 10
params = checkAllParameters(params);
if params.verbose > 1;       params.verbose = 1; end                       % Make sure the verbosity level is at most one
params.n_dict = K;                                                         % Pass this value to the per-patch parameters
params.plot   = false;                                                     % Disable per-patch plotting
patchSizeSel  = extractPatchSize(params);                                  % Get the patch size from the parameters (and do error checking)

memmaped = isobject(data);
if ~memmaped
    Y      = data;
    clear data;  % TODO check if necessary
    sizY   = size(Y);                                                      % Get size of data
    Yr     = reshape(Y,[],sizY(end));
    F_dark = min(Yr(:));
    % create a memory mapped object named data_file.mat and open it read-only
    if params.create_memmap
        save('data_file.mat','Yr','Y','F_dark','sizY','-v7.3');
        data     = matfile('data_file.mat','Writable',false);
        memmaped = true;
    elseif isa(Yr, 'single') || isa(Yr, 'double')
        data = single(Yr);
    else
        data = single(Yr);
    end
else
    sizY = data.sizY;                                                      % Get size of data
end


if nargin < 3 || isempty(patches)
    fprintf('Patch construct not provided; constructing patches...')
    patches = construct_patches(sizY(1:end-1),patchSizeSel);             % TODO fix default for 3d case
    fprintf('done.\n')
end
n_patches = length(patches);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% running GraFT on each patch

RESULTS(n_patches) = struct('S', [],  'D', []);
if memmaped   %%%% TODO :: avoid copying the entire dataset to each worker, for in-memory data
    for i = 1:n_patches                
        RESULTS(i) = runMemmapPatch(data, F_dark, K, p, tau, options, ...
                                                               patches, i);% Run a single instance of memory mapped patch-wise processing
    end
else 
    for i = 1:n_patches
        RESULTS = runSinglePatch(Y, patches, corr_kern, params, RESULTS, i);% Run a single patch
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% combine results into one structure

[D, S, full_corr_kern] = combineDictionaryPatches(RESULTS, n_patches, patches, sizY, K);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% merge results

[Sm, Dm] = mergeWrapper(S, D, params, Yr, full_corr_kern, 2);              % Wrapper for merge code

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Final step: restimate spatial components given merged dictionary
% TODO: merge graphs from patch results instead of re-calculating graph

if isfield(params,'reEst')&&(params.reEst)
    fprintf('Final re-estimation of coefficient maps...')
    ckern = checkCorrKern(Y, corr_kern);                                   % Set up spatial smoothing kernel for full image
    Sm    = dictionaryRWL1SF(Y, Dm, ckern, params, full(Sm));              % Infer coefficients given the data and dictionary
    fprintf('done in %f min.\n', toc/60)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Some post-processing: reorder dictionary by magnitude

[Dm, Sm] = reorderDictionary(Y, Dm, Sm);                                   % Reorders te dictionary elements based on magnitude

if nargout > 3
    RESULTS2.S   = S;
    RESULTS2.D   = D;
    varargout{1} = RESULTS2;
    clear RESULTS2
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extra functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function to check all input parameters

function params = checkAllParameters(params)

    dParams.tau         = 1;                                               % Default tau values to be spatially varying
    dParams.lambda      = 0.6;                                             % Default lambda parameter is 0.6
    dParams.lamForb     = 0;                                               % Default Forbenius norm parameter is 0 (don't use)
    dParams.lamCont     = 0;                                               % Default Dictionary continuation term parameter is 0 (don't use)
    dParams.lamContStp  = 1;                                               % Default multiplicative change to continuation parameter is 1 (no change)
    dParams.lamCorr     = 0;                                               % Default Dictionary correlation regularization parameter is 0 (don't use)
    dParams.beta        = 0.09;                                            % Default beta parameter to 0.09
    dParams.maxiter     = 0.01;                                            % Default the maximum iteration to whenever Delta(Dictionary)<0.01
    dParams.numreps     = 2;                                               % Default number of repetitions for RWL1 is 2
    dParams.tolerance   = 1e-8;                                            % Default tolerance for TFOCS calls is 1e-8
    dParams.verbose     = 10;                                              % Default to full verbosity level
    dParams.likely_form = 'gaussian';                                      % Default to a gaussian likelihood ('gaussian' or 'poisson')
    dParams.step_s      = 1;                                               % Default step to reduce the step size over time (only needed for grad_type = 'norm')
    dParams.step_decay  = 0.995;                                           % Default step size decay (only needed for grad_type = 'norm')
    dParams.max_learn   = 1e3;                                             % Maximum number of steps in learning is 1000 
    dParams.learn_eps   = 1e-3;                                            % Default learning tolerance: stop when Delta(Dictionary)<0.01
%    dParams.n_dict      = selectDictSize(data_obj);                        % Default number of dictionary elements is a function of the data
    dParams.verb        = 1;                                               % Default to no verbose output
    dParams.grad_type   = 'full_ls_cor';                                   % Default to optimizing a full optimization on all dictionary elements at each iteration
    dParams.GD_iters    = 1;                                               % Default to one GD step per iteration
    dParams.bshow       = 0;                                               % Default to no plotting
    dParams.nneg_dict   = 0;                                               % Default to not having negativity constraints
    dParams.nonneg      = false;                                           % Default to not having negativity constraints on the coefficients
    dParams.plot        = false;                                           % Default to not plot spatial components during the learning
    dParams.updateEmbed = false;                                           % Default to not updateing the graph embedding based on changes to the coefficients
    dParams.mask        = [];                                              % for masked images (widefield data)
    dParams.normalizeSpatial = false;                                      % default behavior - time-traces are unit norm. when true, spatial maps normalized to max one and time-traces are not normalized
    dParams.randProj = false;
    dParams.mergeIters  = 3;
    dParams.reEst       = false;
    
    params = setParams(dParams, params);
    
    if ~isempty(params.mask)
        [params.nRows, params.nCols] = size(params.mask);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function to select patch sizes

function patchSizeSel = extractPatchSize(params)

if (~isfield(params, 'patchSize'))||(isempty(params.patchSize))            % Check that the patch size was selected
    params.patchSize  = [50,50];                                           % If not use a default of 50x50
end
if numel(params.patchSize) == 1
    patchSizeSel = params.patchSize*[1,1];                                 % Make sure the patchSizeSelis a 2-vector
elseif numel(params.patchSize) == 2
    patchSizeSel = params.patchSize;                                       % 
else
    error('Too many elements in the specified patch size!')
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Wrapper function to call GraFT on a single patch

function RESULTS = runSinglePatch(Y, patches, corr_kern, params, RESULTS, i)
fprintf(['Running GraFT on patch # ',num2str(i),' out of ',...
                                        num2str(length(patches)), '.\n']);
patchTime            = tic;                                                % Start the clock
patch_idx            = patch_to_indices(patches{i});                       % Get infices of current patch
Yp                   = Y(patch_idx{:},:);                                  % Extract the current patch
[D, S, ~, ~, K]      = GraFT(Yp, [], corr_kern, params);                   % Run GraFT on that patch
RESULTS(i).D         = D;                                                  % Save dictionary
RESULTS(i).S         = S;                                                  % Save the spatial maps
RESULTS(i).patch_idx = i;        
RESULTS(i).K         = K;
% Keep track of the actual patch name
Tend                 = toc(patchTime);                                     % Read out the time
fprintf(['Finished processing patch # ',num2str(i),' out of ', ...
                          num2str(length(patches)), ' in %f s.\n'], Tend);
clear patchTime
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function idx = patch_to_indices(patch)
    % helper function to build indices vector from patch start/stop indices
    idx = arrayfun(@(x,y) x:y, patch(1:2:end), patch(2:2:end), 'un', false);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function to process a single patch object

function CNM = process_patch_object(Y,F_dark,K,p,tau,options)
    CNM = CNMF();
    if ndims(Y) > 3; d3 = size(Y,3); else; d3 = 1; end
    options = CNMFSetParms(options,...
                'd1',size(Y,1),...
                'd2',size(Y,2),...
                'd3',d3,...
                'p',p,...
                'gSig',tau,...
                'temporal_parallel',false,...
                'spatial_parallel',false,...
                'space_thresh',options.patch_space_thresh,...
                'time_thresh',options.patch_time_thresh,...
                'cnn_thr',options.patch_cnn_thr,...
                'min_fitness',options.patch_min_fitness);
    Y           = single(Y) - single(F_dark);
    Y(isnan(Y)) = single(F_dark);
    CNM.fit(Y,options,K);                
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Wrapper for merge code to merge components

function [Sm, Dm] = mergeWrapper(S, D, params, Yr, full_corr_kern, mergeOpt)

fprintf('Running full GraFT to merge components...\n')

tic
if mergeOpt == 1
    Dm = D;
    Sm = S;

    for i = 1:params.mergeIters
        fprintf('Starting iteration %d...\n',i);tic;
        params.verbose = 10;
        gcp
        [Sm, W] = dictionaryRWL1SF(Yr,Dm,full_corr_kern,params,Sm);   % Infer coefficients given the data and dictionary
        fprintf('...Updated coefficients in %f seconds.\n', toc); tic
        Dm = dictionary_update(Yr.',Dm, Sm.',params.step_s,params);  % Take a gradient step with respect to the dictionary
        fprintf('...Updated dictionary in %f seconds.\n', toc);
     end
     if ~params.normalizeSpatial
        if isfield(params,'likely_form')&&strcmp(params.likely_form,'homotopy')
            [Sm, ~] = IRWH(Yr, Dm, full_corr_kern, params, W);
        else
            [Sm, ~] = dictionaryRWL1SF(Yr,Dm,full_corr_kern,params,Sm);   % Infer coefficients given the data and dictionary
        end
    end
else
%     Dm = D;
%     Sm = S;
%     
%     for i = 1:params.mergeIters
%         [Dm,Sm] = mergeGraFTdictionaries(Dm,Sm,[],params.normalizeSpatial);
%     end
%     Sm      = reCalcCoefWithLS(Yr, Dm, Sm);
    [Dm,Sm] = mergeGraFTdictionaries(D,S,[],params.normalizeSpatial);
    Sm      = reCalcCoefSparse(Yr, Dm, Sm, params.lambda, true, 'lasso');      % Recalculate the coefficients with the LASSO
end


Dnorms   = sqrt(sum(Dm.^2,1));                                             % Get norms of each dictionary element
Smax     = max(Sm,[],1);                                                    % Get maximum value of each spatial map
actMeas  = Dnorms(:).*Smax(:);                                             % Total activity metric is the is the product of the above
[~,IX]   = sort(actMeas,'descend');                                        % Get the indices of the activity metrics in descending order
Dm       = Dm(:,IX);                                                       % Reorder the dictionary
Sm       = Sm(:,IX);                                                        % Reorder the spatial maps

toc
fprintf(' done. \n');

end
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function to re-order the dictionary to me more sensible

function [Dm, Sm] = reorderDictionary(Y, Dm, Sm)

Dnorms   = sqrt(sum(Dm.^2,1));                                             % Get norms of each dictionary element
Smax     = max(Sm,[],1);                                                   % Get maximum value of each spatial map
actMeas  = Dnorms(:).*Smax(:);                                             % Total activity metric is the is the product of the above
[~,IX]   = sort(actMeas,'descend');                                        % Get the indices of the activity metrics in descending order
Dm       = Dm(:,IX);                                                       % Reorder the dictionary
Sm       = Sm(:,IX);                                                       % Reorder the spatial maps
Sm       = reshape(full(Sm),size(Y,1),size(Y,2),[]);                       % Reshape presence maps to be an image stack

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

function RESULTS = runMemmapPatch(data, F_dark, K, p, tau, options, patches, i)

fprintf(['Starting processing patch # ',num2str(i),' out of ', ...
                                        num2str(length(patches)), '.\n']);
patch_idx = patch_to_indices(patches{i});
Yp        = data.Y(patch_idx{:},:);
RESULTS   = process_patch_object(Yp,F_dark, K, p, tau, options);
fprintf(['Finished processing patch # ',num2str(i),' out of ',...
                                        num2str(length(patches)), '.\n']);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
function idx = patch_to_linear(patch, sizY)
    % helper function to build linear indices from patch start/stop indices
    slice_idx = patch_to_indices(patch);
    subs_idx = cell(1, numel(slice_idx));
    [subs_idx{:}] = ndgrid(slice_idx{:});
    subs_idx = cellfun(@(x) x(:), subs_idx, 'un', false);
    idx = sub2ind(sizY(1:end-1), subs_idx{:});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

function result = process_patch(Y, F_dark, K, p, tau, options)
    % helper function to apply CNMF to a small patch
    sizY       = size(Y);
    options.d1 = sizY(1);
    options.d2 = sizY(2);
    if ndims(Y) == 3    options.d3 = 1;
    else                options.d3 = sizY(3);
    end
    options.nb = 1;
    options.temporal_parallel = 0;  % turn off parallel updating for temporal components
    options.spatial_parallel  = 0;   % turn off parallel updating for spatial components
    options.space_thresh      = options.patch_space_thresh;    % put a low acceptance threshold initially
    options.time_thresh       = options.patch_time_thresh;

    Y           = double(Y - F_dark);
    Y(isnan(Y)) = F_dark;
    [P,Y]       = preprocess_data(Y,p,options);
    Yr          = reshape(Y,[],sizY(end));

    [Ain,Cin,bin,fin] = initialize_components(Y,K,tau,options,P);
    [A,b,Cin,P]       = update_spatial_components(Yr,Cin,fin,[Ain,bin],P,options);
    P.p               = 0;
    options.p         = 0;
    [C,f,P,S,YrA]     = update_temporal_components(Yr,A,b,Cin,fin,P,options);

    if ~isempty(A) && ~isempty(C)
        [Am,Cm,~,~,P] = merge_components(Yr,A,b,C,f,P,S,options);
        [A,b,Cm,P]    = update_spatial_components(Yr,Cm,f,[Am,b],P,options);
        [C,f,P,S,YrA] = update_temporal_components(Yr,A,b,Cm,f,P,options);
        [rval_space,rval_time,ind_space,ind_time] = classify_comp_corr(Y,A,C,b,f,options); 
        %ind = ind_space & ind_time;        
        ind_corr = ind_space;
        
        try  % matlab 2017b or later is needed
            [ind_cnn,value] = cnn_classifier(A,[options.d1,options.d2],'cnn_model',options.cnn_thr);
        catch
            ind_cnn = true(size(A,2),1);                        % components that pass the CNN classifier
        end     

        fitness = compute_event_exceptionality(C+YrA,options.N_samples_exc,options.robust_std);
        ind_exc = (fitness < options.min_fitness);
        ind     = (ind_corr | ind_cnn) & ind_exc;
        %fitness_delta = compute_event_exceptionality(diff(C+YrA,[],2),0);
        %ind = (ind_space & ind_time) | (fitness < options.patch_max_fit) | (fitness_delta < options.patch_max_fit_delta);
        
        P.rval_space    = rval_space;
        P.rval_time     = rval_time;
        P.ind_space     = ind_space;
        P.ind_time      = ind_time;
        P.fitness       = fitness;
        P.fitness_delta = fitness_delta;
        P.A_throw       = A(:,~ind);
        P.C_throw       = C(~ind,:);
    end
    result.A = A(:,ind);
    result.b = b;
    result.C = C(ind,:);
    result.f = f;
    result.S = S;
    result.P = P;
end

%}

%% EOF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
