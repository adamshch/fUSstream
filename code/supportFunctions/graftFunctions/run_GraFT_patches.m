function [Sm,Dm,RESULTS,varargout] = run_GraFT_patches(data,K,patches,corr_kern, params)
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

params.n_dict = K;

memmaped = isobject(data);
if ~memmaped
    Y      = data;
    clear data;  % TODO check if necessary
    sizY   = size(Y);
    Yr     = reshape(Y,[],sizY(end));
    F_dark = min(Yr(:));

    % create a memory mapped object named data_file.mat and open it read-only
    if params.create_memmap
        save('data_file.mat','Yr','Y','F_dark','sizY','-v7.3');
        data = matfile('data_file.mat','Writable',false);
        memmaped = true;
    elseif isa(Yr, 'single') || isa(Yr, 'double')
        data = Yr;
    else
        data = single(Yr);
    end
else
    sizY = data.sizY;
end

if nargin < 2 || isempty(K)
    K = 10; 
end

if nargin < 3 || isempty(patches)
    fprintf('Patch construct not provided; constructing patches...')
    patches = construct_patches(sizY(1:end-1),[50,50]);  % TODO fix default for 3d case
    fprintf('done.\n')
end
n_patches = length(patches);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% running GraFT on each patch, in parallel

RESULTS(n_patches) = struct('S', [],  'D', []);
tic
if memmaped   %%%% TODO 
    parfor i = 1:n_patches
        fprintf(['Starting processing patch # ',num2str(i),' out of ',num2str(n_patches), '.\n']);
        patch_idx  = patch_to_indices(patches{i});
        Yp         = data.Y(patch_idx{:},:);
        RESULTS(i) = process_patch_object(Yp,F_dark, K, p, tau, options);
        fprintf(['Finished processing patch # ',num2str(i),' out of ',num2str(n_patches), '.\n']);
    end
else  % avoid copying the entire dataset to each worker, for in-memory data
    for i = n_patches:-1:1
        fprintf(['Running GraFT on patch # ',num2str(i),' out of ',num2str(n_patches), '.\n']);
        patch_idx         = patch_to_indices(patches{i});
        Yp                = Y(patch_idx{:},:);
        future_results(i) = parfeval(@GraFT, 2, Yp,  [], corr_kern, params);
        %future_results(i) = parfeval(@process_patch_object, 1, Yp, F_dark, K, p, tau, options);
        %value = GraFT(Yp, [], corr_kern, params);
        fprintf(['Started running GraFT on patch # ',num2str(i),' out of ',num2str(n_patches), '.\n']);
    end
    for i = 1:n_patches
        fprintf(['Fetching GraFT results from patch # ',num2str(i),' out of ',num2str(n_patches), '.\n']);
        [idx, D, S]     = fetchNext(future_results);
        RESULTS(idx).D  = D;
        RESULTS(idx).S  = S;
        RESULTS(idx).patch_idx  = i;
        fprintf(['Finished processing patch # ',num2str(i),' out of ',num2str(n_patches), '.\n']);
    end
end
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% combine results into one structure
fprintf('Combining results from different patches... \n');
cnt  = 0;
d    = prod(sizY(1:end-1));
S    = sparse(d,n_patches*K);
MASK = zeros(sizY(1:end-1));

for i = 1:n_patches
    patch_lin_idx = patch_to_linear(patches{i}, sizY);
    patch_size    = patches{i}(2:2:end) - patches{i}(1:2:end) + 1;
    for k = 1:K
        if k > size(RESULTS(i).S,2)
            break;
        end
        cnt = cnt + 1;
        S(patch_lin_idx,cnt) = reshape(RESULTS(i).S(:,:,k),[],1);
    end
    MASK(patch_lin_idx) = MASK(patch_lin_idx) + 1;
end
%%
S(:,cnt+1:end) = [];
%S              = spdiags(1./MASK(:),0,d,d)*S;
D              = cell2mat({RESULTS(:).D});
ff             = find(sum(S,1)==0);

D(:,ff) = [];                                                              % Initialize the full dictionary
S(:,ff) = [];                                                              % Initialize the full presence maps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% merge results

fprintf('Merging overlaping components...\n')

Sm = S;
Dm = D;
Km = 0;
Kn = size(S,2);

while Km < Kn
    Kn = size(Sm,2);
    [Sm,Dm,~,~] = merge_components(Y,Sm,Dm,params,MASK);
    Km = size(Sm,2);
end
fprintf(' done. \n');

if nargout > 3
    RESULTS2.S   = Sm;
    RESULTS2.D   = Dm;
    varargout{1} = RESULTS2;
    clear RESULTS2
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Final step: restimate spatial components given merged dictionary
% TODO: merge graphs from patch results instead of re-calculating graph


% fprintf('Final re-estimation of coefficient maps...')
% tic
% try
%     ckern = checkCorrKern(Y, corr_kern);                                   % Set up spatial smoothing kernel for full image
%     Sm = dictionaryRWL1SF(Y, Dm, ckern, params, full(Sm));                 % Infer coefficients given the data and dictionary
% catch ME
%     warning('Issue running funal spatial component inference!\n')
% end
% fprintf('done in %f min.\n', toc/60)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Some post-processing: reorder dictionary by magnitude

Dnorms   = sqrt(sum(Dm.^2,1));                                             % Get norms of each dictionary element
Smax     = max(Sm,[],1);                                                   % Get maximum value of each spatial map
actMeas  = Dnorms(:).*Smax(:);                                             % Total activity metric is the is the product of the above
[~,IX]   = sort(actMeas,'descend');                                        % Get the indices of the activity metrics in descending order
Dm       = Dm(:,IX);                                                       % Reorder the dictionary
Sm       = Sm(:,IX);                                                       % Reorder the spatial maps
Sm       = reshape(full(Sm),size(Y,1),size(Y,2),[]);                       % Reshape presence maps to be an image stack

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extra functions

function idx = patch_to_indices(patch)
    % helper function to build indices vector from patch start/stop indices
    idx = arrayfun(@(x,y) x:y, patch(1:2:end), patch(2:2:end), 'un', false);
end

function idx = patch_to_linear(patch, sizY)
    % helper function to build linear indices from patch start/stop indices
    slice_idx = patch_to_indices(patch);
    subs_idx = cell(1, numel(slice_idx));
    [subs_idx{:}] = ndgrid(slice_idx{:});
    subs_idx = cellfun(@(x) x(:), subs_idx, 'un', false);
    idx = sub2ind(sizY(1:end-1), subs_idx{:});
end

function result = process_patch(Y, F_dark, K, p, tau, options)
    % helper function to apply CNMF to a small patch
    sizY       = size(Y);
    options.d1 = sizY(1);
    options.d2 = sizY(2);
    if ndims(Y) == 3
        options.d3 = 1;
    else
        options.d3 = sizY(3);
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
