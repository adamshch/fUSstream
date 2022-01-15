function [S_mat, tau_mat] = IRWH(mov, D, corr_kern, params, tau_mat)
% Iterative reweighted L1 homotopy with spatial filtering
% using l1-homotopy toolbox from https://github.com/sasif/L1-homotopy 
% 2018 - Gal Mishne

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parse inputs

addpath(genpath('l1-homotopy/'))

if (~isfield(params,'tau'))||isempty(params.tau)
    params.tau = 2*sqrt(mean(mov,3));                                      % Default tau values to be spatially varying
else
    % Do nothing
end
if (~isfield(params,'lambda'))||isempty(params.lambda)
    lambda = 0.6;                                                          % Default lambda parameter to 0.6
else
    lambda   = params.lambda;
end
if (~isfield(params,'beta'))||isempty(params.beta)
    beta = 0.09;                                                           % Default beta parameter to 0.09
else
    beta   = params.beta;
end
if (~isfield(params,'maxiter'))||isempty(params.maxiter)
    maxiter = 0.01;                                                        % Default the maximum iteration
else
    maxiter = params.maxiter;
end
if (~isfield(params,'numreps'))||isempty(params.numreps)
    numreps = 2;
else
    numreps = params.numreps;
end
if (~isfield(params,'tolerance'))||isempty(params.tolerance)
    tolerance = 1e-8;
else
    tolerance = params.tolerance;
end
if (~isfield(params,'nonneg'))||isempty(params.nonneg)
    nonneg = false;
else
    nonneg = params.nonneg;
end
if (~isfield(params,'verbose'))||isempty(params.verbose)
    verbose = 10;
else
    verbose = params.verbose;
end
if (~isfield(params,'cardinality'))||isempty(params.cardinality)
    cardinality = 0;
else
    cardinality = params.cardinality;
end

if isempty(corr_kern)                                                      % Check if a spatial kernel is NOT provided
    corr_kern = mkCorrKern([]);                                            % If an empty array, make the kernel using the default parameters
elseif isstruct(corr_kern)
    corr_kern = mkCorrKern(corr_kern);                                     % If it is a struct then parameters were provided: make the kernel using those parameters
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get basic parameters
if size(mov,3)>1
    im_x = size(mov,1);                                                        % Get movie height
    im_y = size(mov,2);                                                        % Get movie width
    nt   = size(mov,3);                                                        % Get number of time-steps
    mov = single(reshape(mov, im_x*im_y, nt));                                 % Reshape the movie to columns for easier processing (each pixel over time is a column)
else
    im_x = params.nRows;
    im_y = params.nCols;
%     nt = size(mov,2);
end
nd   = size(D,2);
mov = mov';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set Initializations

%if nargin < 5 || isempty(S)
S = zeros([im_x, im_y, nd], 'single');                                     % Initialize event data
%end
% Initialize event data
if nargin < 5 || isempty(tau_mat)
    if numel(params.tau) == 1                                                  % Initialize trade-off parameter value/matrix
        tau_mat = lambda*single(params.tau)*ones(im_x, im_y, nd, 'single');           % If only one value is given then all coefficients at all locations have the same values
    else
        tau_mat = lambda*single(bsxfun(@times, params.tau, ones(1,1,nd)));            % Otherwise replicate tau variables along the 3rd dimension (dictioanry coefficients) and make values singles for better memeory
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% init
S_mat   = single(reshape(S, im_x*im_y, nd));                           % Reshape the time-traces to columns
S_mat = S_mat';
tau_mat = single(reshape(tau_mat, im_x*im_y, nd));                     % Reshape tradeoff values to a column
tau_mat = tau_mat';
%pk = zeros(size(S_mat));
prev_iter = struct('pk',{},'gamma',{},'iAtA',{});

parfor ll = 1:im_x*im_y
    in           = [];
    in.tau       = tau_mat(:,ll);
    in.debias    = 0;
    in.verbose   = 0;
    in.plots     = 0;
    in.record    = 0;
    in.delx_mode = 'mil';
    in.nonneg    = true;
    in.cardinality = cardinality; 
    out          = l1homotopy(D,mov(:,ll),in);
    S_mat(:,ll)         = out.x_out;
    prev_iter(ll).gamma = out.gamma;
    if ~isempty(out.gamma)
            prev_iter(ll).pk    = out.pk;
            prev_iter(ll).iAtA  = out.iAtA;
        else
            prev_iter(ll).pk    = [];
            prev_iter(ll).iAtA  = [];
    end
end

S_mat   = S_mat';
S       = reshape(S_mat, [im_x, im_y, nd]);
W_old   = tau_mat;
%% Run loop
for kk = 1:numreps
    
    % Update spike weights:
    %                tau
    % S =  ------------------------
    %        beta + |S| + |W*S|
    
    if numel(params.tau) == 1
        if size(corr_kern,1) == im_x*im_y
            S  = double(reshape(S, im_x*im_y, nd));
            W = lambda*params.tau./(beta + S + corr_kern*S);
        else
            CF = zeros(size(S));
            for i =1:size(S,3)
                CF(:,:,i) = (conv2(S(:,:,i)  ,corr_kern,'same'));
            end
            W = lambda*params.tau./(beta + S + CF);
        end
    else
        W = lambda*bsxfun(@times, params.tau, ones(1,1,nd))./(beta + S + convn(S,corr_kern,'same'));
    end
    verbPrint(verbose, 2, sprintf('Iteration %d weight update finished.',kk)) 
    
    % Reshape arrays to fascilitate parallel processing
    S_mat   = single(reshape(S, im_x*im_y, nd));                           % Reshape the time-traces to columns
    S_mat = S_mat';
    W = single(reshape(W, im_x*im_y, nd));                     % Reshape tradeoff values to a column
    W = W';
    
    % [~] = parfor_progress(im_x*im_y);
    verbPrint(verbose, 2, sprintf('Starting iteration %d inference.',kk)) 
    %tic
    parfor ll = 1:im_x*im_y
        in = [];
        in.delx_mode = 'mil';
        in.nonneg    = true;
        in.debias    = 0;
        in.verbose   = 0;
        in.plots     = 0;
        in.record    = 0;
        in.cardinality = cardinality; 

        y         = mov(:,ll);
        xh_old    = S_mat(:,ll);
        Dtr       = D'*(D*xh_old - y);
        u         = -W(:,ll).*sign(xh_old) - Dtr;
        pk_old    = Dtr + u;
        
        in.xh_old = xh_old;
        in.pk_old = pk_old;
        in.u      = u;
        in.W_old  = W_old(:,ll);
        in.W      = W(:,ll);
        in.gamma  = prev_iter(ll).gamma;
        in.iAtA   = prev_iter(ll).iAtA;
        
        out = l1homotopy(D,y,in);
        
        S_mat(:,ll)         = out.x_out;
        prev_iter(ll).gamma = out.gamma;
        if ~isempty(out.gamma)
            prev_iter(ll).pk    = out.pk;
            prev_iter(ll).iAtA  = out.iAtA;
        else
            prev_iter(ll).pk    = [];
            prev_iter(ll).iAtA  = [];
        end
    end
    %toc
    
    S_mat   = S_mat';
    S       = reshape(S_mat, [im_x, im_y, nd]);                            % Reshape back to a 3D array
    tau_mat = reshape(W', [im_x, im_y, nd]);                               % Reshape back to a matrix
    verbPrint(verbose, 2, sprintf('Iteration %d inference finished.',kk)) 
    W_old   = W;
end
