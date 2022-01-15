function uo = fuPCA(uo, varargin)

% fu = fuPCA(fuObj)
% 
% Function to apply PCA within the functional ultrasound struct framework
% 
% 2020 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

p = inputParser;                                                           % Set up an object to parse all of the various inputs
p.addParameter('baseSelect'    , 'median');                                % Select what criteria to calculate the baseline image as
p.addParameter('basePrctile'   , 25      );                                % For the situation that 'prctile' is chosen, which percentile should the baseline image be calculated using
p.addParameter('numPCs'        , 25      );                                % Select which pixels to select: either an array or 'vary' or 'rand'
p.addParameter('numTraces'     , []      );                                % Select how many traces to select
p.addParameter('normTraces'    , false   );                                % Choose if each trace is normalized to its own max or not
p.addParameter('dispBurstErrs' , false   );                                % Select underlay of burst error locations
p.addParameter('dispMotErrs'   , false   );                                % Select underlay of motion correction errors
p.addParameter('motionCorrect' , true    );                                % Select whether to use the smoothed traces
p.addParameter('weightMotion'  , 0       );                                % Select how much to weight motion error frames
p.addParameter('PCAtype'       , 'normal');                                % Select the flavor of PCA
p.addParameter('lambda_L'      , 5e-5    );                                % Select low-rank constraint for RPCA
p.addParameter('lambda_S'      , 5e-7    );                                % Select sparsity constraint for RPCA
p.addParameter('useMask'       , true    );

parse(p,varargin{:});
p = p.Results;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get the movie to operate on

mov = uo.makeGetMovieFunction(p);                                          % Get the movie given the requested options (denoised/motion corrected)
mov(isnan(mov)) = 0;                                                       % Remove nan values at the padded edges with zeros

if p.useMask;    uo.ensureMask();   end                                    % Make sure that a mask is available
if    p.motionCorrect; padVals = uo.post.motionPad;
else                 ; padVals = [];
end
mov = applyMask(mov, uo.mask, padVals, p);                                 % Apply the pads if motion corrected

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set potential weights

w = setWeights(uo,p);                                                      % Set potential weights based on motion errors etc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find the SVD

[U, S, V] = runPCAwrapper(uo,mov,p,w);                                     % Wrapper to the function that runs the correct PCA

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up output

PCAout.spatial  = reshape(U, [size(mov,1), size(mov,2), size(U,2)]);       % Reshape and store the spatial PCs
PCAout.temporal = V.';                                                     % Transpose V so that the PCs are organized by column
PCAout.singVals = S;                                                       % Save the singular values
PCAout.params   = p;                                                       % Save the parameters used for future reference

uo.PCAout = PCAout;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function mov = fillBurstErrors(fuObj, mov)

if ~(isfield(fuObj.errs,'burst')&&(numel(fuObj.errs.burst)==fuObj.movieLen))
    warning('No burst error detection output found. Running findBurstFrames with default settings...\n')
    fuObj.findBurstFrames(); 
end
errFrames          = fuObj.errs.burst;                                     % Extract the location of burst errors
mov(:,:,errFrames) = NaN;                                                  % Set error frames to NaNs to fill in

for ll = 1:prod(fuObj.frameSize)
    [i1,i2] = ind2sub(fuObj.frameSize,ll);
    mov(i1,i2,:) = fillmissing(mov(i1,i2,:),'spline');                     % Use piecewise cubic spline interpolation to fill in the NaNs
end

end

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function to set weights for motion error frames

function w = setWeights(fuObj,p)
w = ones(1,fuObj.movieLen);                                                % Initialize weights to zero
if p.weightMotion>0
    if ~isfield(fuObj.errs,'motion')                                       % Check to make sure that the motion correction error check was run
        warning('No motion correction error detection output found. Running findMotionCorrectionError with default settings...\n')
        fuObj.findMotionCorrectionError();                                 % Run the motion correction error detection
    end
    motDep = fuObj.errs.motion.depth;                                      % Extract the depth motion offsets
    motLat = fuObj.errs.motion.ap;                                         % Extract the lateral motion offsets

    if numel(motDep)~=fuObj.movieLen                                       % If motDep is not the right size, then the batch version of motion correction must have been used...
        motDep = repelem(motDep(:),fuObj.errs.motion.params.batchSz);      % In that case replecate the elements to reflect the batch processing
        motDep = motDep(1:fuObj.movieLen);                                 % Truncate to accound fot the last block potentially being incomplete
    end
    
    if numel(motLat)~=fuObj.movieLen                                       % If motLat is not the right size, then the batch version of motion correction must have been used...
        motLat = repelem(motLat(:),fuObj.errs.motion.params.batchSz);      % In that case replecate the elements to reflect the batch processing
        motLat = motLat(1:fuObj.movieLen);                                 % Truncate to accound fot the last block potentially being incomplete
    end

    motDep(isnan(motDep)) = 0;
    motLat(isnan(motLat)) = 0;
    
    absErrs     = sqrt(motLat.^2 + motDep.^2);
    highErrs    = absErrs>0.5;
%     w(highErrs) = 1./p.weightMotion;
    w           = 1./(p.weightMotion*highErrs);
    w           = w./max(w);
    w(isnan(w)) = 1;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% wrapper to run PCA or RPCA

function [U, S, V] = runPCAwrapper(fuObj,mov,p,w)

mov = bsxfun(@plus, mov,-mean(mov,3));                                      % Mean subtraction
mov = bsxfun(@times,mov,reshape(w,[1,1,fuObj.movieLen]));                  % Weight down the motion error frames

switch p.PCAtype
    case 'normal'
        verStr = version('-release');
        if str2double(verStr(1:4)) >= 2020
            [U, S, V] = matlab.internal.math.randLowRankSVD(...
                        reshape(double(mov), [], fuObj.movieLen),p.numPCs);
                    % [U,S,V] = svdsketch(A,tol)
        else
            [U, S, V] = svds(reshape(mov,[],fuObj.movieLen),p.numPCs);
        end

    case 'robust'
        [L,~,~,~] = solver_RPCA_Lagrangian( reshape(double(mov), ...
                              [], fuObj.movieLen) ,p.lambda_L,p.lambda_S); % Run robust PCA
        rnk = rank(L,1e-4);
        [U, S, V] = svds(L,rnk);
    otherwise
        warning('Unknown PCA type, running normal PCA')
        p.PCAtype = 'normal';
        [U, S, V] = runPCAwrapper(mov,p);
end    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

