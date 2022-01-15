function uo = correctResidualMotion(uo, varargin)

% function fuObj = findMotionCorrectionError(fuObj, varargin)
% 
% Function that seeks out and quantifies residual sub-pixel motion errors
% in functional ultrasound imaging. 
% 
% 2020 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

p = inputParser;                                                           % Set up an object to parse all of the various inputs
p.addParameter('mcType'      , 'batchMedian');                             % 
p.addParameter('frameSel'    , []);                                        % 
p.addParameter('refImage'    , []);                                        % Option to input a user-defined reference image (useful for multi-dataset computations)
p.addParameter('medBaseLine' , 'allmedian');                               % Select the reference image
p.addParameter('batchSz'     , 100);                                       % Select batch sizes to reduce computation time
p.addParameter('useDenoised' , false);                                     % Option to use the denoised movie   
p.addParameter('reCalc'      , []);                                        % Option to recalculate the motion-corrected movie
parse(p,varargin{:});
p = p.Results;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Some initializations

if p.useDenoised; uo.ensureDenoisedData(); end                             % Ensure that the denoised data exists

if isempty(p.reCalc)
    if ~uo.isMatFileBased()
    else
        if isfile([uo.meta.dataPath,'/', uo.meta.fileName,...
                                                  '_motionCorrected.mat']) % Check that the file containing the motion correction actually exists
            uo.post.movMotionCorrected = matfile([uo.meta.dataPath,'/',...
                                 uo.meta.fileName,'_motionCorrected.mat']);% If so, get matfile pointer to the saved motion corrected data
        end
        if uo.doesMotionCorrectedExist()
            oldP = uo.post.movMotionCorrected.p;
            if isequal(p,oldP)
                p.reCalc = false;
            end
        else
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Onset of motion correction (will commence if MC needs to be recomputed)
if p.reCalc
    [motDep, motLat, ~] = getMotionVectors(uo);                                % Extract the per-frame motion errors from the uo object
    padVals = computePadValuesFromErrs(motDep, motLat);

    if p.useDenoised                                                           % Check if the denoised version of the data is best to use
        if ~uo.isMatFileBased() 
            getDataBlock = @(z) uo.post.movDenoised(:,:,z);                    % Extract a movie block
        else
            tmpMov       = uo.post.movDenoised.uoDenoised;                     % Temporarily load this movie (one movie is not that big.... usually)
            getDataBlock = @(z) tmpMov(:,:,z);                                 % Set up a function that extracts a block of the movie
        end
    else
        if ~uo.isMatFileBased() 
            getDataBlock = @(z) uo.movie(:,:,z);                               % Set up a function that extracts a block of the movie
        else
            tmpMov       = uo.meta.movMatFile.(uo.meta.movVarName);            % Temporarily load this movie (one movie is not that big.... usually)
            getDataBlock = @(z) tmpMov(:,:,z);                                 % Set up a function that extracts a block of the movie
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Correct data over batches

    szCorrMov = uo.frameSize + padVals([1,3]) + padVals([2,4]);            % Add in the padding to get the frame size for the new moveie                                             
%     uo.post.movMotionCorrected = zeros(szCorrMov(1), szCorrMov(2), ...
%                                                    uo.movieLen, 'single'); % Initialize the motion corrected movie
    uoMotionCorrected = zeros(szCorrMov(1), szCorrMov(2), ...
                                                   uo.movieLen, 'single'); % Initialize the motion corrected movie
    fprintf('Correcting residuals')
    for ll = 1:uo.movieLen
        uoMotionCorrected(:,:,ll) = single(realignSingleFrame(...
           feval(getDataBlock,ll), [motLat(ll),motDep(ll)], 'padAmt', padVals)); % Correct the ll^th frame
        if (ll>2)&&(mod(ll,100) == 0);     fprintf('.'); end
        if (ll>2)&&(mod(ll,150*100) == 0); fprintf('\n'); end
    end
    fprintf('done\n')

    uo.post.movMotionCorrected = uoMotionCorrected;
    
    if ~uo.isMatFileBased()                                                % nothing left to do: it's all saved to the object struct!
    else
        uoMotionCorrected = uo.post.movMotionCorrected;
        save([uo.meta.dataPath,'/',uo.meta.fileName,'_motionCorrected.mat'],...
                                        'uoMotionCorrected', 'p', '-v7.3');% Save the motion corrected movie to disk
        uo.post.movMotionCorrected = matfile([uo.meta.dataPath,'/',...
                                 uo.meta.fileName,'_motionCorrected.mat']);% Replace the motion corrected movie with a pointer to the saved motion corrected data
    end
    uo.post.motionPad = padVals;                                          % Save the padding values
else
    fprintf('Motion correction was already run with these parameters!')
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function padVals = computePadValuesFromErrs(motDep, motLat)

padVals = zeros(1,4);                                                      % Initialize the padding to be nothing on all sides
if min(motDep)<0;    padVals(1) = ceil(-min(motDep)); end
if max(motDep)>0;    padVals(2) = ceil(max(motDep));  end
if min(motLat)<0;    padVals(3) = ceil(-min(motLat)); end
if max(motLat)>0;    padVals(4) = ceil(max(motLat));  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%