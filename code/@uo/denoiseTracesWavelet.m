function denoiseTracesWavelet(uo, varargin)

% denoiseTracesWavelet(fuObj, varargin)
% 
% Function to denoise functional ultrasoud data in time using wavelet
% denoising. To avoid the missing data problem, error frames from burst
% errors are filled in using cubic spline interpolation.
% 
% 2020 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

p = inputParser;                                                           % Set up an object to parse all of the various inputs
p.addParameter('smoothLvl'       , 2         );                            % Select at what wavelet level to smooth at
p.addParameter('DenoisingMethod' , 'BlockJS' );                            % For the situation that 'prctile' is chosen, which percentile should the baseline image be calculated using
p.addParameter('wDenoiseFun'     , 'wdenoise');                            % Select which pixels to select: either an array or 'vary' or 'rand'
p.addParameter('Wavelet'         , 'sym4'    );                            % Select hoe many traces to select

parse(p,varargin{:});
p = p.Results;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Some initializations    

if ~isstruct(uo.post);     uo.post = struct; end                           % If post-processing subfield doesn't exist... make it
reCalc = true;

if ~isfield(uo.post,'movDenoised'); uo.post.movDenoised = []; end

if ~uo.isMatFileBased()
    uo.post.movDenoised = uo.movie;                                        % Initialized the denoised movie (from the object if pre-loaded)
else
    if isfile([uo.meta.dataPath,'/', uo.meta.fileName,'_denoised.mat'])
        uo.post.movDenoised = matfile([uo.meta.dataPath,'/',...
                                       uo.meta.fileName,'_denoised.mat']); % Get the matfile pointer to the previously saved denoised data
    end
    if uo.doesDenoiseExist()
        oldP = uo.post.movDenoised.p;
        if isequal(p,oldP)
            reCalc = false;
        end
    else
        uo.post.movDenoised = uo.meta.movMatFile.(uo.meta.movVarName);         % Initialized the denoised movie (from disk if not pre-loaded)
    end
end

if reCalc 

    matSize = [prod(uo.frameSize),uo.movieLen];                                % Get size of reshaped matrix (pix X frames)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Check for error frames

    fprintf('Fixing burst error frames...')
    uo.post.movDenoised = fillBurstErrors(uo, uo.post.movDenoised);   % Use cubic spline interpolation to fill in the burst error frames
    fprintf('done.\n')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Function to denoise time traces

    fprintf('Smoothing...')
    if p.smoothLvl == 0                                                        % No smoothing, nothing to do
    elseif p.smoothLvl > 0
        switch p.wDenoiseFun
            case 'wdenoise'
                uo.post.movDenoised = wdenoise(...
                        double(reshape(uo.post.movDenoised,matSize)'),...
                                           p.smoothLvl,'DenoisingMethod',...
                                        p.DenoisingMethod,'Wavelet',p.Wavelet);
                uo.post.movDenoised = single(reshape(uo.post.movDenoised', ...
                    [uo.frameSize(1), uo.frameSize(2), uo.movieLen]));% Reshape to a movie
            case 'cmddenoise'
                for ll = 1:prod(uo.frameSize)
                    [i1,i2] = ind2sub(uo.frameSize,ll);
                    uo.post.movDenoised(i1,i2,ll) = cmddenoise(...
                                         uo.post.movDenoised(i1,i2,ll), ...
                                                        p.Wavelet, smoothLvl);
                end
            otherwise
                warning('No valid denoising function chosen! Skipping denoising step.\n')
        end
    else
        warning('Bad value for smoothLvl: defaulting to NOT smoothing time traces')
    end
    fprintf('Finished denoising.\n')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Output parsing

    if (uo.meta.loadMov)&&(~isempty(uo.movie))  % nothing left to do: it's all saved to the object struct!
    else
        uoDenoised = uo.post.movDenoised;
        save([uo.meta.dataPath,'/',uo.meta.fileName,'_denoised.mat'],...
                                              'uoDenoised', 'p', '-v7.3'); % Save the motion corrected movie to disk
        uo.post.movDenoised = matfile([uo.meta.dataPath,'/',...
                                       uo.meta.fileName,'_denoised.mat']); % Replace the denoised movie with a pointer to the saved motion corrected data
    end
else
    % Nothing to compute!
    fprintf('Data already denoised with the requested parameters: skipping...\n')
end



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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
