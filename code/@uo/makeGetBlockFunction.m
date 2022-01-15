function getBlock = makeGetBlockFunction(uo,p)

% getBlock = makeGetBlockFunction(uo,p)
%
% Note: denoised over-rides motionCorrect
%
% 2020 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing


if (~isfield(p,'motionCorrect'))||isempty(p.motionCorrect)                 % Select if the motion corrected data should be used
    p.motionCorrect = false;
end
if (~isfield(p,'denoised'))||isempty(p.denoised)                           % Choose to display denoised movie
    p.denoised = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check that the right movies are available

if p.denoised
    uo.ensureDenoisedData();                                               % Ensure that the denoised data exists
elseif p.motionCorrect
    uo.ensureMotionCorrection();                                           % Ensure that the motion corrected data exists
end          

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Make the get trace function

if ~uo.isMatFileBased()
    if p.denoised
        getBlock = @(z) uo.post.movDenoised(:,:,z);
    elseif p.motionCorrect
        getBlock = @(z) uo.post.movMotionCorrected(:,:,z);
    else
        getBlock = @(z) uo.movie(:,:,z);
    end
else
    if p.denoised
        tmpMov = uo.post.movDenoised.uoDenoised;                           % Load denoised movie
    elseif p.motionCorrect
        tmpMov = uo.post.movMotionCorrected.uoMotionCorrected;             % Load motion corrected movie
    else
        tmpMov = uo.meta.movMatFile.(uo.meta.movVarName);                  % Load raw movie (no motion correction)
    end
    getBlock = @(z) tmpMov(:,:,z);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%