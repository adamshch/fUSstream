function getTrace = makeGetTraceFunction(uo,p)

% getTrace = makeGetTraceFunction(uo,p)
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
        getTrace = @(ix,iy) uo.post.movDenoised(ix,iy,:);
    elseif p.motionCorrect
        getTrace = @(ix,iy) uo.post.movMotionCorrected(ix,iy,:);
    else
        getTrace = @(ix,iy) uo.movie(ix,iy,:);
    end
else
    if p.denoised
        tmpMov = uo.post.movDenoised.uoDenoised;                           % Load denoised movie
    elseif p.motionCorrect
        tmpMov = uo.post.movMotionCorrected.uoMotionCorrected;             % Load motion corrected movie
    else
        tmpMov = uo.meta.movMatFile.(uo.meta.movVarName);                  % Load raw movie (no motion correction)
    end
    getTrace = @(ix,iy) tmpMov(ix,iy,:);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%