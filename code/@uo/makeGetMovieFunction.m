function movieOut = makeGetMovieFunction(uo,p)

% getMovie = makeGetBlockFunction(uo,p)
%
% Note: denoised over-rides motionCorrect.
% This function basically relies on the fact that movieOut will NEVER be
% modified (aside from "reshape") so that the memory requirements are not
% doubled (i.e., this is meant to be a "shared data" copy). 
%
% 2020 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

if isempty(p); p = struct(); end

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
        movieOut = uo.post.movDenoised;
    elseif p.motionCorrect
        movieOut = uo.post.movMotionCorrected;
    else
        movieOut = uo.movie;
    end
else
    if p.denoised
        movieOut = uo.post.movDenoised.uoDenoised;                         % Load denoised movie
    elseif p.motionCorrect
        movieOut = uo.post.movMotionCorrected.uoMotionCorrected;           % Load motion corrected movie
    else
        movieOut = uo.meta.movMatFile.(uo.meta.movVarName);                % Load raw movie (no motion correction)
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%