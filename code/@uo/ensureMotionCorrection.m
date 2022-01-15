function uo = ensureMotionCorrection(uo)

% uo = ensureDenoisedData(uo)
%
% Function to ensure that the denoised version of the data exists
%
% 2020 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run check

uo.ensureMotionErrsComputed();                                             % Make sure that the motion errors were computed

if ~uo.doesMotionCorrectedExist()
    warning('Requested motion corrected movie, but no motion corrected movie available. Running rigid motion correction with base parameters...')
    uo.correctResidualMotion();
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%