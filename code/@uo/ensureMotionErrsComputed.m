function uo = ensureMotionErrsComputed(uo)

% function uo = ensureMotionErrsComputed(uo)
%
% Function to make sure that the motion correction has been run.
%
% 2020 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run the check 

if (~isfield(uo.errs,'motion'))||isempty(uo.errs.motion)                    % Check to make sure that the motion correction error check was run
    warning('No motion correction error detection output found. Running findMotionCorrectionError with default settings...\n')
    uo.findMotionCorrectionError();                                        % Run the motion correction error detection
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%