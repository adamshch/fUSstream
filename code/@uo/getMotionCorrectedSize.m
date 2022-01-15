function motCorrSize = getMotionCorrectedSize(uo)

% motCorrSize = getMotionCorrectedSize(uo)
% 
% Get the size of the motion corrected movie
% 
% 2020 - Adam Charles

if ~uo.isMatFileBased()
    motCorrSize = size(uo.post.movMotionCorrected);                        % Get the size of the motion corrected movie
else
    motCorrSize = size(uo.post.movMotionCorrected.uoMotionCorrected);      % Get the size of the motion corrected movie
end



end