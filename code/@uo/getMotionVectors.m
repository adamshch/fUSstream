function [motDep, motLat, maxMot] = getMotionVectors(uo)

% [motDep, motLat, maxMot] = getMotionVectors(uo)
% 
% Function to extract the frame-by-frame residual motion errors given the
% batch motion residual calculations.
% 
% 2020 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Make sure the residual motion errors were computed

uo.ensureMotionErrsComputed();                                             % Check to make sure that the motion correction error check was run
motDep = uo.errs.motion.depth;                                             % Extract the depth motion offsets
motLat = uo.errs.motion.ap;                                                % Extract the lateral motion offsets

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Use the block length to get the per-frame residual motion

if numel(motDep)~=uo.movieLen                                              % If motDep is not the right size, then the batch version of motion correction must have been used...
    motDep = repelem(motDep(:),uo.errs.motion.params.batchSz);             % In that case replecate the elements to reflect the batch processing
    motDep = motDep(1:uo.movieLen);                                        % Truncate to accound fot the last block potentially being incomplete
end

if numel(motLat)~=uo.movieLen                                              % If motLat is not the right size, then the batch version of motion correction must have been used...
    motLat = repelem(motLat(:),uo.errs.motion.params.batchSz);             % In that case replecate the elements to reflect the batch processing
    motLat = motLat(1:uo.movieLen);                                        % Truncate to accound fot the last block potentially being incomplete
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Some post-processing and stat computations

motDep(isnan(motDep)) = 0;                                                 % NaNs usually indicate a zero shift
motLat(isnan(motLat)) = 0;                                                 % NaNs usually indicate a zero shift

maxMot = max(max(abs(motDep)), max(abs(motLat)));                          % Compute the max values

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%