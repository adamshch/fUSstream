function TFout = doesMotionCorrectedExist(uo)

% Outputs TRUE if motion corrected movie already exists and FALSE if not 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Perform check

if uo.isMatFileBased()
    try 
        TFout = isa(uo.post.movMotionCorrected,'matlab.io.MatFile');
    catch
        TFout = false;
    end
else
    TFout = (isfield(uo.post,'movMotionCorrected')&&(numel(uo.post.movMotionCorrected)==numel(uo.movie)));
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%