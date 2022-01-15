function TFout = doesDenoiseExist(uo)

% Outputs TRUE if the denoised movie already exists and FALSE if not 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Perform check

if uo.isMatFileBased()
    try 
        TFout = isa(uo.post.movDenoised,'matlab.io.MatFile');
    catch
        TFout = false;
    end
else
    TFout = (isfield(uo.post,'movDenoised')&&(numel(uo.post.movDenoised)==numel(uo.movie)));
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%