function mov = applyMask(mov, mask, pads, p)

% mov = applyMask(mov, mask, pads, p)
% 
% This function applies a mask to a movie. 
%
% Inputs:
%  - mov    - 3D movie array
%  - mask   - 2D image with the mask as saved in the ultrasound object
%  - pads   - the padding that is needed to fit the mask to the movie. This
%               is necessary if motion correction has been applied
%  - p      - struct of parameters for fUSi processing. Important here that
%               p includes p.useMask (T/F) and p.motionCorrect (T/F)
%
% 2020 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Apply the mask if needed                                                

if p.useMask
    if p.motionCorrect                                                     % If motion corrected, the mask needs to be padded
        mask = padarray(mask, pads([1,3]), 0, 'pre');                      % Pad front of the mask to get the right motion correction padding
        mask = padarray(mask, pads([2,4]), 0, 'post');                     % Pad back of the mask to get the right motion correction padding
    else
        % do nothing
    end
    mov = bsxfun(@times, mov, mask);                                       % Apply the mask
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
