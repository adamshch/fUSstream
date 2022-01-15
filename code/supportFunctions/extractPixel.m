function movSel = extractPixel(uo, getTrace, pixSelect)

% function movSel = extractPixel(uo, getTrace, pixSelect)
%
% Function to extract pixels from the ultrasound object uo. 
%
% Inputs are:
%  - uo         - Ultrasound object (see @uo folder for details)
%  - getTrace   - Anonymous function that extracts a single trace from the
%                   Ultrasound movie. Takes in an [x,y] coordinate of the
%                   pixel.
%  - pixSelect  - an Nx2 array of x,y coordinates of pixels to extract
%
% 2020 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extract the requested pixels

movSel = zeros(uo.movieLen, size(pixSelect,1));                            % Preallocate array of time traces

for ll = 1:size(pixSelect,1)                                               % Loop over pixels that are requested
    movSel(:,ll) = getTrace(pixSelect(ll,1),pixSelect(ll,2));              % Extract the requested time traces
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
