function uom = findMotionCorrectionErrorAll(uom, varargin)

% function uom = findMotionCorrectionErrorAll(uom, varargin)
% 
% Function that seeks out and quantifies residual sub-pixel motion errors
% in functional ultrasound imaging. 
% 
% 2021 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

% p = inputParser;                                                           % Set up an object to parse all of the various inputs
% p.addParameter('mcType'      , 'batchMedian');                             % Option of what motion correction type to use ('median', 'batchMedian'). This list will expand
% p.addParameter('frameSel'    , []);                                        % Option to select which frames to motion correct
% p.addParameter('refImage'    , []);                                        % Option to input a user-defined reference image (useful for multi-dataset computations)
% p.addParameter('globalOrient', false);                                      % Option to globally orient the data
% p.addParameter('medBaseLine' , 'firstmedian');                             % Select the reference image
% p.addParameter('batchSz'     , 5);                                         % Select batch sizes to reduce computation time
% p.addParameter('useDenoised' , false);                                     % Option to motion correct based on the denoised data
% p.addParameter('searchBlock' , 5);                                         % Set the limit (in pixels) that the motion correction searches for shifts over
% p.addParameter('useMask'     , true);                                      % Option to use the mask in order to compute shifts only based on brain motion (ignore out-of-brain motion)                           
% parse(p,varargin{:});
% p = p.Results;
% 

for ll = 1:size(uom.uo)
    uom.uo{ll}.findMotionCorrectionError(varargin{:});
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
