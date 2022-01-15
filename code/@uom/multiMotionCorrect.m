function uom = multiMotionCorrect(uom, varargin)

% uom = multiMotionCorrect(uom, varargin)
%
%
%
% 2020 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loop through motion correct:

for ll = 1:numel(uom.uo)                                                   % Iterate over all the data objects
    fprintf('Motion correcting dataset %d of %d...\n', ll, numel(uom.uo))
    uom.uo{ll}.correctResidualMotion(varargin{:});                         % Correct residual motion errors for llth dataset, passing through all parameters
    fprintf('done.\n')
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%