function denoiseMulti(muo, varargin)

% denoiseMulti(muo, varargin)
%
% Function to denoise all the datasets within the multi-ultrasound object
%
% 2020 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for ll = 1:numel(muo.fuObj)                                                % Iterate over all the data objects
    fprintf('Denoising dataset %d of %d...', ll, numel(muo.fuObj))
    muo.fuObj{ll}.denoiseTracesWavelet(varargin{:});                       % Denoise the ll^th dataset
    fprintf('done.\n')
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%