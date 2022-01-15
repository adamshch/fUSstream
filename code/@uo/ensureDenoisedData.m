function uo = ensureDenoisedData(uo)

% uo = ensureDenoisedData(uo)
%
% Function to ensure that the denoised version of the data exists
%
% 2020 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run check

if ~uo.doesDenoiseExist()
    warning('Requested denoised movie for PCA, but no denoised movie available. Running wavelet denoising with basic features...')
    uo.denoiseTracesWavelet();
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%