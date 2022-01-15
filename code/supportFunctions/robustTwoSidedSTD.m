function [STDs, varargout] = robustTwoSidedSTD(data)

% [STDs, STD_lims] = robustTwoSidedSTD(data)
%
% Function to compute the robust two-sided standard deviation, as based
% on the mode of the data. 
% 
% 2020 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Some initial calculations

STDs      = NaN(1,2);                                                      % Initialize standard deviations
dataMode  = halfSampleMode(data);                                          % Get the mode of the distribution
idxOver   = data>dataMode;                                                 % Isolate the data greater than the mode
idxUnder  = data<dataMode;                                                 % Isolate the data less than the mode

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Next work with data over the mode

dataAugUnder = dataMode - (vec(data(idxUnder)) - dataMode);                % Create a summetric distribition for the under data
STDs(1)      = robustSTD([dataAugUnder(:);vec(data(idxUnder))]);           % Get rubust STD for augmented data
clear dataAugunder

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% First work with data over the mode

dataAugOver = dataMode - (vec(data(idxOver)) - dataMode);                  % Create a symmetric distribution for the over data
STDs(2)     = robustSTD([dataAugOver(:);vec(data(idxOver))]);              % Get rubust STD for augmented data
clear dataAugOver

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output parsing

if nargout > 1
    varargout{1} = dataMode + STDs.*[-1,1];
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
