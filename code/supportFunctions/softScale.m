function data = softScale(data, lims, cent)

% data = softScale(data, lims, cent)
%
% Scales down the data in 'data' using a Sigmiodal nonlinearity that 
% preserves values around the central point 'cent' and flattens out close
% to the data limits 'lims'.
%
% 2021 - Adam Charles


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a    = max(abs(lims - cent));                                              % Calculate a sigmoid parameter based on the farthest points (lims) to be retained
data = sigmoid(data, cent, 1./a);                                          % Scale down the data using a sigmoid centered around 'cent'  using the parameter 'a'

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extra functions

function X = sigmoid(X,c,a)                                                % Function to compute a sigmoid of data

X = 1./(1+exp(-(X*a - c)));                                                % This is a sigmoid

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
