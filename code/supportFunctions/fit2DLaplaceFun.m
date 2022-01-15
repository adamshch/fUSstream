function z = fit2DLaplaceFun(data, midPt, nonOff, varargin)

% function z = fit2DLaplaceFun(data, midPt, nonOff, [offMax])
%
% Function to fit a parametric 2D Laplace function to data. 
% Inputs:
%  - data    - 2D array of values that the Laplace function needs to be
%                fit to.
%  - midPt   - coordinates within the "image" given by data where the zero
%                point lies (for making sure that the learned parameters
%                are correctly centered).
%  - nonOff  - Center of the data to keep for the fit so that the Laplace 
%                fit does not overly bias towards fitting noise in the 
%                data tail (high frequency noise). 
%  - offMax  - maximum offset to enable
%
% 2020 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin > 3
    offMax = varargin{1};
else
    offMax = 10;
end

if numel(offMax) == 1
    offMax = offMax*[1,1];
end

costType = 'logDomain';                                                    % For now the cost type is set to "log-domain" since fitting in the log domain seems to be the most robust
optType  = 'fmincon';                                                      % For now fminunc is the chose optimizer (MATLAB's built-in) but this might change in later versions

if numel(nonOff)==1
    [nonOff(1), nonOff(2)] = ind2sub(size(data),nonOff);                   % if the fit radius offset (focal point for the radius in which to fit the data) is a scalar, make it a 2x1 vector 
end

data(data<0.1*max(data(:))) = nan;                                         % Set data that's too small (less than 30% of the peak) to NaNs as it's basically noise
[X,Y] = meshgrid((1:size(data,2))-midPt(2), (1:size(data,1))-midPt(1));    % Get a meshgrid of x/y values to evaluate the Laplacian fit on
data(((X+midPt(2)-nonOff(2)).^2 + (Y+midPt(1)-nonOff(1)).^2)>50) = nan;    % Set data too far away from the data to NaNs

switch lower(costType)
    case 'logdomain'
        fLap = @(z) lapLogFunc(z, data, X, Y);                             % Log-domain function
    otherwise
        fLap = @(z) lapFunc(z, data, X, Y);                                % Standard function fit
end

switch lower(optType)
    case 'fmincon'
        opts = optimoptions('fmincon','OptimalityTolerance',1e-25, ...
               'StepTolerance',1e-10, 'Algorithm','interior-point',...
               'Display','none','SpecifyObjectiveGradient',true);
        z = fmincon(fLap, [1,0,0,1,1],[],[],[],[],[0,-offMax(1),-offMax(2),0,0],...
                                              [10,offMax(1),offMax(2),Inf,Inf],[],opts); % Solve for the parameters of the Laplace function (double exponential)
    case 'fminunc'
        opts = optimoptions('fminunc','OptimalityTolerance',1e-25, ...
                   'StepTolerance',1e-10, 'Algorithm','quasi-newton','Display',...
                                 'none','SpecifyObjectiveGradient', true); % Unconstrained version - only to be used with log-domain fit!
        z = fminunc(fLap, [1,0,0,1,1], opts);
    otherwise
        error('unknown optimization method!')
end
if any(isnan(z))
    z = zeros(size(z));
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
