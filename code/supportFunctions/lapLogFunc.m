function [fval, grad] = lapLogFunc(z, data, X, Y)


% [fval, grad] = lapLogFunc(z, data, X, Y)
%
% Function to compute the value and gradient (if requested) for a function 
% of the form 
%
% f(x; u,v,a,b,c) =  a*exp(-|x-u|*b)*exp(-|y-v|*c)
%
% in the LOG-domain.  z = [a,u,v,b,c] are the parameters. data is the data 
% and [X,Y] are the grid of points over data (need to be the same size). 
% Note that the data must be moved to the log-domain, but otherwise this
% is a more stable fitting procedure.
% 
% fval is the returned function value and grad is the returned gradient. 
% Note that grad is not computed if not requested (for speed). 
%
% 2020 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute the necessary values
%
data = log(data);                                                          % Transform the data to the log domain to turn exponential parameters into algebraic parameters

llFit = log(z(1))-abs(X-z(2))*z(4)-abs(Y-z(3))*z(5);                       % Calculate the current fit given the parameters (in the log domain)
erVal = data - llFit;                                                      % Calculate the l2 error signal of the current parameters in the log domain (log(data) - log(model fit))
fval  = double(nansum(reshape(erVal,[],1).^2));                            % Cost function is the sum-squares of the log-domain data

if nargout > 1
    grad(1) = -2*nansum(reshape(erVal./z(1),[],1));                        % Gradient with respect to z(1): the offset (multiplicative scale in the double exponential)
    grad(2) = -2*nansum(reshape(erVal.*(sign(X-z(2))*z(4)),[],1));         % Gradient with respect to z(2):
    grad(3) = -2*nansum(reshape(erVal.*(sign(Y-z(3))*z(5)),[],1));         % Gradient with respect to z(3):
    grad(4) = 2*nansum(reshape(erVal.*abs(X-z(2)),[],1));                  % Gradient with respect to z(4):
    grad(5) = 2*nansum(reshape(erVal.*abs(Y-z(3)),[],1));                  % Gradient with respect to z(5):
    grad    = double(grad);                                                % Make sure that the gradient is a double
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
