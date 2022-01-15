function [fval, grad] = lapFunc(z, data, X, Y)

% [fval, grad] = lapFunc(z, data, X, Y)
%
% Function to compute the value and gradient (if requested) for a function 
% of the form 
%
% f(x; u,v,a,b,c) = a*exp(-|x-u|*b)*exp(-|y-v|*c) 
%
% z = [a,u,v,b,c] are the parameters. data is the data and [X,Y] are the 
% grid of points over data (need to be the same size). 
% 
% fval is the returned function value and grad is the returned gradient. 
% Note that grad is not computed if not requested (for speed). 
%
% 2020 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lapFit = exp(-abs(X-z(2))*z(4)-abs(Y-z(3))*z(5));                          % Compute the actual function fit given the parameters
fval   = double(nansum(reshape(data - z(1)*lapFit,[],1).^2));              % Compute the l2 error of (data - fit)      

if nargout > 1
    grad(1) = -2*nansum(reshape( (data - z(1)*lapFit).*lapFit,[],1));      % Gradient with respect to z(1)
    grad(2) = -2*nansum(reshape( (data ...
                 - z(1)*lapFit).*(z(1)*lapFit.*sign(X-z(2))*z(4)),[],1));  % Gradient with respect to z(2)
    grad(3) = -2*nansum(reshape( (data ...
                 - z(1)*lapFit).*(z(1)*lapFit.*sign(Y-z(3))*z(5)),[],1));  % Gradient with respect to z(3)
    grad(4) = -2*nansum(reshape( (data ...
                 - z(1)*lapFit).*(z(1)*lapFit.*abs(X-z(2))),[],1));        % Gradient with respect to z(4)
    grad(5) = -2*nansum(reshape( (data ...
                 - z(1)*lapFit).*(z(1)*lapFit.*abs(Y-z(3))),[],1));        % Gradient with respect to z(5)
    grad    = double(grad);                                                % Make sure the gradient vector is a double
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
