function D = deviateDictionary(D, mag, varargin)

% Find an E such that ||D - Dnew||_2^2 = ||E||_2^2 = mag

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

if nargin > 2
    nonneg = varargin{1};
else
    nonneg = true;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 

mag0 = mag;
Dold = D;

mag = mag*ones(1,size(D,2));
while max(mag) > 1e-3*mag0
    E = randn(size(D));                                                    % Create a random purturbation
    E = E*(diag(mag./sqrt(sum(E.^2,1))));                                % Normalize the perturbation to be of size "mag"
    D = D + E;
    if nonneg
        D(D<0) = 0;
    end
    mag = mag0 - sum((Dold-D).^2,1);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%