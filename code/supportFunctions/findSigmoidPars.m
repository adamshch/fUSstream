function [a,u] = findSigmoidPars(pt1, pt2)

% [a,u] = findSigmoidPars(pt1, pt2)
% 
% Closed form parameters for fitting a sigmoidal to two points (with no
% vertical offset/scaling). Fit is for 
%
%  f(x) = 1./(1+exp(-a*(x-u))
%
%  with f(c1) = b1 and f(c2) = b2
%  and pt1 = [c1,b1], pt2 = [c2,b2]
% 
% 2020 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input checking

if pt1(1)>pt2(1)
    tmp = pt1; pt1 = pt2; pt2 = tmp;                                       % Switch points to make them be in order
    clear tmp
end

if pt1(2) > pt2(2)                                                         % If the sigmoid has the first point higher than the second
    warning('First point above second. Fit will be for 1-sigmoid()')       % Call out a warning and flip the heights to fit 1-sigmoid()
    pt1(2) = 1 - pt1(2);
    pt2(2) = 1 - pt2(2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Y'all can do the algebra. It ends up like this

a = (log((pt2(2)-pt2(2)*pt1(2))./(pt1(2)-pt2(2)*pt1(2))))./(pt2(1)-pt1(1));% Find scale value   
u = pt1(1) + log(1./pt1(2) - 1)./a;                                        % Find offself value 


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
