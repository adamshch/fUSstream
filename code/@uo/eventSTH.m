function eSTH = eventSTH(uo, varargin)

%
%
%
%
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

p = inputParser;                                                           % Set up an object to parse all of the various inputs
p.addParameter('kernLen'    , 40);
p.addParameter('kernParams' , [2,2]);
p.addParameter('kernel'     , 'gauss');
p.addParameter('baseLine'   , 'median');
p.addParameter('eventVec'   , []);
p.addParameter('motionCorrect' , false   );                                % Select if the motion corrected data should be used
p.addParameter('denoised'      , false   );                                % Choose to display denoised movie
parse(p,varargin{:});
p = p.Results;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

p.kernel = checkKernel(p);
movieOut = uo.makeGetMovieFunction(p);                                     % Get the movie given the requested options (denoised/motion corrected)

if isempty(p.eventVec);    eventVec = uo.event2timeseries()==0;
else;                      eventVec = p.eventVec;
end

if max(eventVec) == 1
    eventVec = eventVec == 1;
elseif max(eventVec) > 1
    TMP           = eventVec(eventVec>0);
    eventVec      = zeros(1, uo.movieLen);
    eventVec(TMP) = 1;
    eventVec      = eventVec == 1;
    clear TMP
end

eSTH = sum(bsxfun(@times, movieOut, reshape(conv(eventVec(:)', p.kernel(:)','same'),[1,1,uo.movieLen])),3);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function kernel = checkKernel(p)

switch lower(p.kernel)
    case {'gauss','gaussian'}
        kernel = round(-p.kernLen/2):round(p.kernLen/2);
        kernel = exp(-(kernel.^2)./(2*p.kernParams(1)));
    case 'asymgauss'
        kernel1 = round(-p.kernLen/2):0;
        kernel1 = exp(-(kernel1.^2)./(2*p.kernParams(1)));
        kernel2 = 0:round(p.kernLen/2);
        kernel2 = exp(-(kernel2.^2)./(2*p.kernParams(2)));
        kernel = [kernel1(:)', kernel2(:)'];
    case 'box'
        kernel = ones(1,p.kernLen);
    otherwise
        p.kernel = 'gauss';
        kernel   = checkKernel(p);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

