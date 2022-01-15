function plotWaveImages(varargin)

% plotWaveImages()
% 
% This function plots 2D wavefronts to demonstrate the path of a Dopplar
% signal.
%
% Optional parameters:
%   freqPlane - plane wave frequency (default 3)
%   freqCirc  - circular wave frequency (default 6)
%   figNo     - figure to plot to (default 101)
%   xLims     - range of x-axis (default [-1,1])
%   yLims     - range of y-axis (default [-1,1])
%   xDens     - density of points in x (default 1000)
%   yDens     - density of points in y (default 1000)
%   alpha     - Alpha value for semi-transparent plotting (default 0.75)
%   cmap      - Colormap selection (default 'r')                                   
%   power     - Power to sharpen color map wavefronts (default 1)
%
% 2020 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;                                                           % Set up an object to parse all of the various inputs
p.addParameter('freqPlane'   , 3);                                         % 
p.addParameter('freqCirc'    , 5);                                         % 
p.addParameter('figNo'       , 101);                                       %
p.addParameter('xLims'       , [-1,1]);                                    % 
p.addParameter('yLims'       , [-1,1]);                                    % 
p.addParameter('xDens'       , 1000);                                      % 
p.addParameter('yDens'       , 1000);                                      % 
p.addParameter('alpha'       , 0.75);                                      % 
p.addParameter('cmap'        , 'r');                                       % 
p.addParameter('power'       , 1);                                         % 
parse(p,varargin{:});
p = p.Results;

cmap  = chooseCmap(p.cmap);                                                % Choose the color map
nX    = p.xDens*(p.xLims(2)-p.xLims(1));
nY    = p.yDens*(p.yLims(2)-p.yLims(1));
xSamp = linspace(p.xLims(1),p.xLims(2),nX);
ySamp = linspace(p.yLims(1),p.yLims(2),nY);
[X,Y] = meshgrid(xSamp,ySamp); 

planeWav = (1+cos(2*pi*p.freqPlane*Y)).^p.power;
circWav  = (1+cos(2*pi*p.freqCirc*sqrt(X.^2+Y.^2))).^p.power;
circWav(sqrt(X.^2+Y.^2)>min(min(abs(p.xLims)), min(abs(p.yLims)))) = nan;

figure(p.figNo); cla;

subplot(1,2,1), imagesc(planeWav,'AlphaData',p.alpha*(~isnan(planeWav)))
axis image; axis off;
% colormap(flipud(fireprint))
colormap(cmap);

subplot(1,2,2), imagesc(circWav,'AlphaData',p.alpha*(~isnan(circWav)))
axis image; axis off; 
% colormap(flipud(fireprint))
colormap(cmap);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cmap = chooseCmap(cmapOpt)

switch cmapOpt
    case 'r'
        cmap = AdvancedColormap('wr');
    case 'g'
        cmap = AdvancedColormap('wg');
    case 'b'
        cmap = AdvancedColormap('wb');
    case 'k'
        cmap = flipud(1-0.5*((255:-1:0)'/255)*[1,1,1]);
    case 'fireprint'
        cmap = flipud(fireprint);
    otherwise
        cmap = AdvancedColormap(cmapOpt);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
