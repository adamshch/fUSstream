function varargout = correlationMap(uo, varargin)

% function uo = correlationMap(uo, varargin)
% 
% Function that creates a correlation map of the fUs data across space
% 
% 2020 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

p = inputParser;                                                           % Set up an object to parse all of the various inputs
p.addParameter('numPts'        , 3);                                       % 
p.addParameter('pixSelect'     , 'vary'  );                                % Select which pixels to select: either an array or 'vary' or 'rand'
p.addParameter('numTraces'     , 3       );                                % Select hoe many traces to select
p.addParameter('motionCorrect' , false   );                                % Select if the motion corrected data should be used
p.addParameter('denoised'      , false   );                                % Choose to display denoised movie
p.addParameter('plotOpt'       , true    );                                % Choose to display denoised movie
p.addParameter('useMask'       , true    );
parse(p,varargin{:});
p = p.Results;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

frameSize   = uo.frameSize+(p.motionCorrect==1)*sum(reshape(uo.post.motionPad,2,2),1);
getTrace    = uo.makeGetTraceFunction(p);                                  % Make the anonymous function that extracts single traces
p.pixSelect = selectFusiTraces(uo,getTrace,p.numTraces,p.pixSelect,...
                                                         p.motionCorrect); % Select the pixels to plot
traceSel    = extractPixel(uo, getTrace, p.pixSelect);                     % Extract the appropriate pixel time-traces to plot

if p.useMask;    uo.ensureMask();   end                                    % Make sure that a mask is available
if    p.motionCorrect; padVals = uo.post.motionPad;
else                 ; padVals = [];
end
maskAug = applyMask(ones(frameSize), uo.mask, padVals, p);                     % Get a modified mask

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute correlation map

corrMap = zeros(frameSize(1),frameSize(2),size(traceSel,2),'single');% Pre-allocate the correlation map
fprintf('Computing correlations')
for ll = 1:size(traceSel,2)
    for kk = 1:uo.frameSize(1)
        for nn = 1:uo.frameSize(2)
            if maskAug(kk,nn) == 0
                corrMap(kk,nn,ll) = NaN;                                   % If not in the mask, don't compute
            else
                TMPcorr           = corrcoef(traceSel(:,ll),...
                                                         getTrace(kk,nn)); % Get the correlation coefficient between pixel(kk,nn) and test pixel ll
                corrMap(kk,nn,ll) = TMPcorr(1,2);                          % Store that correlation coefficient
            end
        end 
        fprintf('.')
    end
    fprintf(':\n')
end
fprintf('done.\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Optional plotting code

if p.plotOpt
    figure(2)
    if size(corrMap,3)==1 
        imagesc(corrMap,[0,1])                                             % Plot the image
        pbaspect([uo.scanArea(2),uo.scanArea(1),1]);                       % Set the aspect ratio
        axis off; colormap gray;                                           % Remove axes and set colormap
        title(sprintf('Correlation with pixel (%d, %d)', p.pixSelect(1),...
                                                           p.pixSelect(2)))% Give a title to the image
    elseif size(corrMap,3)==2
        imagesc(cat(3,corrMap,0*corrMap(:,:,1)),[0,1])                     % Plot the image
        pbaspect([uo.scanArea(2),uo.scanArea(1),1]);                       % Set the aspect ratio
        axis off;                                                          % Remove axes 
        title(sprintf('Correlation map'))                                  % Give a title to the image
    elseif size(corrMap,3)==3
        figure; hold on
        cla
        imagesc(corrMap,[0,1])                                             % Plot the image
        plot(p.pixSelect(:,2), p.pixSelect (:,1), 'rx')
        pbaspect([uo.scanArea(2),uo.scanArea(1),1]);                       % Set the aspect ratio
        axis off;                                                          % Remove axes 
        title(sprintf('Correlation map'))                                  % Give a title to the image
        hold off
    else
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output parsing

if nargout > 0;     varargout{1} = corrMap     ; end                       % If required, output the correlation map
if nargout > 1;     varargout{2} = p.pixSelect ; end                       % If required, output the correlation map
end