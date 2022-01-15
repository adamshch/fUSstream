function uo = drawROI(uo, varargin)

% function uo = drawROI(uo, varargin)
%
% Function to draw an ROI around the brain in a functional ultrasound image
%
% 2020 - Adam Charles & Ahmed El Hady 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing 

p = inputParser;                                                           % Set up an object to parse all of the various inputs
p.addParameter('baseSelect'  , 'mean' );                                   % Choose what image to draw on
p.addParameter('scaleFactor' , 0.5    );                                   % Select scaling factor for the image
p.addParameter('saveMask'    , true   );                                   % Option to save the mask for later use
p.addParameter('figNo'       , 104    );                                   % Select figure number to plot to
p.addParameter('redraw'      , false  );                                   % Choose to redraw the mask
parse(p,varargin{:});
p = p.Results;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Display the frame

if ~p.redraw
   if isfile([uo.meta.dataPath,'/', uo.meta.fileName,'_mask.mat'])
        fprintf('Mask file detected, loading mask...')
        TMP = matfile([uo.meta.dataPath,'/',uo.meta.fileName,'_mask.mat']);% Get the matfile pointer to the previously saved denoised data
        uo.mask = TMP.maskToSave;                                          % Load up the old mask
        fprintf('done.\n')
        maskLoaded = true;                                                 % Set a flag to prevent redrawing the mask
    else
        fprintf('No prior mask file detected: Please redraw mask.\n')
        maskLoaded = false;                                                % Set a flag to indicating that the mask must be redrawn
    end 
end

if p.redraw||(~maskLoaded)
    
    baseImg = getBaseImage(uo, p);                                         % Select the base image to draw the ROI on
    maxVal  = max(abs(baseImg(:)));                                        % Compute the maximum value

    figure(p.figNo)                                                        % Open a new figure
    imagesc(baseImg/(p.scaleFactor*maxVal), [0,1])                         % Display the image
    pbaspect([uo.scanArea(2),uo.scanArea(1),1]);                           % Set the right aspect ratio
    axis off; colormap gray;                                               % Remove axes and set colormap
    title(sprintf('Please draw the mask around the brain'))                % Give a title to the image

    uo.mask = roipoly;                                                     % Ask the user to draw the ROI

    if p.saveMask
        fprintf('Saving mask...')
        maskToSave = uo.mask;
        save([uo.meta.dataPath,'/',uo.meta.fileName,'_mask.mat'],...
                                                  'maskToSave', 'p', '-v7.3'); % Save the motion corrected movie to disk
       fprintf('done.\n')
    end
else
    
end



end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function baseImg = getBaseImage(uo, p)

switch p.baseSelect
    case 'median'                                                          % If the median is selected...
        if isempty(uo.stats.medImg); uo = uo.calcBasicStats(); end         %  ...Make sure the median image was computed
        baseImg = uo.stats.medImg;                                         %  ...Extract the median image
    case 'var'                                                             % If the variance image is selected... 
        if isempty(uo.stats.varImg); uo = uo.calcBasicStats(); end         %  ...Make sure the variance image was computed
        baseImg = uo.stats.varImg;                                         %  ...Extract the variance image
    case 'mean'                                                            % If the mean is selected...
        if isempty(uo.stats.medImg); uo = uo.calcBasicStats(); end         %  ...Make sure the mean image was computed
        baseImg = uo.stats.meanImg;                                        %  ...Extract the variance image
    case 'max'                                                             % If a max image is selected
        if isempty(uo.stats.maxImg); uo = uo.calcBasicStats(); end         %  ...Make sure the max image was computed
        baseImg = uo.stats.maxImg;                                         %  ...Extract the max image
    otherwise                                                              % If nothing is selected...
        p.baseSelect = 'mean';                                             %  ...Default to using the mean image
        baseImg = getBaseImage(uo, p);                                     %  ...And get that image
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%