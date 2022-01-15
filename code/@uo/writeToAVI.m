function writeToAVI(uo, video_name, varargin)

% writeToAVI(uo, video_name, contrast_val)
%
% Function to write the movie Fmov to an avi movie file. The inputs to this
% function are:
%   - Fmov         - 3D array of data to save (mov(:,:,kk) is the kk^th
%                    frame of the movie)
%   - video_name   - String containing the file-name to save the data as
%   - contrast_val - OPTIONAL input that indicates the contrast with which
%                    to plot each frame
%
% 2017 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

if ~strcmp(video_name(end-3:end),'.avi')
    error('File name should have a .avi extension!')                       % Make sure to write to a .avi file
end

p = inputParser;                                                           % Set up an object to parse all of the various inputs
p.addParameter('contrast_val' , 0.9);                                      % Pick a contrast value
p.addParameter('avg_val'      , 5);                                        % an average value number
p.addParameter('dff'          , false);                                    % Choose whether to plot the dff or not
p.addParameter('motionCorrect', false);                                    % Choose whether to plot the motion corrected data
p.addParameter('denoised'     , false);                                    % Choose whether to plot the denoised data
parse(p,varargin{:});
p = p.Results;

if isempty(p.avg_val);    p.avg_val = 5;       end                         % Make sure average value isn't empty
if isempty(p.dff);        p.dff     = false;   end                         % Make sure dff isn't empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Write video

writerObj = VideoWriter(video_name);                                       % Set up the video object to work with
open(writerObj)                                                            % Open movie object 

Fmov = makeGetMovieFunction(uo,p);                                         % Get movie

% if ~uo.isMatFileBased();    Fmov = uo.movie;                               % Get movie directly
% else;                       Fmov = uo.meta.movMatFile.(uo.meta.movVarName);% Get movie from MAT file
% end

if p.dff
    Fmed = median(Fmov,3);
    maxF = max(vec(bsxfun(@times,bsxfun(@plus,Fmov,-Fmed),1./Fmed)));      % Get the maximum fluorescence value
else
    maxF = max(Fmov(:));                                                   % Get the maximum fluorescence value
end

clims = [0*min(Fmov(:)), p.contrast_val*maxF];                             % Set constant color limits for the video frames

for kk = 1:(size(Fmov,3)-p.avg_val)
    if p.dff
        imagesc((mean(Fmov(:,:,kk:(kk+p.avg_val-1)),3) - Fmed)./Fmed,clims)% Make the image of the kk^th frame
    else
        imagesc(mean(Fmov(:,:,kk:(kk+p.avg_val-1)),3),clims)               % Make the image of the kk^th frame
    end
    pbaspect([uo.scanArea(2),uo.scanArea(1),1]);                           % Make sure the axis sizes reflect the true aspect ratio
    colormap gray                                                          % Colorscale should be gray                                       
    axis off                                                               % Remove axis numbering
    title(sprintf('Time: %3.2f', (kk-1)*uo.dt),'FontSize',20)              % Set the title
    set(gcf,'color',[1,1,1])
    drawnow
    writeVideo(writerObj,getframe(gcf));                                   % Write movie frame
    pause(0.0001)
end

close(writerObj);                                                          % Close movie object

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
