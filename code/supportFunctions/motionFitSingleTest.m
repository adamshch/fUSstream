function frameOffsets = motionFitSingleTest(refImg, testImg, CalcType, varargin)

% frameOffsets = motionFitSingleTest(refImg, testImg, midPt, CalcType, ['searchBlock', searchBlock])
%
% Function to test a single frame (testImg) for motion, as compared to a
% reference image (refImg). CalcType (either 'cont' or 'disc') denotes if
% the offset should be calculated as either a continuous or discrete
% offset. 'searchBlock' denotes the distance (in pixels) to search for 
% offsets over. 
%
% 2020 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

p = inputParser;                                                           % Set up an object to parse all of the various inputs
p.addParameter('searchBlock' , 50);                                        %
parse(p,varargin{:});
p = p.Results;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up the problem: normalize & compute correlation

p.searchBlock = ceil(0.5*p.searchBlock*[1,1]);                             % Turn a scalar block size to a 2d vector useful for padarray

if (norm(refImg(:))==0)||(norm(testImg(:))==0)
    frameOffsets = [0,0];                                                  %  ... Save the temporary location
else
    refImg   = refImg/norm(refImg(:));                                     % Normalize the reference image
    testImg  = testImg/norm(testImg(:));                                   % Normalize the text image
    testImg  = padarray(testImg, p.searchBlock, 0, 'both');                % Padding array allows for 'valid' to be used in the convolution
    tmpXCorr = conv2(testImg, rot90(refImg,2),'valid');                    % Compute the 2-D cross correlation function: use conv2 instead of xcorr2 since 'valid' option reduces computation

    midPt = p.searchBlock+1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%

    switch CalcType
        case 'cont'                                                        % The 'cont' function looks for a non-integer shift...
            tmpLoc = find(tmpXCorr==max(tmpXCorr(:)));                     %  ... Get the temporary maximum value
            if isempty(tmpLoc)
                tmpLoc;
            end
            z = fit2DLaplaceFun(tmpXCorr, midPt, tmpLoc, p.searchBlock);   %  ... Fit a parametrized 2D Laplace function to the cross-correlation function computed above
            frameOffsets(1) = -z(3);                                       %  ... Store the depth offset
            frameOffsets(2) = -z(2);                                       %  ... Store the AP offset
        otherwise                                                          % Otherwise just look for the max value and take that as the offset
            tmpLoc = find(tmpXCorr==max(tmpXCorr(:)));                     %  ... Get the temporary maximum value
            [x,y]  = ind2sub(size(tmpXCorr),tmpLoc(1));                    %  ... Convert the location to indices
            frameOffsets = midPt - [x,y] - offsetAdd;                      %  ... Save the temporary location
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
