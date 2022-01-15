function fOut = realignSingleFrame(fIn, offset, varargin)

% fOut = realignSingleFrame(fIn, offset, ['transMethod', transMethod, 'padAmt', padAmt, 'padVal', padVal])
%
% Wrapper for a single image translation. fIn is the image to be 
% translated, and offset is the offset by which to translate. 
%
% Optional inputs:
%  - transMethod - method by which to translate (default 'linear'). See
%                     imtranslate.m for more options. 
%  - padAmt      - 4-vector denoting the amount of padding to include on 
%                     all sides of the data after translation (default 
%                     [0,0,0,0]). 
%  - padVal      - Value with which to pad with (default NaN)
%
% 2020 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

p = inputParser;                                                           % Set up an object to parse all of the various inputs
p.addParameter('transMethod' , 'linear');                                  % Translation interpolation method
p.addParameter('padAmt'     , [0,0,0,0]);                                  % Optional padding for the array
p.addParameter('padVal'     , NaN);                                        % Optional padding for the array
parse(p,varargin{:});
p = p.Results;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Translate the image

fOut = single(fIn);                                                        % Make sure that the data is in single form
fOut = padarray(fOut,p.padAmt([1,3]),0,'pre');                             % Pad with zeros in the front
fOut = padarray(fOut,p.padAmt([2,4]),0,'post');                            % Pad with zeros in the back
fOut = imtranslate(fOut, offset, p.transMethod,'FillValues',0);            % Actually do the translation (fill in extra values with NaNs)

newPad = computeNewPadding(p.padAmt, offset);                              % Compute the padding after the offset to correctly place NaNs where no infor is to be had
fOut   = rePadImage(fOut, newPad, p.padVal);                               % Pad the new image

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function to update the padding to correctly fit the data

function newPad = computeNewPadding(padAmt, offset)

newPad = [padAmt(1:2) + offset(1)*[1,-1], padAmt(3:4) + offset(2)*[1,-1]]; % Add/subtract the offset to get the new pads
newPad = ceil(newPad);                                                     % Make sure the padding is not fractional
newPad(newPad<0) = 1;                                                      % Padding less than 0 is not allowed!

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function to pad the data after translation

function fOut   = rePadImage(fOut, newPad, padVal)

fOut(1:newPad(1),:) = padVal;                                              % Pad on the left
fOut(:,1:newPad(3)) = padVal;                                              % Pad on the top
fOut((end-newPad(2)+1):end,:) = padVal;                                    % Pad on the right
fOut(:,(end-newPad(4)+1):end) = padVal;                                    % Pad on the bottom

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
