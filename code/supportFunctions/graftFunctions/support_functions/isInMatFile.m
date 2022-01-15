function bool = isInMatFile(fileName, varName)

% isInMatFile(fileName, varName)
%
% Chech if a variable called "varName" is in the file "fileName"
%
% 2018 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input checking

if ~ischar(fileName)
    error('The file name must be a string!')
end

if ~isequal(fileName(end-3:end),'.mat')
    error('The file name must p[oint to a mat file!')
end

if ischar(varName)
    vN{1} = varName;
elseif iscell(varName)
    vN = varName;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check if the variable name exists

variableInfo = who('-file', fileName);
bool         = zeros(numel(vN),1);
for kk = 1:numel(bool)
    bool(kk) = ismember(vN{kk}, variableInfo); 
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%