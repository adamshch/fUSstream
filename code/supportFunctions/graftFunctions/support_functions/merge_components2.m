function [A,C,nr,merged_ROIs] = merge_components2(A,C,options,normalizeSpatial,merged_ROIs)

% merging of spatially overlapping components that have highly correlated tmeporal activity
% The correlation threshold for merging overlapping components is user specified in P.merge_thr (default value 0.85)
% Inputs:
% A:            matrix of spatial components
% C:            matrix of temporal components
% options:      struct for algorithm parameters

% Outputs:
% A:            matrix of new spatial components
% C:            matrix of new temporal components
% nr:           new number of components
% merged_ROIs:  list of old components that were merged

% Adapted from: CNMF code by: 
%           Eftychios A. Pnevmatikakis, Simons Foundation, 2015
% this version: 
%           Gal Mishne (UCSD) & Adam Charles (Johns Hopkins), 2019


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing and initializations
% 
% aRatio = getCrossDictWeights(C,A);
% A      = A*diag(aRatio);

C = C';                                                                    % Transpose the time-courses for easier processing
if nargin < 3; options = []; end                                           % Make sure there is an options variable
[thr, mx]  = parseMergeParameters(options);                                % parse out the merging parameters (max & threshold)
[d, T, nr] = getProblemDims(A, C);                                         % Get dimensions of the problem
if nr == 0                                                                 % If there are no components, nothing to merge and so return an empty matrix
    merged_ROIs = [];                                                      % No components, so nothing to merge: return an empty array
    return                                                                 % Nothing else to do here
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 

if nargin < 5
    [A_corr, FF2]     = getSpatialOverlapGraph(A, nr);                     % Find graph of overlapping spatial components (positive inner products)
    C_corr            = getTemporalCorrelations(C, A_corr, nr);            % Calculate the temporal correlations betwwen overlapping spatial componnets
    MC                = getMergeIndices(C_corr, FF2, thr);                 % Find if any componentes are strongly connected and to who
    cor               = sumMergeComponentCorrelations(MC, C_corr);         % Add up all the correlations for the graph of indices to potentially merge
    [nm, merged_ROIs] = organizeIndicesToMerge(cor, MC, mx);               % Create an array where each element is a set of components to merge into a single component
else                                                                       % If merged_ROIs is provided, use those (allows for custom merging criteria)
    nm = length(merged_ROIs);                                              %  - In this case only the number of merges is needed
end

A_merged = zeros(d,nm);                                                    % Initialize the merged spatial profiles
C_merged = zeros(nm,T);                                                    % Initialize the merged temporal profiles

for i = 1:nm
    MASK = createMergingMask(A, merged_ROIs, i);                           % Create a mask over the spatial area of the profiles to merge
    nC   = calculatePatchNormalizations(C, merged_ROIs, i);                % 
    [A_merged(:,i), C_merged(i,:)] = calculateMergedROIs(C, A, MASK, ...
                                    nC, merged_ROIs, i, normalizeSpatial); % Merge all ther ROIs
                                
end

[A, C, nr] = replaceMergedComponents(A, C, merged_ROIs, nr, nm, ...
                                                      A_merged, C_merged); % Replace the merged components with the new combination

C = C';                                                                    % Transpose the time-courses for back to how they were
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extra functions 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 

function [thr, mx] = parseMergeParameters(options)

if ~isfield(options, 'merge_thr') || isempty(options.merge_thr)            % merging threshold 
    thr = 0.85; 
else
    thr = options.merge_thr;
end     
if ~isfield(options,'max_merg')                                            % maximum merging operations
    mx = Inf; 
else 
    mx = options.max_merg; 
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

function [d, T, nr] = getProblemDims(A, C)

d  = size(A,1);                                                            % Get number of pixels
T  = size(C,2);                                                            % Get number of time steps
nr = size(A,2);                                                            % Get the total number of components

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

function C_corr = getTemporalCorrelations(C, A_corr, nr)
%C_corr = corr(full(C(1:nr,:)'));
C_corr = zeros(nr);
for i = 1:nr
    overlap_indeces = find(A_corr(i,:));
    if ~isempty(overlap_indeces)
        corr_values = corr(C(i,:)',C(overlap_indeces,:)');
        C_corr(i,overlap_indeces) = corr_values;
        C_corr(overlap_indeces,i) = corr_values;
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

function [A_corr, FF2] = getSpatialOverlapGraph(A, nr)

A_corr              = triu(A'*A);                                          % Calculate all inner products between spatial components
A_corr(1:nr+1:nr^2) = 0;                                                   % Set diagonal to zero
FF2                 = A_corr > 0;                                          % Find graph of overlapping spatial components (positive inner products)
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

function MC = getMergeIndices(C_corr, FF2, thr)
FF1    = triu(C_corr)>= thr;                                               % find graph of strongly correlated temporal components
FF3    = and(FF1,FF2);                                                     % intersect the two graphs
[l, c] = graph_connected_comp(sparse(FF3+FF3'));                           % extract connected components
MC     = [];                                                               % Start with no components to merge

for i = 1:c                                                                % Loop through all connections
    if length(find(l==i))>1                                                % Find if any componentes are strongly connected and to who
%             MC = [MC,(l==i)'];
        MC = cat(2,MC,(l==i)');                                            % If there are strong connections, add that merge to the list of components to merge
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

function cor = sumMergeComponentCorrelations(MC, C_corr)

cor = zeros(size(MC,2),1);
for i = 1:length(cor)
    fm = find(MC(:,i));
    for j1 = 1:length(fm)
        for j2 = j1+1:length(fm)
            cor(i) = cor(i) + C_corr(fm(j1),fm(j2));
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

function [nm, merged_ROIs] = organizeIndicesToMerge(cor, MC, mx)

[~,ind]     = sort(cor,'descend');
nm          = min(length(ind),mx);                                         % number of merging operations
merged_ROIs = cell(nm,1);

for i = 1:nm
    merged_ROIs{i} = find(MC(:,ind(i)));
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

function MASK = createMergingMask(A, merged_ROIs, i)

MASK             = sum(A(:,merged_ROIs{i})>0,2);
MASK(MASK(:)==0) = 1;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

function nC = calculatePatchNormalizations(C, merged_ROIs, i)

nC = sqrt(sum(C(merged_ROIs{i},:).^2,2));  % Should this be the product of norm C & norm A?
nC = nC / max(nC);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

function [aa, cc] = calculateMergedROIs(C, A, MASK, nC, merged_ROIs, i, normalizeSpatial)

% Trying to find aa & cc such that aa*cc ~ A_subset*C_subset via
%    argmin_{aa,cc} ||aa*cc - A_subset*C_subset||_F^2
% Alternating approach: pseudoinverse over aa, cc
% 

A_subset = A(:,merged_ROIs{i});
C_subset = C(merged_ROIs{i},:);

A_subset = bsxfun(@times,bsxfun(@rdivide,A_subset,MASK),nC');
aa       = sum(A_subset,2);             % sum(A_subset*spdiags(nC,0,length(nC),length(nC)),2);
for iter = 1:10
    cc = (aa'*A_subset)*C_subset/sum(aa.^2); % sum(aa)? 
    aa = A_subset*(C_subset*cc')/norm(cc)^2;
end

if normalizeSpatial
    na = max(aa(:));
else
    na = 1./norm(cc(:));
end
aa = aa/na;
cc = na*cc';

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

function [A, C, nr] = replaceMergedComponents(A, C, merged_ROIs, nr, nm, A_merged, C_merged)

neur_id = unique(cell2mat(merged_ROIs));

A = [A(:,1:nr),A_merged,A(:,nr+1:end)];
C = [C(1:nr,:);C_merged;C(nr+1:end,:)];
A(:,neur_id) = [];
C(neur_id,:) = [];
nr = nr - length(neur_id) + nm;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% RELIC CODE

%S_merged = zeros(nm,T);

% MASK = reshape(MASK,[],1);
% 
% for i = 1:nm
%     nC       = sqrt(sum(C(merged_ROIs{i},:).^2,2));
%     A_subset = bsxfun(@rdivide,A(:,merged_ROIs{i}),MASK);
%     aa       = sum(A_subset*spdiags(nC,0,length(nC),length(nC)),2);
%     for iter = 1:10 
%         cc = (aa'*A_subset)*C(merged_ROIs{i},:)/sum(aa.^2);
%         aa = A_subset*(C(merged_ROIs{i},:)*cc')/norm(cc)^2;
%     end
%     na = sqrt(sum(aa.^2)/max(sum(A(:,merged_ROIs{i}).^2)));
%     aa = aa/na;
%     cc = na*cc';
%     ss = cc;
%     A_merged(:,i) = aa;    
%     C_merged(i,:) = cc;
%     S_merged(i,:) = ss;
% end

