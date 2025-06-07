function [p_values, observedClusters, observedStats, permutedStats, clusterMask, fdr_p, zStats2D, scores] = ...
    clusterPermutationTest3D_GEE(D, Y, ID, adjacencyMatrix, varargin)
%
% Calculates the cluster permutation effect for ERP data for two
% conditions. either within subject effect (is G is the same size as D)
% 
% D: per subject ERPs (channel x time x subject) matrix
% G: either 1) a matrix of size(D) that is used to subtract from D to get the
%              within subject effect
%           2) and empty ([]) indicating that D is already the difference
%              matrix (D-G)
%           3) a group vector for a between subject effect
%
% Additional parameters for varargin
% - adjacencyMatrix: matrix describing which channels are connected (numChannels x numChannels)
% - numPermutations: how many permutations to perform
% - pThreshold:      threshold at which clusters are defined (def=0.01)
% - direction:       'pos' or 'neg': positive or negative effects (def =
%                    'pos')
% - numClusters:     number of clusters to extract (def = 1)

    % Initialize input parser
    p = inputParser;

    % Add required parameters
    addRequired(p, 'D');
    addRequired(p, 'Y');
    addRequired(p, 'ID');
    addRequired(p, 'adjacencyMatrix');

    % Add optional parameters with default values
    addParameter(p, 'numPermutations', 1000);
    addParameter(p, 'pThreshold', 0.01);
    addParameter(p, 'direction', 'pos');
    addParameter(p, 'numClusters', 1);
    addParameter(p, 'doImpute', true);
    addParameter(p, 'covariates', []);
    addParameter(p, 'parallel', []);

    % Parse input arguments
    parse(p, D, Y, ID, adjacencyMatrix, varargin{:});

    
    % Extract values from parser
    numPermutations = p.Results.numPermutations;
    pThreshold = p.Results.pThreshold;
    direction = p.Results.direction;
    numClusters = p.Results.numClusters;
    doImpute = p.Results.doImpute;
    Cov = p.Results.covariates;
    tmp = p.Results.parallel;
    doParallel = true;
    if ~isempty(tmp) && (strcmpi(tmp,'off') || tmp<=1)
        doParallel = false;
    end
    
    % check direction
    if ismember(lower(direction), {'pos','+'})
        dir = 1;
    elseif ismember(lower(direction), {'neg','-'})
        dir = -1;
    else
        error('direction must be ''pos'',''neg'' ')
    end
    
    % assume a single data matrix, change is necessary later
    X=D;

    % get the right data to work on and put this in X (and G for between)
    if isempty(Y)
        error('must pass Y variable')
    end
    
    if length(Y) ~= size(D,3)
            error('Size of Y does not match up to data')
    end
    if ~any(size(Y)==1)
        error('Size of Y does not match up to data')
    end
    if all(ismember(unique(Y), [0 1]))
        warning('Will use "bonomial" link function for logistic regression')
        family = 'binomial';
    else
        family = 'gaussian';
    end
    if dir<0
        % testing can only be done for 
        Y = 1-Y;
    end

    
    % get some sizes, note that numSubjects may not be the actual number of
    % subjects (subject numbers can be repeated in ID)
    X = permute(X,[3 1 2]);
    [numSubjects, numLocations, numTimePoints] = size(X);

    % check sizes
    if size(adjacencyMatrix,1) ~= numLocations
        error('adjacencyMatrix size does not match up to the data.')
    end

    if ~isempty(Cov) && size(Cov,1) ~= numSubjects
        error('Covariates size does not match up with data size (n_rows differs)')
    end

    % remove the covariates.
    if ~isempty(Cov)
        reg = regstats(Y, Cov);
        Y = reg.r;
    end
    
    % test for positive effects only! if ID's are assumed not unique, use GEE!
        
    [~, ~, zStats, pStats] = fastGEE(Y, X(:,:), ID);
    zStats = reshape(zStats, [1 size(X,2) size(X,3)]);
    pStats = reshape(pStats, [1 size(X,2) size(X,3)]);

    % get p-values and FDR based on t and df
    %pStats = pnorm(zStats);  % Two-tailed p-values from t
    zStats2D = squeeze(zStats); 
    pStats2D = squeeze(pStats); 
    fdr_p = fdr(pStats(:));
    fdr_p = reshape(fdr_p, numLocations, numTimePoints);


    % Initialize permutation distribution of max cluster stats and details
    maxClusterStats = zeros(numPermutations, numClusters);
    % permutedClusterDetails = cell(numPermutations, numClusters);

    [sortedStats, sortedInfo] = RunPermutationGEE(X, Y, ID, family, pThreshold, numLocations, numTimePoints, adjacencyMatrix);
    
    if length(sortedInfo) < numClusters
        numClusters = length(sortedInfo);
        warning('reducing the number of clusters! Less than the requested number found')
    end
    
    % save these ordered, unpermuted data.
    observedClusters = sortedInfo;
    observedStats = sortedStats;

    % score subjects on the clusters if so requested
    if nargout>7
        scores = {};
        for clust = 1:numClusters
            P = size(sortedInfo{clust}, 1);
            
            % Convert (row, col) to linear indices for a single 2D matrix
            lin_idx = sub2ind([numLocations, numTimePoints], sortedInfo{clust}(:,1), sortedInfo{clust}(:,2));  % P × 1 linear indices
        
            % Preallocate result in scores
            scores{clust} = zeros(1, numSubjects);
            
            % Loop across subjects
            for n = 1:numSubjects
                slice = squeeze(D(:,:, n));      % extract individual matrix first from the original data
                values = slice(lin_idx);         % Mask the elements
                scores{clust}(n) = mean(values); % Average them
            end
        end
    end
    

    maxClusterStats = cell(1, numPermutations+1); % Preallocate cell array

%     if doParallel
%         % START OF PARALLEL LOOP ---------------------------------------------
%         parfor perm=1:numPermutations+1
%             rnd = randperm(numSubjects);
%             Y_permuted = Y(rnd);
%             sortedStats = RunPermutationGEE(X, Y_permuted, ID, family, Cov, pThreshold, numLocations, numTimePoints, adjacencyMatrix);
%                 
%             % Sort and store top cluster statistics and details
%             if length(sortedStats)<numClusters
%                 sortedStats = [sortedStats repelem(0,numClusters-length(sortedStats))];
%             end
%             maxClusterStats{perm}(1:numClusters) = sortedStats(1:numClusters);
%         end

    % START SERIAL LOOP --------------------------------------------------
%    else
    for perm=1:numPermutations+1
        Y_permuted = permuteClustersBySize(Y, ID);
        sortedStats = RunPermutationGEE(X, Y_permuted, ID, family, pThreshold, numLocations, numTimePoints, adjacencyMatrix);
            
        % Sort and store top cluster statistics and details
        if length(sortedStats) < numClusters
            sortedStats = [sortedStats repelem(0,numClusters-length(sortedStats))];
        end
        maxClusterStats{perm}(1:numClusters) = sortedStats(1:numClusters);
    end

    % END OF SERIAL LOOP -------------------------------------------------
%    end
    

    maxClusterStatsMatrix = cell2mat(maxClusterStats');
    
    % Ensure there are observed clusters to sort and compare
    if isempty(observedStats)
        p_values = ones(numClusters, 1); % Default p-values to 1 if no clusters found
    else
        p_values = nan(1,numClusters);
        for cl = 1:numClusters
            % note that maxClusterStats has one too many! Due to the
            % nonrandmized first iteration
            p_values(cl) = sum([maxClusterStatsMatrix(:, cl); Inf] >= observedStats(cl))/length(maxClusterStatsMatrix(:, cl));
        end
    end

    permutedStats = maxClusterStats;

    % finally create a mask for the numLocations xnumTimepoints that shows
    % where / when the clusters are.
    clusterMask = zeros(numLocations, numTimePoints);
    for cl = 1:numClusters
        clusterMask(sub2ind(size(clusterMask),observedClusters{cl}(:,1),observedClusters{cl}(:,2))) = cl;
    end
end



%-------------------------------------------------------------------------
% Function to find all connected significant indices. Call recursively.
% Added a backward in time search.

function [clusterIndices, visited] = findCluster(startIndexLoc, startIndexTime, adjacencyMatrix, visited, numLocations)
    if visited(1, startIndexLoc, startIndexTime)
        return
    end

    try
        clusterIndices = [startIndexLoc, startIndexTime];
    catch
        pause
    end
    visited(1, startIndexLoc, startIndexTime) = true; 

    mask = zeros(1, numLocations);
    mask(startIndexLoc) = 1;

    % spatially adjacent
    step = double(mask * adjacencyMatrix);
    if any(step)
        ndx = find(step);
        for loc2=ndx
            if ~visited(1, loc2, startIndexTime)
                [tmp2, visited] = findCluster(loc2, startIndexTime, adjacencyMatrix, visited, numLocations);
                clusterIndices = [clusterIndices; tmp2];
            end
        end
    end  
    % temporally adjacent. first forward in time, the  back in time (if
    % that exists, so mind the boundaries)
    if startIndexTime<size(visited,3) && ~visited(1, startIndexLoc, startIndexTime+1)
        [tmp2, visited] = findCluster(startIndexLoc, startIndexTime+1, adjacencyMatrix, visited, numLocations);
        clusterIndices = [clusterIndices; tmp2];
    end
    if startIndexTime>1 && ~visited(1, startIndexLoc, startIndexTime-1)
        [tmp2, visited] = findCluster(startIndexLoc, startIndexTime-1, adjacencyMatrix, visited, numLocations);
        clusterIndices = [clusterIndices; tmp2];
    end
end




%---------------------------------------------------------------------------------------------------------------
% Function to run each permutation: get t-values, p-values and call 
% findCluster. Return the sorted info.
% note that is G is empty a repeated test is used assuming X-permuted 
% is a difference score with randomly flipped sign
    
function [sortedStats, sortedInfo] = RunPermutationGEE(X_permuted, Y, ID, family, pThreshold, numLocations, numTimePoints, adjacencyMatrix)

    % test for positive effects only! if ID's are assumed not unique, use GEE!
    [~, ~, permutedZStats, permutedPStats] = fastGEE(Y, X_permuted(:,:), ID);
    permutedZStats = reshape(permutedZStats, [1 size(X_permuted,2) size(X_permuted,3)]);
    permutedPStats = reshape(permutedPStats, [1 size(X_permuted,2) size(X_permuted,3)]);
    
    % Determine significant indices based on the p-value threshold. Only in
    % a specific direction! To test negative effects, use a different flag.
    permutedSignificantIndices = permutedPStats < pThreshold & permutedZStats>0;

    % Cluster significant permuted stats
    visited = false(1, numLocations, numTimePoints);
    visited(~permutedSignificantIndices) = true;         
    visited(permutedZStats<=0) = true;   % fake it being visited. 
    clusterInfo = {};
    clusterCnt = 0;
    for time=1:numTimePoints
        for loc=1:numLocations
            if ~visited(1, loc, time) 
                % start looking for cluster at this point. Note that
                % this point MUST be significant.
                clusterCnt = clusterCnt + 1;
                [clusterInfo{clusterCnt}, visited] = findCluster(loc, time, adjacencyMatrix, visited, numLocations);
                % get rid of duplicates: convert to table, use unique,
                % and convert back. This may be superfluous! But the
                % function does not seem to suffer much (5% in time) by
                % leaving it in.
                T = array2table(clusterInfo{clusterCnt},"variableNames",{'Key1','Key2'});
                [~, ia] = unique(T(:, {'Key1', 'Key2'}), 'rows', 'stable');
                Tcln = T(ia,:);
                clusterInfo{clusterCnt} = table2array(Tcln);
            end
        end
    end

    % get the sum statistics for the clusters
    permutedZStats2D = squeeze(permutedZStats(1,:,:));
    permutedClusters = nan(1,length(clusterInfo)); 
    for i=1:length(clusterInfo)
        permutedClusters(i) = sum(permutedZStats2D(sub2ind(size(permutedZStats2D), clusterInfo{i}(:,1), clusterInfo{i}(:,2))));
    end
    
    % sort all clusters to size and return the info
    [sortedStats, idx] = sort(permutedClusters, 'descend');
    if nargout==2
        sortedInfo = clusterInfo(idx);
    end
end

%-------------------------------------------------------------------------------------------------------------
% function to do fast GEE for continuous variables

function [Beta, SE, Z, P] = fastGEE(Y, X, ID)

% GEE with identity link and independence correlation
% One model per X column, with robust SEs

    X = X - mean(X, 1);  % subtract column-wise mean from each predictor
    Y = Y - mean(Y);     % subtract overall mean from response

    [~, p] = size(X);
    [~, ~, clusterID] = unique(ID);
    m = max(clusterID);
    
    % Vectorized OLS estimates
    XX = sum(X.^2, 1)';      % p x 1
    XY = X' * Y;             % p x 1
    Beta = XY ./ XX;         % p x 1
    
    % Residuals: R(i,j) = Y(i) - X(i,j)*Beta(j)
    R = Y - X .* Beta';      % n x p
    
    % Accumulate sandwich middle term: sum over clusters of (sum(Xij * rij))^2
    S = zeros(p, 1);
    for i = 1:m
        idx = (clusterID == i);
        Xi = X(idx, :);        % n_i x p
        Ri = R(idx, :);        % n_i x p
    
        sumXr = sum(Xi .* Ri, 1);   % 1 x p
        S = S + (sumXr').^2;        % accumulate per-predictor
    end
    
    SE = sqrt(S ./ (XX.^2));    % robust SE per predictor
    Z = Beta./SE;
    P = pnorm(Z);
end


%-----------------------------------------------------------------------------------------------
% cluster-wise randomizer

function Y_permuted = permuteClustersBySize(Y, ID)
% ID is the cluster ID. Returns swapped data in Y

    % Setup
    [~, ~, clusterID] = unique(ID);
    all_clusters = unique(clusterID);
    
    % Get cluster sizes
    cluster_size_per_cluster = accumarray(clusterID, 1, [], @sum);  % per cluster
    
    % Group clusters by size
    max_size = max(cluster_size_per_cluster);

    % initialize Y
    Y_permuted = Y;

    for s = 1:max_size
        clusters_s = all_clusters(cluster_size_per_cluster == s);
        
        if numel(clusters_s) < 2
            continue  % need at least 2 to permute
        end
    
        % Permute cluster order
        permuted = clusters_s(randperm(numel(clusters_s)));
    
        % Now swap the entire cluster values
        for k = 1:numel(clusters_s)
            from = clusters_s(k);
            to   = permuted(k);
    
            idx_from = (clusterID == from);
            idx_to   = (clusterID == to);
    
            % Swap rows of X, or Y, or both
            Y_permuted(idx_from, :) = Y(idx_to, :);
            % If you want to permute Y instead:
            % Y_perm(idx_from) = Y(idx_to);
        end
    end
end
