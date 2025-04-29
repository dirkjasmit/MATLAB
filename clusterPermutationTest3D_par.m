function [p_values, observedClusters, observedStats, permutedStats, clusterMask, fdr_p, tStats2D, scores] = clusterPermutationTest3D_par(D, G, adjacencyMatrix, varargin)
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
    addRequired(p, 'G');
    addRequired(p, 'adjacencyMatrix');

    % Add optional parameters with default values
    addParameter(p, 'numPermutations', 1000);
    addParameter(p, 'pThreshold', 0.01);
    addParameter(p, 'direction', 'pos');
    addParameter(p, 'numClusters', 1);
    addParameter(p, 'doImpute', true);
    addParameter(p, 'covariates', []);

    % Parse input arguments
    parse(p, D, G, adjacencyMatrix, varargin{:});

    
    % Extract values from parser
    numPermutations = p.Results.numPermutations;
    pThreshold = p.Results.pThreshold;
    direction = p.Results.direction;
    numClusters = p.Results.numClusters;
    doImpute = p.Results.doImpute;
    C = p.Results.covariates;
    
    % check direction
    if ismember(lower(direction), {'pos','+'})
        dir = 1;
    elseif ismember(lower(direction), {'neg','-'})
        dir = -1;
    else
        error('direction must be ''pos'',''neg'' ')
    end
    
    % get the right data to work on
    if isempty(G)
        cDoRepeated = true;
        X = D;
        
    elseif all(size(D)==size(G))
        cDoRepeated = true;
        if dir>0
            X = D - G;
        else
            X = G - D;
        end

    else
        cDoRepeated = false;
        if length(G) ~= size(D,3)
            error('Size of G does not match up to data')
        end
        if ~any(size(G)==1)
            error('Size of G does not match up to data G (between subjects effect) or is not the same size as the other data array (within effect)')
        end
        if ~all(unique(G)==[0 1])
            error('G as grouping variable must hold only 0,1')
        end
        if dir<0
            % reverse the direction. swap group labels 1-->0 and 0-->1
            G = 1-G;
        end
    end
    
    % get some sizes, note that numSubjects may not be the actual number of
    % subjects (subject numbers can be repeated in ID)
    X = permute(X,[3 1 2]);
    [numSubjects, numLocations, numTimePoints] = size(X);

    % check sizes
    if size(adjacencyMatrix,1) ~= numLocations
        error('adjacencyMatrix size does not match up to the data.')
    end

    if ~isempty(C) && size(C,1) ~= numSubjects
        error('Covariates size does not match up with data size (n_rows differs)')
    end


    % impute missing data is requested
    if any(isnan(X(:))) && doImpute
        warning('imputing missing data in difference scores of the data (data1-data2)')
        X = pcambtri(X);
        valid = ~isnan(sum(X,2));
        X = X(valid,:);
        numSubjects = sum(valid);
        X = reshape(X, numSubjects, numLocations, numTimePoints);
    end

    % test for positive effects only! if ID's are unique for each line,
    % then use a simple correlatioon. Otherwise you need to use GEE!
    % flatten D into X
    if isempty(C)
        if cDoRepeated
            [~,~,~,mod] = ttest(X);
            df = size(X,1)-1;
        else
            [~,~,~,mod] = ttest2(X(G==0,:,:),X(G~=0,:,:));
            df = size(X,1)-2;
        end
        tStats = mod.tstat;  % T-statistics for each predictor
    else
        % run the covariates model
        if cDoRepeated
            tStats = intercept_tvals_vectorized(X(:,:), C);
            tStats = reshape(tStats, size(X));
            df = size(X,1)-1-size(C,2);
        else
            tStats = ttest2_with_covariates(X(:,:), G, C);
            tStats = reshape(tStats, size(X));
            df = size(X,1)-2-size(C,2);
        end
    end

    % get p-values and FDR based on t and df
    pStats = 2 * tcdf(-abs(tStats), df);  % Two-tailed p-values from t
    tStats2D = squeeze(tStats(1,:,:)); 
    fdr_p = fdr(pStats(:));
    fdr_p = reshape(fdr_p, numLocations, numTimePoints);


    % Initialize permutation distribution of max cluster stats and details
    maxClusterStats = zeros(numPermutations, numClusters);
    % permutedClusterDetails = cell(numPermutations, numClusters);

    if cDoRepeated
        [sortedStats, sortedInfo] = RunPermutation(X, [], pThreshold, numLocations, numTimePoints, adjacencyMatrix);
    else
        [sortedStats, sortedInfo] = RunPermutation(X, G, pThreshold, numLocations, numTimePoints, adjacencyMatrix);
    end
    
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
        for clust = 1:length(sortedInfo)
            P = size(sortedInfo{clust}, 1);
            
            % Convert (row, col) to linear indices for a single 2D matrix
            lin_idx = sub2ind([numLocations, numTimePoints], sortedInfo{clust}(:,1), sortedInfo{clust}(:,2));  % P × 1 linear indices
        
            % Preallocate result in scores
            scores{clust} = zeros(1, numSubjects);
            
            % Loop across subjects
            for n = 1:numSubjects
                slice = squeeze(D(:,:, n));           % extract indovidual matrix first
                values = slice(lin_idx);       % Mask the elements
                scores{clust}(n) = mean(values);      % Average them
            end
        end
    end
    
    
    % START OF PARALLEL LOOP -----------------------------------------------------

    maxClusterStats = cell(1, numPermutations+1); % Preallocate cell array
    % run permutation in a parallel loop
    parfor perm=1:numPermutations+1

        % Permute the outcome: swap sign of X values for random subjects.
        % Swap for all channels and time points.
        if cDoRepeated
            rnd = logical(randi(2,1,numSubjects)-1);
            X_permuted = X;
            %subplot(1,2,1); imagesc(squeeze(nanmean(X)),[-4 4]);
            X_permuted(rnd,:,:) = -X_permuted(rnd,:,:);
            %subplot(1,2,2); imagesc(squeeze(nanmean(X_permuted)),[-4 4]); drawnow
            sortedStats = RunPermutation(X_permuted, [], pThreshold, numLocations, numTimePoints, adjacencyMatrix);
        else
            rnd = randperm(numSubjects);
            G_permuted = G(rnd);
            sortedStats = RunPermutation(X, G_permuted, pThreshold, numLocations, numTimePoints, adjacencyMatrix);
        end
            
        % Sort and store top cluster statistics and details
        if length(sortedStats)<numClusters
            sortedStats = [sortedStats repelem(0,numClusters-length(sortedStats))];
        end
        maxClusterStats{perm}(1:numClusters) = sortedStats(1:numClusters);
    end
    
    % END OF PARALLEL LOOP -----------------------------------------------------
    
    
    
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




%-------------------------------------------------------------------------
% Function to run each permutation: get t-values, p-values and call 
% findCluster. Return the sorted info.
% note that is G is empty a repeated test is used assuming X-permuted 
% is a difference score with randomly flipped sign
    
function [sortedStats, sortedInfo] = RunPermutation(X_permuted, G, pThreshold, numLocations, numTimePoints, adjacencyMatrix)

     % Calculate permuted correlations and t-statistics
    if isempty(G)
        [~,~,~,mod] = ttest(X_permuted);
        df = size(X_permuted,1)-1;
    else
        [~,~,~,mod] = ttest2(X_permuted(G==0,:,:),X_permuted(G~=0,:,:));
        df = size(X_permuted,1)-2;
    end
    permutedTStats = mod.tstat;
    permutedPStats = 2*tcdf(-abs(permutedTStats), df);
    
    % Determine significant indices based on the p-value threshold. Only in
    % a specific direction! To test negative effects, use a different flag.
    permutedSignificantIndices = permutedPStats < pThreshold & permutedTStats>0;

    % Cluster significant permuted stats
    visited = false(1, numLocations, numTimePoints);
    visited(~permutedSignificantIndices) = true;         
    visited(permutedTStats<=0) = true;   % fake it being visited. 
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
    permutedTStats2D = squeeze(permutedTStats(1,:,:));
    permutedClusters = nan(1,length(clusterInfo)); 
    for i=1:length(clusterInfo)
        permutedClusters(i) = sum(permutedTStats2D(sub2ind(size(permutedTStats2D), clusterInfo{i}(:,1), clusterInfo{i}(:,2))));
    end
    
    % sort all clusters to size and return the info
    [sortedStats, idx] = sort(permutedClusters, 'descend');
    if nargout==2
        sortedInfo = clusterInfo(idx);
    end
end


function t_values = intercept_tvals_vectorized(X, C)
    % INPUT:
    %   X: (n × m) matrix of responses (each column is a Y)
    %   C: (n × p) matrix of covariates
    % OUTPUT:
    %   t_values: (1 × m) vector of t-values for intercepts

    [n, m] = size(X);
    Z = [ones(n, 1), C];              % Design matrix (n × (p+1))
    ZtZ_inv = inv(Z' * Z);            % (p+1) × (p+1)
    H = Z * ZtZ_inv * Z';             % Hat matrix
    M = eye(n) - H;                   % Residual-maker matrix

    % Project X onto residual space
    R = M * X;                        % Residuals after regressing each Y on Z
    sigma2 = sum(R.^2, 1) / (n - size(Z, 2));  % 1 × m vector of residual variances

    % Compute betas all at once
    B = ZtZ_inv * (Z' * X);           % (p+1) × m matrix of coefficients

    % Extract intercepts
    b0 = B(1, :);                     % 1 × m vector of intercepts

    % Standard errors of intercepts
    se_b0 = sqrt(sigma2 * ZtZ_inv(1,1));  % 1 × m (broadcasted scalar * vector)

    % Compute t-values
    t_values = b0 ./ se_b0;
end

function t_values = ttest2_with_covariates(X, G, C)
    % INPUTS:
    %   X: (n × m) matrix of dependent variables
    %   G: (n × 1) binary group label (0/1)
    %   C: (n × p) matrix of covariates
    % OUTPUT:
    %   t_values: (1 × m) vector of t-statistics for G

    [n, m] = size(X);
    
    % Design matrix: intercept, G, and covariates
    Z = [ones(n,1), G, C];            % n × (p+2)
    ZtZ_inv = inv(Z' * Z);            % (p+2) × (p+2)
    H = Z * ZtZ_inv * Z';             % hat matrix
    M = eye(n) - H;                   % residual-maker matrix

    % Residuals of each model
    R = M * X;                        % n × m residuals
    sigma2 = sum(R.^2, 1) / (n - size(Z,2));  % 1 × m vector of residual variances

    % Coefficients for all X columns
    B = ZtZ_inv * (Z' * X);           % (p+2) × m

    % Extract coefficient and standard error for G (2nd row)
    b_G = B(2, :);                    % 1 × m
    se_G = sqrt(sigma2 * ZtZ_inv(2,2));  % scalar × 1 × m broadcast

    % t-values
    t_values = b_G ./ se_G;
end

