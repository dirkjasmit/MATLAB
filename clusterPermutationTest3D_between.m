function [p_values, observedClusters, observedStats, permutedStats, clusterMask, stats] = clusterPermutationTest3D_between(D, G, adjacencyMatrix, varargin)

% Calculates the cluster permutation effect for ERP data for two conditions.
% 
% D              : per subject ERPs (channel x time x subject) matrix,
%                  either values or difference scores between conditions.
% G              : grouping variable. Only two values allowed for now (0 and 1).
% adjacencyMatrix: matrix describing which channels are connected
%                  (numChannels x numChannels). Take a look at
%                  GetConnectMask(), a handy tool to get an adjacency mask
%                  for EEGLAB data.
%
% optional input (format "key", value )
%  "covariates"     : covariates. will be used to remove variance from the
%                     unshuffled D data before starting permutations.
% "numPermutations" : how many permutations to perform
% "pThreshold"      : threshold at which clusters are defined (def=0.01)
% "direction"       : 'pos' or 'neg': positive or negative effects (def =
%                     'pos')
% "numClusters"     : number of clusters to extract (def = 1)

    % Initialize input parser
    p = inputParser;

    % Add required parameters
    addRequired(p, 'D');
    addRequired(p, 'G');
    addRequired(p, 'adjacencyMatrix');

    % Add optional parameters with default values
    addParameter(p, 'covariates', []);
    addParameter(p, 'numPermutations', 1000);
    addParameter(p, 'pThreshold', 0.01);
    addParameter(p, 'direction', 'pos');
    addParameter(p, 'numClusters', 1);
    addParameter(p, 'doImpute', true);

    % Parse input arguments
    parse(p, D, G, adjacencyMatrix, varargin{:});

    % Extract values from parser
    numPermutations = p.Results.numPermutations;
    pThreshold = p.Results.pThreshold;
    direction = p.Results.direction;
    numClusters = p.Results.numClusters;
    covariates = p.Results.covariates;
    doImpute = p.Results.doImpute;


    % initial work on the data.
    D = permute(D,[3 1 2]);

    % get some sizes, note that numSubjects may not be the actual number of
    % subjects (subject numbers can be repeated in ID)
    [numSubjects, numLocations, numTimePoints] = size(D);  % Size of the data matrix

    if size(G,1) ~= numSubjects
        G=G';
    end
    if size(G,1) ~= numSubjects
        error('Cannot match sise of covariates to the data.')
    end
    if ~isvector(G)
        error('Y must be a vector');
    end
    if length(unique(G))>2
        error('For the moment only two values for Y are allowed.')
    end

    if ~isempty(covariates)
        if size(covariates,1) ~= numSubjects
            covariates=covariates';
        end
        if size(covariates,1) ~= numSubjects
            error('Cannot match sise of covariates to the data.')
        end
    end

    if size(adjacencyMatrix,1) ~= numLocations
        error('adjacencyMatrix size does not match up to the data.')
    end
    if ~issymmetric(adjacencyMatrix)
        error('adjacencyMatrix must be square and symmetric.')
    end

    if ~exist('direction') || isempty(direction)
        dir = 1;
    elseif ismember(lower(direction), {'pos','+'})
        dir = 1;
    elseif ismember(lower(direction), {'neg','-'})
        dir = -1;
    else
        error('direction must be ''pos'',''neg'',or [].')
    end

    
    % flatten D and D2 into X, and impute missing values. Remove empties.
    % Apply test direction here! X is really a misnomer, as it is the
    % dependent variable! G is the grouping variable (which is the
    % predictor).
    if dir>0
        X = D(:,:);
    else
        X = -D(:,:);
    end
    if any(isnan(X(:)))
        if doImpute
            warning('Imputing missing data in D.')
            try
                X = pcambtri(X);
            catch E
                warning('Could not impute the data. Is ''MDI Toolbox 4 installed''?');
            end
        end
        valid = ~isnan(sum(X,2));
        X = X(valid,:);
        D = D(valid,:,:);
        G = G(valid);
        numSubjects = sum(valid);
        if ~isempty(covariates)
            covariates = covariates(valid,:);
        end
    end

    if ~isempty(covariates)
        residuals = zeros(size(X));
        C = [ones(numSubjects,1) covariates];
        % remove covariates from the dependent data (D)
        for i = 1:size(X,2)
            % Fit linear model to remove covariates
            b = C \ X(:, i);
            % Calculate residuals with mean added back
            residuals(:, i) = X(:, i) - C * b;
        end
        X = residuals;
    end

    % get degrees of freedom for the ttest
    df = size(X,1) - 2;

    % reshape back into three dimensions
    X = reshape(X, numSubjects, numLocations, numTimePoints);

    % test for positive effects only! if ID's are unique for each line,
    % then use a simple correlation. Otherwise you need to use GEE!
    % flatten D into X
    [~,~,~,mod] = ttest2(X(G==0,:,:),X(G~=0,:,:));
    tStats = mod.tstat;  % T-statistics for each predictor
    pStats = 2 * tcdf(-abs(tStats), df);  % Two-tailed p-values
    tStats2D = squeeze(tStats(1,:,:)); 

    fdr_p = fdr(pStats(:));
    fdr_p = reshape(fdr_p,numLocations,numTimePoints);

    stats = struct;
    stats.tstats = tStats2D;
    stats.pstats = squeeze(pStats(1,:,:));
    stats.fdr = fdr_p;
    stats.df = df;

    % Function to find all connected significant indices. Use the global
    % visited variable! This is a leak but required to pass back the
    % information during the self-call
    function [clusterIndices, visited] = findCluster(startIndexLoc, startIndexTime, adjacencyMatrix, visited)
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
                    [tmp2, visited] = findCluster(loc2, startIndexTime, adjacencyMatrix, visited);
                    clusterIndices = [clusterIndices; tmp2];
                end
            end
        end  
        % temporally adjacent. first forward in time, the  back in time (if
        % that exists, so mind the boundaries)
        if startIndexTime<size(visited,3) && ~visited(1, startIndexLoc, startIndexTime+1)
            [tmp2, visited] = findCluster(startIndexLoc, startIndexTime+1, adjacencyMatrix, visited);
            clusterIndices = [clusterIndices; tmp2];
        end
        if startIndexTime>1 && ~visited(1, startIndexLoc, startIndexTime-1)
            [tmp2, visited] = findCluster(startIndexLoc, startIndexTime-1, adjacencyMatrix, visited);
            clusterIndices = [clusterIndices; tmp2];
        end
    end

    % Initialize permutation distribution of max cluster stats and details
    maxClusterStats = zeros(numPermutations, numClusters);
    % permutedClusterDetails = cell(numPermutations, numClusters);

    for perm=1:numPermutations+1
        % Permute the outcome: swap sign of X values for random subjects.
        % Swap for all channels and time points.
        rnd = randn(numSubjects,1)<.5;
        if (perm>1) % do not permute in round 1
            Y_permuted = G(rnd);
        else
            Y_permuted = G(1:length(G));
        end

        % Calculate permuted correlations and t-statistics
        [~,~,~,tmp] = ttest2(X(Y_permuted==0,:,:),X(Y_permuted~=0,:,:));
        permutedTStats = tmp.tstat;
        permutedPStats = 2* tcdf(-abs(permutedTStats), df);
        % Determine significant indices based on the p-value threshold. Only in
        % a specific direction! To test negative effects, use a different flag.
        permutedSignificantIndices = permutedPStats < pThreshold & permutedTStats>0;

        % Cluster significant permuted stats
        permutedClusters = [];
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
                    [clusterInfo{clusterCnt} visited] = findCluster(loc, time, adjacencyMatrix, visited);
                    % get rid of duplicates: convert to table, use unique,
                    % and convert back
                    T = array2table(clusterInfo{clusterCnt},"variableNames",{'Key1','Key2'});
                    [Tcln, ia] = unique(T(:, {'Key1', 'Key2'}), 'rows', 'stable');
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
        % sort all clusters to size
        [sortedStats, idx] = sort(permutedClusters, 'descend');

        % Sort and store top cluster statistics and details
        if perm==1
            sortedInfo = clusterInfo(idx);
            if length(sortedInfo) < numClusters
                numClusters = length(sortedInfo);
            end
            % set length 
            observedClusters = sortedInfo;
            observedStats = sortedStats;
        else
            if length(sortedStats)<numClusters
                sortedStats = [sortedStats repelem(0,numClusters-length(sortedStats))];
            end
            maxClusterStats(perm,1:numClusters) = sortedStats(1:numClusters);
        end
    end
    
    % Ensure there are observed clusters to sort and compare
    if isempty(observedStats)
        p_values = ones(numClusters, 1); % Default p-values to 1 if no clusters found
    else
        p_values = [];
        for idx = 1:numClusters
            % note that maxClusterStats has one too many! Due to the
            % nonrandmized first iteration
            p_values = [p_values; sum([maxClusterStats(2:end, idx);Inf] >= observedStats(idx))/(length(maxClusterStats(:, idx))-1)];
        end
    end

    permutedStats = maxClusterStats(2:end,:);

    % finally create a mask for the numLocations x numTimepoints that shows
    % where / when the clusters are.
    clusterMask = zeros(numLocations, numTimePoints);
    for idx = 1:numClusters
        clusterMask(sub2ind(size(clusterMask),observedClusters{idx}(:,1),observedClusters{idx}(:,2))) = idx;
    end

end

