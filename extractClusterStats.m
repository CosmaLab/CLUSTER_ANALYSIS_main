% calculates the in-island & global statistics from cluster analysis as
% performed by 'FindClusters.m'
%
% Inputs
%   filexyn(obligatory) - full file path for the .xyn file saved by FindClusters.m
%   params - parameter list object created by FindClustersStruct.m, and 
%       used by FindClusters.m; if not input then it is obtained by loading
%       the .xyn file indicated in 'filexyn' variable
%   ClusterResults - Output from FindClusters.m; if not input then it is
%       obtained by loading the .xyn file indicated in 'filexyn' variable
%
% Outputs
%   InIslandData - structure having fields:
%       numberLocs = number of localizations per cluster for all clusters: 
%           ClusterResults(:,3)
%       clusterArea = area per cluster for all clusters:  
%           (params.original_pixel_size*ClusterResults(:,6)).^2.*pi % units = nm^2
%       nndXY = Nearest Neighbor Distance in nm for in-island clusters ONLY:  
%           ClusterResults(:,10) % units = nm
%       singleClusterIdx = binary array indicating solitary clusters:  
%           ClusterResults(:,11)==1; % isinf(InIslandData.nndXY);
%       nLocsSingleCluster = number of localizations for solitary clusters:
%           InIslandData.numberLocs( InIslandData.singleClusterIdx )
%       nLocsIslandCluster = number of localizations for in-island clusters: 
%           InIslandData.numberLocs( ~InIslandData.singleClusterIdx )
%       singleClusterArea = area per cluster for solitary clusters:
%           InIslandData.clusterArea( InIslandData.singleClusterIdx )
%       islandClustersArea = area per cluster for in-island clusters:
%           InIslandData.clusterArea( ~InIslandData.singleClusterIdx )
%       OccupiedArea = array of cluster-filled areas [total, solitary , in-island]: 
%           [sum(InIslandData.clusterArea), sum(InIslandData.singleClusterArea), sum(InIslandData.islandClustersArea)]*1e-6 % um^2
%   GlobalData - structure having fields:
%       nndXY = nearest global neighbor using X,Y position only, units = nm
%       nndXYZ = nearest global neighbor using X,Y & Z positions, units = nm
%       clusterArea = area per cluster
%       numberLocs = number of localizations per cluster
%   ClustStats - structure summarizing cluster statistics having fields:
%       nClusters = total number of clusters:
%           size(ClusterResults,1)
%       nIslands = total number of islands:
%           size(unique(ClusterResults(:,12)),1)
%       nSingleClusters = number of solitary clusters:
%           sum( InIslandData.singleClusterIdx )
%       nInIslandClusters = number of in-island clusters:
%           ClustStats.nClusters-ClustStats.nSingleClusters
%       Frac_Single = fraction of solitary clusters:
%           ClustStats.nSingleClusters/ClustStats.nClusters

function [InIslandData,GlobalData,ClustStats] = extractClusterStats(filexyn,params,ClusterResults)

%%
if ~exist('filexyn','var') || isempty(filexyn)
    [fname,fpath,~] = uigetfile([pwd '/*.xyn'],'Choose a .xyn file');
    filexyn = fullfile(fpath,fname);
elseif iscell(filexyn)
    filexyn = filexyn{1};
end
if ~exist(filexyn,'file')
    error('input file does not exist')
end
% binary var to optionally load XYN file
load_xynInfo = false(3,1);

reqfields = {'original_pixel_size','minimum_molecules_per_cluster'};
if ~exist('params','var') || isempty(params) 
    load_xynInfo([1 2]) = true;
elseif isstruct(params) 
    % ensure 'params' has the necesasry fields
    for r = 1:length(reqfields)
        if ~isfield(params,reqfields{r})
            disp('not found')
            params.(reqfields{r}) = input(['Required input: ' reqfields{r} ' = '] );
        end
    end
end

if ~exist('ClusterResults','var') || isempty(ClusterResults) || ~isnumeric(ClusterResults) ...
        || size(ClusterResults,2)~=12
    load_xynInfo([1 3]) = true;
end

if load_xynInfo(1)
    % load the xyn file
    xynInfo = XYN(filexyn);
    StrOut = 'Loaded XYN info';
    if load_xynInfo(2)
        params = xynInfo.params;
        StrOut = [StrOut ', & parameter info'];
    end
    if load_xynInfo(3)
        ClusterResults = xynInfo.data;
        StrOut = [StrOut ', & Cluster Results'];
    end
    disp(StrOut)
end
        
%%
InIslandData.numberLocs = ClusterResults(:,3);
InIslandData.clusterArea = (params.original_pixel_size*ClusterResults(:,6)).^2.*pi; % units = nm^2
InIslandData.nndXY = ClusterResults(:,10); % units = nm
InIslandData.singleClusterIdx = ClusterResults(:,11)==1; % isinf(InIslandData.nndXY);
% remove single clusters from nnd data, since they have 'Inf' value
InIslandData.nndXY = InIslandData.nndXY( ~InIslandData.singleClusterIdx );
% separate data for single vs. In-Island clusters
InIslandData.nLocsSingleCluster = InIslandData.numberLocs( InIslandData.singleClusterIdx );
InIslandData.nLocsIslandCluster = InIslandData.numberLocs( ~InIslandData.singleClusterIdx );
InIslandData.singleClusterArea = InIslandData.clusterArea( InIslandData.singleClusterIdx );
InIslandData.islandClustersArea = InIslandData.clusterArea( ~InIslandData.singleClusterIdx );
InIslandData.OccupiedArea = [sum(InIslandData.clusterArea), ... total area filled with clusters
    sum(InIslandData.singleClusterArea), ... area filled by single clusters
    sum(InIslandData.islandClustersArea)... area filled by islands of clusters
    ]*1e-6; % um^2


% Extract info on cluster & island numbers
ClustStats.nClusters = size(ClusterResults,1);
ClustStats.nIslands = size(unique(ClusterResults(:,12)),1);
ClustStats.nSingleClusters = sum( InIslandData.singleClusterIdx );
ClustStats.nInIslandClusters = ClustStats.nClusters-ClustStats.nSingleClusters;
ClustStats.Frac_Single = ClustStats.nSingleClusters/ClustStats.nClusters;


%% obtain global NND & cluster info
fileddc = [filexyn(1:end-3) 'ddc'];
performDDC = false;
if exist(fileddc,'file') % DDC previously calculated
    % load the .ddc file
    ddcInfo = DDC( fileddc );
    % ensure the DDC calculation was performed using a single .xyn file to
    % obtain global metrics
    if strcmp(ddcInfo.file1,ddcInfo.file2)
        cat = ddcInfo.data;
        % first row is expected to have 'NaN' values
        nanidx = isnan( cat{1} );
        if isempty(find(nanidx, 1)) % should not be empty
            error('check .ddc file reading, first row expected to have NaN values')
        end
        % remove NaN and check there are no remaining NaN values
        cat{1} = cat{1}(2:end,:);nanidx = isnan( cat{1} );
        if ~isempty(find(nanidx, 1)) % should be empty
            error('Error in .ddc file, found NaN values beyond first row')
        end
    else
        performDDC = true;
    end
else
    performDDC = true;
end

if performDDC
    [cat, ~] = DistanceDualColor( filexyn,  filexyn, Inf, params.original_pixel_size,...
        params.minimum_molecules_per_cluster-1, false);
    % NEED TO INCLUDE -1 FOR LOC/CLUSTER B/C DDC IS CALCULATED USING STRICT > NOT >=
    % this function will return the GLOBAL nnd, even to clusters
    % not within its same island
end

% orgainze data for output
GlobalData.nndXY = cat{1}(:,1); % units = nm
GlobalData.nndXYZ = cat{1}(:,2); % units = nm
GlobalData.clusterArea = cat{1}(:,3); % units = nm^2
GlobalData.numberLocs = cat{1}(:,4);

end % of fn