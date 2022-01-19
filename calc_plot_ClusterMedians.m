% Function to compile the median values of cluster metrics, available in the
% .xyn files output by the 'FindClusters.m' function, for one or more
% datatypes.  Metrics are plotted in a notched bar plot using the function
% 'plotCompareData', which also measures for statistically significant
% differences between data types.  The calculated metrics are output
% directly and saved to an excel file for ease of export.
%
% Users must select either .xyn files that are an output from the
% FindClusters algorithm, or a .mat file that contains the following
% variables:
%   MLfileName 
%   mask 
%   xynData 
%   ClustStats 
%   ddcData 
% the easiest method to create such a .mat file is to save these variables 
% immediately after the FindClusters function has been called 
%
% IMPORTANT FILE STRUCTURE LIMITATION WHEN USING .xyn FILES:  
%   The .xyn file format contains a full file location of the Insight3 .bin
%   file used for cluster identification. This function will fail if the
%   localition of this .bin file is moved after .xyn file creation.
%   This issue arises because the .bin file is used in this function to 
%   calculate the area of the conventional and super-resolved images for 
%   density calcualtions.
%
% Inputs
%   type = structure where each field is a string corresponding to one data
%       type to be plotted
%   startlocn = (optional) starting folder location where for selection of 
%       either .xyn or .mat files 
%       Default = working directory
%   selectFileExtension = (optional) string being either 'xyn', '.xyn', 'mat' or
%       '.mat' and indicates the type of file to be selected for cluster
%       metric analysis
%       Default = 'xyn'
%   saveFileLoc = (optional) 2-element cell array of strings
%       saveFileLoc{1} = folder location where the .xlsx excel file is saved
%       saveFileLoc{2} = name of the file to be saved
%       Default = {startlocn, 'Compare Cluster Metrics.xlsx'}
%   barPosn = (optional) string having values of either 'top' or 'bottom' 
%       and indicate where the horizontal bars indicating a statistically 
%       significant difference between data types is plotted
%       Default = 'top'
%   alphaValue = (optional) alpha value required for statistically 
%       significant differences according to kstest2
%       Default = 0.05

function export = calc_plot_ClusterMedians(type,startlocn,selectFileExtension,saveFileLoc,barPosn,alphaValue)

if ~exist('type','var') || isempty(type) || ~isstruct(type)
    error('Function requires the "type" input structure variable')
end

if ~exist('barPosn','var') || isempty(barPosn) || ~ischar(barPosn) ||...
    ~sum(strcmpi(barPosn,{'top','bottom'}))
    barPosn = 'top';
elseif strcmpi(barPosn,'top')
    barPosn = 'top'; % ensure it's lower case
elseif strcmpi(barPosn,'bottom')
    barPosn = 'bottom'; % ensure it's lower case
end
   
if ~exist('alphaValue','var') || isempty(alphaValue)
    alphaValue = 0.05;
end

if ~exist('startlocn','var') || isempty(startlocn)
    startlocn = {pwd};
elseif ~iscell(startlocn) && ischar(startlocn) 
    if exist(startlocn,'dir')
        startlocn = {startlocn};
    else
        warning('input starting directory not found to exist')
        startlocn = {pwd};
    end
end

if ~exist('selectFileExtension','var') || isempty(selectFileExtension) || ...
        ~ischar(selectFileExtension) || ...
        ~sum( strcmp(selectFileExtension,{'xyn','.xyn','mat','.mat'}) )
    selectFileExtension = 'xyn';
elseif strcmp(selectFileExtension(1),'.')
    selectFileExtension = selectFileExtension(2:end);
end
% choose cluster-info files for data extraction
datatypes = fieldnames(type);
ntypes = size(datatypes,1);
if length(startlocn) ~= ntypes
    startlocn = repmat(startlocn(1),ntypes,1);
end
files = cell(ntypes,1);
typenames = struct2cell(type)';
maxNfiles = 0;
for t = 1:ntypes
    % First, ask to select .mat files saved during cluster analysis
    files{t} = Select1DataGroup([typenames{t} ' cluster.mat'],...
        ['*.' selectFileExtension],startlocn{t});
    if ~isempty(files{t})
        f = size(files{t}.data,1);
        if f>maxNfiles, maxNfiles = f; end
    end
end

if ~exist('saveFileLoc','var') || isempty(saveFileLoc) || ~iscell(saveFileLoc)...
        || length(saveFileLoc)~=2
    saveFileLoc{1} = startlocn{1};
    saveFileLoc{2} = 'Compare Cluster Metrics.xlsx';
else
    [~,fname,ext] = fileparts(saveFileLoc{2});
    if isempty(ext) || ~strcmpi(ext,'.xlsx')
        saveFileLoc{2} = [fname '.xlsx'];
    end
end


% for i = 1:4, typenames{i} = ['Type ' num2str(i)]; end
%% prepare the output data structure
init = nan(maxNfiles,ntypes);
standard = struct(...
    'numberLocs',init,...   Median number of localizations per cluster
    'clusterArea',init,...  Median area per cluster
    'locDensity',init,...   Localization density within the cluster
    'nnd',init);          % Median nearest neighbor distance

% %%
dataFields = ...
    {... FieldName      Description                 units or subfield description + units
    'island',           'Standard cluster metrics',{'Localizations/Cluster','Area/Cluster','Localization Density/cluster', 'In-Island Nearest Neighbor';...
                                                     'Count',               'Area [nm^2]', 'Loc Density [1/nm^2]',         'NND [nm]'};...  1
    'globalNND',        'Global Nearest Neighbor',  'NND [nm]';...	2
    'SRfillConv',       'S.R. filled area',         'Percent (%)';... 3
    'ClustFillSR',      'Cluster filled areas',     {'SR Area Clustered: All', 'SR Area Clustered: Solo', 'SR Area Clustered: In-Island';...
                                                     'Percent (%)',            'Percent (%)',             'Percent (%)'};... 4
    'Ratio_StoI',       'Single:In-Island ratio',   {'Solo:Island Cluster Area','Solo:Island Num Clusters';...
                                                     'Ratio','Ratio'};...       5
    'LocDensity',       'Localization Density',     'Density [\mum^{-2}]';... 6
    'ClustDensity',     'Cluster Density',          'Density [\mum^{-2}]';... 7
    'IslandDensity',    'Island Density',           'Density [\mum^{-2}]';... 8
    'SingleClustpct',   'Single Clusters',          'Percent (%)';...        9
    'ClustPerIsland',   'Clusters per Island',      'Clusters/Island';...   10
    'Dataset',          'Datasets',                 ''};%                    11

% %%
export = struct( ...
    dataFields{1,1},standard, ...	in-island cluster metrics
    dataFields{2,1},init, ...       global NND
    dataFields{3,1},init,...        SR area / Conv area
    dataFields{4,1},struct( ...     area cluster / SR area
        'total',init,...                  all clusters
        'single',init,...               solitary clusters
        'island',init),...              in-island clusters
    dataFields{5,1},struct( ...     solitary / in-island
        'ClustArea',init,...            area_s / area_ii
        'Nclust',init),...              nSingleClusters / nInIslandClusters
    dataFields{6,1},init,...        numOrigLocs / SR Area
    dataFields{7,1},init,...        nClusters / SR Area
    dataFields{8,1},init,...        nIslands / SR Area
    dataFields{9,1},init,...        nSingleClusters / nClusters
    dataFields{10,1},init,...        nInIslandClusters / nIslands, islands have more than 1 cluster 
    dataFields{11,1},{cell(maxNfiles,ntypes)} );

% ask user which fields to calculate for export
calcVal = false(size(dataFields,1),1);
calcVal(end) = true; % Dataset is exported by default

longestStr = max(max(cell2mat(cellfun(@size,dataFields(:,2),'UniformOutput',false))));
% %%
[selection,button] = listdlg(...
    'ListString',dataFields(1:end-1,2),... Dataset is exported by default
    'SelectionMode','multiple',...     
    'ListSize',[1 1]*300*longestStr/30,...
    'InitialValue',1,...    
    'Name','Cluster Metrics',...     
    'PromptString','Choose metrics to be analyzed');%,...    'OKString','OK',...    'CancelString','Choose All' )
if button % cancelling the selection termiunates analysis
% reset metrics to be calculated
calcVal( selection ) = true;
%%
% initialize user choice to save .mat file
saveMat = 'No'; 
% load cluster data and calculate desired metrics
for t = 1:ntypes
    if ~isempty(files{t})
        nfiles = size(files{t}.data,1);
        for f = 1:nfiles
            dataloaded = false;
            % First, check the .mat file to see if it has the needed
            % information saved
            currfilename = fullfile(files{t}.data{f,2},files{t}.data{f,1});
            [~,~,currext] = fileparts(currfilename);
            if strcmp(currext,'.mat')
                matDataTmp = load( currfilename );
                select2files = false;
                if isfield(matDataTmp,'xynData') && isfield(matDataTmp,'ddcData') ...
                        && isfield(matDataTmp,'ClustStats') && isfield(matDataTmp,'mask')
                    
                    MLfileName = matDataTmp.MLfileName;
                    if isfield(matDataTmp,'numOrigLocs')
                        numOrigLocs = matDataTmp.numOrigLocs;
                    else
                        % load localization list to get number of
                        % localizations, since it wasn't previously saved
                        LL = Insight3( MLfileName );
                        numOrigLocs = LL.numMolecules;
                    end
                    mask = matDataTmp.mask;
                    xynData = matDataTmp.xynData;
                    ClustStats = matDataTmp.ClustStats;
                    ddcData = matDataTmp.ddcData;
                    dataloaded = true;
                    
                elseif isfield(matDataTmp,'MLfileName');
                    % needed info missing, but the original localization list full
                    % file path is included.
                    MLfileName = matDataTmp.MLfileName;
                    % Now, look for its associated.xyn file, then load it
                    [fpath,fname,ext] = fileparts(MLfileName);
                    xyndir = [fpath filesep 'cluster_analysis'];
                    if exist(xyndir,'dir')
                        xynfile = recursiveFindFile(xyndir,['*' fname '*.xyn']);
                        if size(xynfile,1) > 1
                            % presume didn't find it
                            select2files = true;
                        else
                            % found unique XYN file
                            xynfile = xynfile{1};
                        end
                    else
                        select2files = true;
                    end
                    
                    if select2files
                        error(['Need to make select2files code: t = ' num2str(t) ', f = ' num2str(f)])
                        % select the original .bin file (needed for mask calc) and the
                        % associated .xyn file resulting from cluster analysis
                        
                        % load the two files and calculate the needed info for xynData
                        
                    end
                    
                    % now have the localization list & xyn file locations
                    LLxyn = XYN(xynfile);
                    % ensure that we have the right localization list
                    if strcmp( MLfileName, LLxyn.params.i3file )
                        % along with it's localization list
                        LL = Insight3( MLfileName );
                    else
                        error(['couldn''t find corresponding .xyn & .bin files,'...
                            't = ' num2str(t) ', f = ' num2str(f)])
                    end
                    
                end
                
            elseif strcmp(currext,'.xyn')
                LLxyn = XYN( currfilename );
                MLfileName = LLxyn.params.i3file;
                % ensure that the localization list exists
                if ~exist( MLfileName,'file' )
                    % try looking this folder and ones above
                    [~,MLname,MLext] = fileparts( MLfileName );
                    [xynPath,xynFile,xynExt] = fileparts( currfilename );
                    xynPath = regexp(xynPath,filesep,'split');
                    remEntry = 0;
                    while ~exist( MLfileName,'file' )
                        if remEntry==length(xynPath)
                            errStruct = ...
                                struct(...
                                'message',...
                                ['Working on xyn file: ' xynFile xynExt ...
                                ' localization list: ' LLxyn.params.i3file ...
                                ' could not be located'],...
                                'identifier',...
                                'calc_plot_ClusterMedians:xynInfo:I3_binFileNotFound');
                            error( errStruct )
                        end
                        MLfileName = fullfile(xynPath{1:end-remEntry},[MLname MLext]);
                        remEntry = remEntry+1;
                    end
                end
                LL = Insight3( MLfileName );
                
            else
                error(['Selected filetype is unrecognized: ' currext ]);
            end
            
            if ~dataloaded
                % obtain original number of localizations used for clustering
                numOrigLocs = LL.numMolecules;
                mask = calcMaskAreas( LL, LLxyn.params );
                % obtain the xyn & ddc information
                params = LLxyn.params;
                ClusterResults = LLxyn.data;
                [xynData,ddcData,ClustStats] = extractClusterStats(currfilename,params,ClusterResults);
                
                % ask user if they want to save the info to a .mat file for
                % faster future data extraction
                if t == 1 && f == 1
                    % only ask on the first loop
                    qstring = 'Save Cluster Metric data to .mat file for faster loading in future?';
                    saveMat = questdlg(qstring,...
                        'Save Data as .mat?',...
                        'Yes','No','No');
                    if isempty(saveMat), saveMat = 'No'; end
                end
                switch saveMat
                    case 'Yes'
                        %                     filexyn = create_xynFilename( params );
                        [sPath,sName,~] = fileparts( LLxyn.filename );
                        savefile = fullfile(sPath, [sName '.mat']);
                        % save the info in a .mat file
                        save(savefile,...
                            'MLfileName', 'mask', 'params', 'ClusterResults',...
                            'numOrigLocs', 'xynData', 'ddcData', 'ClustStats')
                end
            end
            
            % --- End data loading, start calculating metrics ---
            
            
            currField = dataFields{1,1};
            if calcVal(~cellfun(@isempty,strfind(dataFields(:,1),currField)))
                % Median number of localizations per cluster for in-island clusters
                export.(currField).numberLocs(f,t) = median(xynData.numberLocs);
                % Median area per cluster for in-island clusters
                export.(currField).clusterArea(f,t) = median(xynData.clusterArea);
                % Median localization density within each cluster
                export.(currField).locDensity(f,t) = median(xynData.numberLocs./xynData.clusterArea);
                % Median in-island NND
                export.(currField).nnd(f,t) = median(xynData.nndXY);
            else
                if isfield(export,currField)
                    export = rmfield(export,currField);
                end
            end
            
            currField = dataFields{2,1};
            if calcVal(~cellfun(@isempty,strfind(dataFields(:,1),currField)))
%                 % Median number of localizations per cluster globally
%                 export.(currField).numberLocs(f,t) = median(ddcData.numberLocs);
%                 % Median area per cluster for all clusters
%                 export.(currField).clusterArea(f,t) = median(ddcData.clusterArea);
%                 % Median localization density within each cluster
%                 export.(currField).locDensity(f,t) = median(ddcData.numberLocs./ddcData.clusterArea);
                % Median global NND
                export.(currField)(f,t) = median(ddcData.nndXY);
                %             export.(currField).nnd(f,t) = median(ddcData.nndXY);
            else
                if isfield(export,currField)
                    export = rmfield(export,currField);
                end
            end
            
            currField = dataFields{3,1};
            if calcVal(~cellfun(@isempty,strfind(dataFields(:,1),currField)))
                % percent of conventional area filled with localizations
                export.(currField)(f,t) = 100*mask.SRArea/mask.convArea;
            else
                if isfield(export,currField)
                    export = rmfield(export,currField);
                end
            end
            
            currField = dataFields{4,1};
            if calcVal(~cellfun(@isempty,strfind(dataFields(:,1),currField)))
                % percent of SR area filled with clusters
                export.(currField).total(f,t) =  100*xynData.OccupiedArea(1)/mask.SRArea;
                % percent of SR area filled with solitary clusters
                export.(currField).single(f,t) =  100*xynData.OccupiedArea(2)/mask.SRArea;
                % percent of SR area filled with in-island clusters
                export.(currField).island(f,t) =  100*xynData.OccupiedArea(3)/mask.SRArea;
            else
                if isfield(export,currField)
                    export = rmfield(export,currField);
                end
            end
            
            currField = dataFields{5,1};
            if calcVal(~cellfun(@isempty,strfind(dataFields(:,1),currField)))
                % Ratio of area in single vs. island clusteres
                export.(currField).ClustArea(f,t) = xynData.OccupiedArea(2)/xynData.OccupiedArea(3);
                % Ratio of numbers of solitary clusters / in-island clusters
                export.(currField).Nclust(f,t) = ClustStats.nSingleClusters/ClustStats.nInIslandClusters;
            else
                if isfield(export,currField)
                    export = rmfield(export,currField);
                end
            end
            
            currField = dataFields{6,1};
            if calcVal(~cellfun(@isempty,strfind(dataFields(:,1),currField)))
                % Average Localization density
                export.(currField)(f,t) = numOrigLocs/mask.SRArea;
            else
                if isfield(export,currField)
                    export = rmfield(export,currField);
                end
            end
            
            currField = dataFields{7,1};
            if calcVal(~cellfun(@isempty,strfind(dataFields(:,1),currField)))
                % Total cluster density
                export.(currField)(f,t) = ClustStats.nClusters/mask.SRArea;
            else
                if isfield(export,currField)
                    export = rmfield(export,currField);
                end
            end
            
            currField = dataFields{8,1};
            if calcVal(~cellfun(@isempty,strfind(dataFields(:,1),currField)))
                % Density of islands
                export.(currField)(f,t) = ClustStats.nIslands/mask.SRArea;
            else
                if isfield(export,currField)
                    export = rmfield(export,currField);
                end
            end
            
            currField = dataFields{9,1};
            if calcVal(~cellfun(@isempty,strfind(dataFields(:,1),currField)))
                % Percent single clusters / all clusters
                export.(currField)(f,t) = 100*ClustStats.Frac_Single;
            else
                if isfield(export,currField)
                    export = rmfield(export,currField);
                end
            end
            
            currField = dataFields{10,1};
            if calcVal(~cellfun(@isempty,strfind(dataFields(:,1),currField)))
                % Number of clusters per island
                export.(currField)(f,t) = ClustStats.nClusters / ClustStats.nIslands;
                % alternative if we want to exclude islands having only one cluster
%                 export.(currField)(f,t) = ClustStats.nInIslandClusters/(ClustStats.nIslands - ClustStats.nSingleClusters);
            else
                if isfield(export,currField)
                    export = rmfield(export,currField);
                end
            end
            
            % Dataset name
            [~,export.Dataset{f,t},~] = fileparts(MLfileName);
            
        end % files for-loop
    end % files if ~isempty statement
end% types for-loop

%% Organize data structure for saving into excel format
% collapse into a cell
plotData = struct2cell( export );
% remove Datasets
plotData = plotData(1:end-1);
for i = 1:length( plotData )
    if isstruct( plotData{i} );
        plotData{i} = struct2cell( plotData{i} )';
    end
end

% find the total number of metrics 
nel = cell2mat(cellfun(@size,plotData,'UniformOutput',false));
cellEntry = cell2mat(cellfun(@iscell,plotData,'UniformOutput',false));
nDatCol = sum([ nel(cellEntry,2); sum(~cellEntry)]);

% expand each data entry into a single data cell
output = cell(1, nDatCol);
m = 0;
for i = 1:length( plotData )
    
    if iscell( plotData{i} )
        % if data cell has multiple entries, save each in a separate cell
        for j = 1:size( plotData{i},2)
            m = m+1;
            output{m} = plotData{i}{j};
        end
    else
        % otherwise, there's only a single data entry
        m = m+1;
        output{m} = plotData{i};
    end
end
% convert to a matrix, then to a cell for ease of concatenation
output = num2cell( cell2mat( output ) );
% add the dataset information
output = [output, export.Dataset];
% remove nan entries
for r = 1:size(output,1)
    for c = 1:size(output,2)
        if isnan(output{r,c})
            output{r,c} = [];
        end
    end
end

% make headers for each data column
header = cell(1,nDatCol+1);
m = 0;
totDataFields = 0;
for i = 1:size(dataFields,1)
    % see if the header info should be included
    currField = dataFields{i,1};
    if calcVal(~cellfun(@isempty,strfind(dataFields(:,1),currField)))
        % count the number of subfields
        try nSubFld = size( fieldnames( export.(currField) ), 1);
        catch, nSubFld = 1; end
        totDataFields = totDataFields+nSubFld;
        if iscell( dataFields{i,3} )
            for j = 1:nSubFld
                m = m+1;
                header{m} = [repmat(dataFields{i,3}(:,j),1,ntypes);...
                             typenames];
            end
        else
            m = m+1;
            header{m} = [repmat(dataFields(i,2),1,ntypes);...
                repmat(dataFields(i,3),1,ntypes);...
                repmat(typenames,1,nSubFld) ];
        end
    end
end

% correct data field count for the datasets info columns
totDataFields = totDataFields-1;
% Re-format header info
headerOut = cell(3,ntypes*(nDatCol+1));
m = 0;
for c = 1:nDatCol+1
    for cc = 1:ntypes
        m = m+1;
        for r =1:3
            headerOut{r,m} = header{c}{r,cc};
        end
    end
end

%% send data for plotting
% get info for titles & y-labels
labels = dataFields(calcVal(1:end-1),:);
% make figure and prepare sub-plot axis handles
screenSize =ScreenUtilities.getMonitorPositions;
screenUnits = get(0,'units');
screenSize = screenSize(1,3:4);
fgpos = [screenSize(1)*0.035, screenSize(2)*0.065, ...
    screenSize(1)*0.95 screenSize(2)*0.85];
fg = figure('units',screenUnits,'position',fgpos);
% set number of rows & columns in subplot
nrc(1) = 1; nrc(2) = 1; % [nrows ncols]
while totDataFields > nrc(1)*nrc(2)
    if nrc(1) == nrc(2), nrc(2)=nrc(2)+1;
    else nrc(1) = nrc(2);
    end
end

% package info for plotting
ax_h = nan(totDataFields,1);
m = 0;
strIn{1} = typenames;
pvals = cell(nDatCol,1);
for d = 1:size(plotData,1)
    if iscell(plotData{d})
        for subd = 1:length(plotData{d})
            m = m+1;
            strIn{2} = labels{d,3}{1,subd}; % title
            strIn{3} = labels{d,3}{2,subd}; % y-label
            plotdat = plotData{d}{subd};
            [pvals{m}, ax_h, fg] = plotCompareData(plotdat,barPosn,alphaValue,strIn,fg,ax_h,nrc,m);
            if ~isempty(pvals{m})
                pvals{m} = [ strIn{2},repmat({[]},1,length(typenames)); ...
                             [{[]},typenames]; ...
                             [typenames',num2cell(pvals{m})] ];%;...                             repmat({[]},1,length(typenames)+1) ];
            end
        end
        
    else
        m = m+1;
        strIn{2} = labels{d,2}; % title
        strIn{3} = labels{d,3}; % y-label
        plotdat = plotData{d};
        [pvals{m}, ax_h, fg] = plotCompareData(plotdat,barPosn,alphaValue,strIn,fg,ax_h,nrc,m);
        if ~isempty(pvals{m})
            pvals{m} = [ strIn{2},repmat({[]},1,length(typenames)); ...
                         [{[]},typenames]; ...
                         [typenames',num2cell(pvals{m})] ];%;...                         repmat({[]},1,length(typenames)+1) ];
        end
    end
    
    
end

%% Format the p-value matrix for writing to excel file

if ~isempty(pvals{1})
    [nRow,nCol] = size(pvals{1});
    pvalsWrite = [['p-Values from Two-sample Kolmogorov-Smirnov test (kstest2)',repmat({[]},1,nCol-1)] ;...
        repmat({[]},1,nCol);...
        cell( nDatCol*(1+nRow), nCol ) ];
    prevR = 3;
    for d = 1:nDatCol
        st = prevR;
        ed = st+nRow-1;
        pvalsWrite(st:ed,1:nCol) = pvals{d};
        prevR = ed+2;
    end
else
    pvalsWrite = [];
end

%% write Excel-readable file into a two sheets
% saveFileLoc{2} = 'Compare Cluster Metrics.xlsx';
save_xls_filename = fullfile(saveFileLoc{1},saveFileLoc{2});
if exist(save_xls_filename,'file')
    qstring = {['Overwrite exising file: ' saveFileLoc{2}]; ...
        ['in folder: ' saveFileLoc{1}]};
    choice = questdlg(qstring,...
        'Overwrite File?',...
        'Yes','No','No');
    if isempty(choice), choice = 'No'; end
    switch choice
        case 'Yes'
            delete(save_xls_filename)
        case 'No'
            [fpath,fname,fext] = fileparts( save_xls_filename );
            verN = '1';
            save_xls_filename = fullfile(fpath,[fname ' ver' verN fext]);
            % see if there's already a version number in file name
            while exist( save_xls_filename,'file' )
                verN = num2str( str2double(verN)+1 );
                save_xls_filename = fullfile(fpath,[fname ' ver' verN fext]);
            end
    end
end

% write file to excel format
warning('off','MATLAB:xlswrite:AddSheet')
xlswrite(save_xls_filename,[headerOut;output],1);
if ~isempty(pvalsWrite)
    xlswrite(save_xls_filename,pvalsWrite,2);
end
warning('on','MATLAB:xlswrite:AddSheet')

disp(['Wrote data to ' save_xls_filename])

% now save the figure generated as both fig & png
[fpath,fname,~] = fileparts(save_xls_filename);
saveas(fg,fullfile(fpath,[fname '.fig']),'fig')
saveas(fg,fullfile(fpath,[fname '.png']),'png')
disp(['Saved figure as .fig and .png to ' fullfile(fpath,fname)])

disp(' ')


end % of button selection choice
end % of function