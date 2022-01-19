% Function organizing parameters required for comparing metrics obtained
% from .xyn files produced by FindClusters.m
%
% inputs = 'name', 'value' pairs
%   'datatypes'
%       character string with separate data types separated by commas
%   'startlocn' 
%       directory to start selecting .xyn files or .mat files saved
%       immediately following cluster analysis
%   'savelocn'
%       directory for saving output of cluster metric comparison
%   'savefile'
%       filename of the saved file, must be .xlsx extension
%   'alphaValue'
%       value for determining statistically significant differences between
%       data types
%   'fileExt'
%       file extension of the input data to be selected: either .xyn or
%       .mat
%   'barPosn'
%       position of the horizontal bar indicating statistically significant
%       differences between data groups, either 'top' or 'bottom' are
%       permitted
%   

function export = ... 
            compareClusterMetrics(varargin)
%%

% initialize variables
startlocn = {pwd};
startlocnDisp = startlocn{1};
startlocnIn = [];
dataTypes = 'Type 1, Type 2';
saveFileLoc{1} = startlocn{1};
saveFileLoc{2} = 'Compare Cluster Metrics.xlsx';
alphaValue = 0.05;
FileExtension = 'xyn';
barPosn = 'top';
for i = 1:2:nargin
    
    if strcmpi(varargin{i},'datatypes') || ~isempty(strfind( varargin{i},'type' ))
        uiType = varargin{i+1};
        if ischar(uiType) 
            dataTypes = uiType;
        end
    
    elseif strcmpi(varargin{i},'startlocn') || ~isempty(strfind( varargin{i},'start' ))
        uiPath = varargin{i+1};
        if ischar(uiPath) && exist(uiPath,'dir')
            startlocn = {uiPath};
        elseif iscell(uiPath)
            startlocn = cell(length(uiPath),1);
            direxist = 0;
            for n = 1:length(uiPath)
                if ischar(uiPath{n}) && exist(uiPath{n},'dir')
                    startlocn{n} = uiPath{n};
                    direxist = direxist+1;
                end
            end
            if direxist == length(uiPath) % all directories exist
                startlocnIn = startlocn;
            end
        else
            startlocn = {pwd};
        end
        if ~isempty(startlocnIn)
            startlocnDisp = 'Input cell array of folders accepted';
        end        
        
        
    elseif ~isempty(strfind( varargin{i},'save' ))
        uiSave = varargin{i+1};
        if ischar(varargin{i}) && ischar(uiSave)
            if ~isempty(strfind( varargin{i},'loc' ) ) % presume save location
                if exist(uiSave,'dir')
                    saveFileLoc{1} = uiSave;
                end
            elseif ~isempty(strfind( varargin{i},'file' ) ) ||...
                    ~isempty(strfind( varargin{i},'File' ) ) % presume save file name
                saveFileLoc{2} = uiSave;
            end
        end
            
    elseif strcmpi(varargin{i},'alphaValue') || ~isempty(strfind( varargin{i},'alpha' ))
        uiAlpha = varargin{i+1};
        if isa(uiAlpha,'double') 
            if uiAlpha > 1
                alphaValue = uiAlpha;
            else
                alphaValue = uiAlpha/100;
            end
        end
        
    elseif strcmpi(varargin{i},'fileExt') || ~isempty(strfind( varargin{i},'ext' ))...
            || ~isempty(strfind( varargin{i},'Ext' ))
        uiExt = varargin{i+1};
        if ischar(uiExt) && sum( strcmpi(uiExt,{'xyn','.xyn','mat','.mat'}) )==1
            FileExtension = uiExt;
        end
        
    elseif strcmpi(varargin{i},'barPosn') || ~isempty(strfind( varargin{i},'bar' ))
        uiBar = varargin{i+1};
        if ischar(uiBar) && sum( strcmpi(uiBar,{'top','bottom','bot'}) )==1
            barPosn = uiBar;
        end
        
    end % varargin if statement
end % nargin for-loop




dlgTitle = 'Comparing Clusters Metadata';
defaults = ...
    {dataTypes,	 'Type the Data Type names separated by a comma';... 1
    FileExtension,     'What type of file will you select, .xyn or .mat?';...             2
    barPosn, 'Type ''top'' or ''bottom'' for showing Stat Sig differences when ploting';... 3
    num2str(alphaValue), 'p-value for considering statistically siginificant differences';...	4           5
    startlocnDisp, 'Folder for data selection';... 5
    saveFileLoc{1},'Folder to save data & figures';... 6
    saveFileLoc{2},'Name of the .xlsx file for saving data'};%            7
W = 70;
nLines = ...
    [repmat([1 W],size(defaults,1)-3,1); ...
    repmat([2 W],2,1); 1 W];

inputs = inputdlg(defaults(:,2),dlgTitle,nLines,defaults(:,1));


%% parce input values

if ~isempty(inputs)
    % type definition
    clear type
    n = 1;
    commaIdx = strfind(inputs{n},',');
%     if ~ismember(commaIdx,1), commaIdx = [1,commaIdx];end
    commaIdx = [commaIdx, size(inputs{n},2)];
    prevIdx = 1;
    ncommas =length(commaIdx);
    for i = 1:ncommas
        val = inputs{n}(prevIdx:commaIdx(i));
%         wspIdx = ~isstrprop( val,'wspace');
        val = val(find(~isstrprop( val,'wspace'),1,'first'):end);% remove white space at start
        val = val( val~=',' ); % remove commas 
        type.(['t' num2str(i)]) = val;
        if i < ncommas
            prevIdx = commaIdx(i)+1;
        end
    end
    
    
    % extension of files to be selected
    n = 2;
    FileExtension = inputs{n};
    if ~sum( strcmp(FileExtension,{'xyn','.xyn','mat','.mat'}) )
        FileExtension = 'xyn';
    elseif strcmp(FileExtension(1),'.')
        FileExtension = FileExtension(2:end);
    end
    
    
    % position of the Stat Sig bar
    n = 3;
    barPosn = inputs{n};
    if strcmpi(barPosn,'top')
        barPosn = 'top'; % ensure it's lower case
    elseif strcmpi(barPosn,'bottom') || ~isempty(strfind(barPosn,'b'))
        % accept 'bottom' if they hit 'b', since it's easier to mis-spell
        barPosn = 'bottom'; % ensure it's lower case
    else
        barPosn = 'top';
    end
    
    
    % alpha value
    n = 4;
    alphaValue = str2double(inputs{n});
    if alphaValue > 1, 
        disp(['Input alpha value = ' num2str(alphaValue) ...
            ', reset to ' num2str(alphaValue/100)])
        alphaValue = alphaValue/100; 
    end
    
    
    % starting directory
    n = 5;
    if isempty(startlocnIn) % nothing recognized input
        if exist(inputs{n},'dir')
            startlocn{1} = inputs{n};
        else
            startlocn{1} = {pwd};
        end
    else
        startlocn = startlocnIn;
    end
    
    
    % save file directory
    n = 6;
    if exist(inputs{n},'dir')
        saveFileLoc{1} = inputs{n};
    else
        saveFileLoc{1} = startlocn{1};
    end
    
    % save file name
    n = 7;
    [~,fname,ext] = fileparts( inputs{n} );
    if strcmp(ext,'.xlsx')
        saveFileLoc{2} = inputs{n};
    else
        saveFileLoc{2} = [fname '.xlsx'];
    end

    
    % send parameters into cluster comparison calculator
    export = ...
        calc_plot_ClusterMedians(type,startlocn,FileExtension,saveFileLoc,barPosn,alphaValue);
else
    export = [];
    disp('Calculation cancelled')
end
% type,startlocn,FileExtension,saveFileLoc,barPosn,alphaValue
