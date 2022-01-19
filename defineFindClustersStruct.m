% function to facilitate user visualization of input parameters for the
% FindClusters.m algorithm
%
%
% OUTPUT
%   params - the same output as created by FindClustersStruct.m
%

function params = defineFindClustersStruct(pix2nm)

% width of the dialog fields
W = 50;

% define the parameters
params = FindClustersStruct();

% explicity alter default values
params.use_iterative_segmentation = true;
if nargin > 0 && ~isempty(pix2nm)
    params.original_pixel_size = pix2nm;
end
% params.original_pixel_size = 155; % STORM3 
% params.original_pixel_size = 160; % NSTORM 
% params.original_pixel_size = 107; % Martin 
% params.original_pixel_size = 98.69; %Vutara  


% get the field information for display & altering
fields = fieldnames(params);
% skip the first field, which is 'i3file'
fields = fields( cellfun(@isempty,strfind(fields,'i3file')) ); 
nfields = length(fields);

% check which inputs should be logical values
logicalParams = false(nfields,1);
for f = 1:nfields
    if islogical( params.(fields{f}) )
        logicalParams(f) = true;
    end
end

dlgTitle = 'Input Parameters for FindClusters.m';

defaults = cell(nfields,2);

for f = 1:nfields
    if logicalParams(f) 
        appx = '  (1/0 true/false)';
    elseif sum( strcmp( fields{f}, ...
            {'original_pixel_size', 'analysis_pixel_size', 'localization_precision'} ...
            ) )
        appx = '  [nm]';
    elseif sum( strcmp( fields{f}, ...
            {'image_width','image_height','sum_roi_size'} ...
            ) )
        appx = '  [pixels]';
    elseif strcmp( fields{f}, ...
            'use_channels' ...
            )
        appx = '  [all channels: -1]';
    elseif strcmp( fields{f}, ...
            'factor' ...
            )
        appx = '  [use: analysis_pixel / factor]';
    else
        appx = '';
    end
    defaults(f,:) = { [fields{f} appx], num2str(params.(fields{f})) };
end


nLines = [1 W];

inputs = inputdlg(defaults(:,1),dlgTitle,nLines,defaults(:,2));

for f = 1:nfields
    if logicalParams(f) 
        % input should be a logical true/false value
        params.(fields{f}) = logical( str2num( inputs{f} ) );
    else
        params.(fields{f}) = str2num( inputs{f} );
    end
        
end

