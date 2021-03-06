% function that will save information from FindClusters function for rapid
% analysis in compareClusterMetrics function
%
% INPUTS
%   params - output of FindClusterStruct function and used also for cluster
%       identification
%   ClusterResults - output of FindClusters function

function [xynData, mask, ClustStats, ddcData, savefile] = saveClusterMetricData( params, ClusterResults, pathstr, subFolder, fileName )

% calculate the area of the conventional & S.R. images generated by
% the localization list
LL = Insight3(params.i3file);
numOrigLocs = LL.numMolecules;
mask = calcMaskAreas( LL, params, true );

% perform both island-based .xyn and global .ddc nearest neighbor calculation
filexyn = create_xynFilename( params );
[xynData,ddcData,ClustStats] = extractClusterStats(filexyn,params,ClusterResults);
allResults = [ClusterResults, ddcData.nndXY];


%2021-12-17 Laura Martin added column 13 'NND[nm] Global' in the
%ClusterResults variable, now called 'allResults'
%%%%%%%%%%%%%%%%%%%%%%%%%%%
allResults = [ClusterResults, ddcData.nndXY];

%Variables
useChannel = params.use_channels;
use_iterative_segmentation = params.use_iterative_segmentation;
name = fileName;
roi = params.sum_roi_size;
threshold = params.sum_threshold;
factor = params.factor;
analysis_pixel_size = params.analysis_pixel_size;
precision = params.localization_precision;
minCluster = params.minimum_molecules_per_cluster;
max_area = params.max_segmentation_area;
NN = allResults(:,3);
SIG = allResults (:,6);
original_pixel_size = params.original_pixel_size;
findClustersStruct = params;


% [pathstr,name,~] = fileparts(i3.getFilename());
%    
% if ~isequal(exist(strcat(pathstr, '/', subFolder, '/'), 'dir'),7)
%     mkdir(strcat(pathstr, '/', subFolder, '/'));
% end

channels = sprintf('%d,',useChannel);

fprintf('Saving results to %s\\%s\\\n', pathstr,subFolder);
if use_iterative_segmentation
    fname = sprintf('%s_roi%d_th%d_fac%g_ch%s_pix%d_prec%g_min%d_area%d.xyn', name, roi, threshold, factor, channels(1:end-1), analysis_pixel_size, precision, minCluster, max_area);
else        
    fname = sprintf('%s_roi%d_th%d_fac%g_ch%s_pix%d_prec%g_min%d.xyn', name, roi, threshold, factor, channels(1:end-1), analysis_pixel_size, precision, minCluster);
end
txtFile = fullfile(pathstr, subFolder, fname);

% fopen seems to only be able to open filenames with < 256 characters
% append XXX so that it indicates that the filename is not complete
if 256 - length(txtFile) < 0
    [p, n, e] = fileparts(txtFile);
    txtFile = [p '\' n(1:end-abs(256 - length(txtFile)-3)) 'XXX' e];
end

fprintf('Writing %s... ', fname);
% write the parameters to the file
findClustersStruct.write(txtFile, '\t');



try
    % write the median N, area and density to the file        
    fle = fopen(txtFile, 'a');        
    medN = median(NN);
    original_pixel_size = original_pixel_size;
    fprintf(fle, 'medianNum\t%d\n', medN);
    medArea = pi*(median(SIG)*params.original_pixel_size)^2;
    fprintf(fle, 'medianArea[nm^2]\t%f\n', medArea);
    fprintf(fle, 'medianDensity\t%f\n', medN/medArea);
    fprintf(fle, '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', 'X[pix]','Y[pix]','NumLocalizations','sigX[pix]','sigY[pix]','sigAve[pix]','sigQuad[pix]','Z[nm]','sigZ[nm]','NND[nm] InIsland','NumClustersInIsland', 'IslandIndex', 'NND[nm] Global');
    fclose(fle);
    dlmwrite(txtFile, allResults, '-append', 'precision', '%.10g', 'delimiter', '\t')
    fprintf('DONE\n');
catch
    fprintf('\n');
    error('ERROR! Cannot save the results. Is the following file open?\n\t%s',txtFile);
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% save input & output
[sPath,sName] = fileparts(filexyn);
savefile = fullfile(sPath, [sName '.mat']);
% save the info in a .mat file
MLfileName = params.i3file;
save(savefile,...
    'MLfileName', 'mask', 'params', 'ClusterResults',...
    'numOrigLocs', 'xynData', 'ddcData', 'ClustStats')

end