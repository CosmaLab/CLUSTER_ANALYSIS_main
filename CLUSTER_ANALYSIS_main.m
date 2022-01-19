clear, clc % close all
% This parameter can either be a default folder to use to select a 
% .bin molecule list (a new Window will open) or a specific Insight3 file
rootdir = 'Directory';
plotstats = true; 
pix2nm = 160; % NSTORM 
% pix2nm = 157; % STORM3
% pix2nm = 107; % Martin 
% pix2nm = 98.69; %Vutara 

% set parameters
% Define the parameters to use. For a description of each parameter see the
% type: 'help FindClustersStruct' in the command window
params = defineFindClustersStruct(pix2nm);
% %%
files = Select1DataGroup('LocLists','*.bin',rootdir);
% give the computer a moment to close the GUI
pause(0.5)


% Move through all folders & perform cluster analysis
for fl = 1:size(files.data,1)
    
    params.i3file = fullfile( files.data{fl,2},files.data{fl,1} );
    
    if exist(params.i3file, 'file') == 2
        
        % perform the cluster analysis
          % ClusterResults column names:
          % 1 X(pix), 2 Y(pix), 3 Number_of_Loc_per_Cluster, 4 SigX(pix), 5 SigY(pix),
          % 6 Mean(sigX,sigY), 7 sqrt(sigX^2+sigY^2), 8 Z(nm), 9 sigZ(nm), 10 NND(nm)
          % 11 NumClusterInIsland, 12 IslandIndex
        [ClusterResults,i3,pathstr, subFolder, fileName] = FindClusters(params);
          % I3:
          % channel 0 = centroid of each cluster
          % channel 1-9 = assignment of each cluster to a channel for multi-color visualization in Insight3
        
          if ~isempty(ClusterResults)
              % Plot statistics
              if plotstats
                  
                  % save the cluster data for future comparison using
                  % callCompareClusterMetrics function
                  
                  %%2021-12-17 in saveClusterMetricData, Laura added variable 'allResults' 
                  [xynData, mask, ClustStats, ddcData, savefile] = saveClusterMetricData( params, ClusterResults, pathstr, subFolder, fileName );
                  saveClusterMetricData( params, ClusterResults, pathstr, subFolder, fileName );
                  [~, filename, ~] = fileparts(params.i3file);
                  
                  fg = plotClusterStat(filename, xynData, mask, ClustStats, ddcData, params);
                  % save the plotted data
                  saveas(fg,[savefile(1:end-4) 'Plot.fig'],'fig')
                  saveas(fg,[savefile(1:end-4) 'Plot.png'],'png')
                  
              else
                  saveClusterMetricData( params, ClusterResults, pathstr, subFolder, fileName );
                  
              end
          end
        
        disp(' ')
    end
end

disp(' ')
disp(['Finished cluster analysis for ' num2str(size(files.data,1)) ' files'])
disp(' ')