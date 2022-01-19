% XYN(path)
%
% The FindClusters program creates a XYN file. 
% Use this function to read a XYN file.
%
% Input
% -----
% path : string
%   The path to where the XYN file is located
% 
classdef XYN < handle
    properties
        filename; % the full path of the XYN file
        header; % the names of each column in 'data'
        data; % the results from the FindClusters program
        medianNum; % the median number of clusters found
        medianArea; % the median area, area = pi * median(sigAve)^2
        medianDensity; % equal to medianNum / medianArea
        medianUnit; % unit that medianArea and medianDensity are in (either nm or pixel)
        params; % the paramters used in FindClusters.m
    end
    methods
        function self = XYN(filename)
        % obj = XYN(filename)
        %
        % Load XYN files which are created from the FindClusters program
        %
        % Inputs
        % ------
        % filename : string
        %   the path of the XYN file
        %
            if ~ischar(filename)
                error('The input must be a string');
            elseif isempty(filename);
                error('The input filename is empty');
            elseif exist(filename, 'file') ~= 2
                error('"%s" does not exist', filename);
            elseif ~strcmp(filename(end-3:end), '.xyn')
                error('"%s" is not a XYN file', filename);
            end
            
            self.params = {};
            self.filename = filename;
            
            fid = fopen(filename, 'r');
            temp = strsplit(fgetl(fid),'\t');
            if strcmp(temp(1), 'X[pix]')
                self.header = temp;
                fclose(fid);
                self.data = dlmread(filename, '\t', 1, 0);
                self.medianNum = median(self.data(:,3));
                self.medianArea = pi*median(self.data(:,5))^2;
                self.medianDensity = self.medianNum/self.medianArea;
                self.medianUnit = 'pixel';
            elseif strcmp(temp(1), 'medianNum')
                self.medianNum = str2double(temp(2));
                temp = strsplit(fgetl(fid),'\t');
                self.medianArea = str2double(temp(2));
                temp = strsplit(fgetl(fid),'\t');
                self.medianDensity = str2double(temp(2));
                self.medianUnit = 'nm';
                self.header = strsplit(fgetl(fid),'\t');
                fclose(fid);
                self.data = dlmread(filename, '\t', 4, 0);
            elseif strcmp(temp(1), 'i3file')
                self.params = FindClustersStruct();
                offset = 1;
                while ~strcmp(temp(1), 'medianNum')
                    if strcmp(temp(1), 'i3file')
                        self.params.(temp{1}) = temp{2};
                    else
                        self.params.(temp{1}) = str2double(temp{2});
                    end
                    temp = strsplit(fgetl(fid),'\t');
                    offset = offset + 1;
                end
                self.medianNum = str2double(temp(2));
                temp = strsplit(fgetl(fid),'\t');
                self.medianArea = str2double(temp(2));
                temp = strsplit(fgetl(fid),'\t');
                self.medianDensity = str2double(temp(2));
                self.medianUnit = 'nm';
                self.header = strsplit(fgetl(fid),'\t');
                fclose(fid);
                self.data = dlmread(filename, '\t', offset+3, 0);
            else
                fclose(fid);
                error('Function not yet set up to handle this XYN format');
            end
        end
        
        function out = getNumClustersInEachIsland(self)
        % function counts = getNumClustersInEachIsland()
        %
        % Returns a Nx2 array of the [islandIndex NumClustersInIsland]
        % values. This function removes the duplicate values that
        % are in the xyn file.
            
            out = [];
            
            islandIndex = self.getColumnIndex('IslandIndex');
            clusterIndex = self.getColumnIndex('NumClustersInIsland');
            
            islands = self.data(:,islandIndex);
            clusters = self.data(:,clusterIndex);
            
            % walk through the islands in the forward direction
            for i = 1:length(islands)-1
                if islands(i) ~= islands(i+1)
                  out = [out; [islands(i) clusters(i)]];  
                end
            end

            % walk through in the backward direction in order to
            % get the value for the last island
            for i = length(islands):-1:2
                if islands(i) ~= islands(i-1)
                  out = [out; [islands(i) clusters(i)]];
                  break;
                end
            end
            
        end

        function idx = getColumnIndex(self, columnName)
        % function idx = getColumnIndex(columnName)
        %
        % Returns the index of the specified column name.
        %
        % Raises an error if the column name is not found within the
        % header.
        %
        % Inputs
        % ------
        % columnName : string
        %  the column name in the header to return the index of
        %
            assert(isa(columnName, 'char'), 'columnName must be a string')
            
            idx = find(strcmp(self.header, columnName));
            if isempty(idx)
                error('There is no ''%s'' column in the header of the ddc file', columnName);
            end
        end

        
    end
end