% DDC(path)
%
% The DistanceDualColor program creates a DDC file. 
% Use this function to read a DDC file.
%
% Input
% -----
% path : string
%   The path to where the DDC file is located
% 
classdef DDC < handle
    properties
        filename; % the full path of this DDC file
        file1; % the full path of the first XYN file
        file2; % the full path of the second XYN file
        range; % the range, in nm, that was used to separate the different distance categories
        rangeHeader;
        numClusters1; % the total number of clusters in file1
        numClusters2; % the total number of clusters in file2
        percentages; % the percentage of the number of clusters in each distance range
        header; % the names of each column in 'data' and 'medians'
        medians; % 2D array of the median values in each distance range
        data; % struct of the data for each distance range
    end
    methods
        function self = DDC(filename)
        % obj = DDC(filename)
        %
        % Load DDC files which are created from the DistanceDualColor program
        %
        % Inputs
        % ------
        % filename : string
        %   the path of the DDC file
        %
            if ~ischar(filename)
                error('The input must be a string');
            elseif isempty(filename);
                error('The input filename is empty');
            elseif exist(filename, 'file') ~= 2
                error('"%s" does not exist', filename);
            elseif ~strcmp(filename(end-3:end), '.ddc')
                error('"%s" is not a DDC file', filename);
            end
            
            % the DDC filename
            self.filename = filename;
            
            fid = fopen(filename, 'r');
            
            % get the first file name
            temp = strsplit(fgetl(fid),'\t');
            self.file1 = temp{2};

            % get the second file name
            temp = strsplit(fgetl(fid),'\t');
            self.file2 = temp{2};

            % get the distance range
            temp = strsplit(fgetl(fid),'\t');
            self.range = [];
            self.rangeHeader = {};
            for i=1:length(temp)-2
                self.rangeHeader{i} = temp{i};
                r = textscan(temp{i}, '%% in %d-%d nm');
                if r{2} ~= 2147483647
                    self.range = [self.range, r{2}];
                end
            end
            
            % get the percentages and the total number of clusters
            temp = strsplit(fgetl(fid),'\t');
            self.percentages = [];
            for i=1:length(self.range)+3
                val = str2double(temp{i});
                if i <= length(self.range) + 1
                    self.percentages = [self.percentages, val];
                elseif i == length(self.range) + 2
                    self.numClusters1 = val;
                else
                    self.numClusters2 = val;
                end
            end
            
            % get the header
            temp = strsplit(fgetl(fid),'\t');
            self.header = temp;
            
            % get the median data for each range
            self.medians = [];
            for i=1:length(self.range)+2
                temp = strsplit(fgetl(fid),'\t');
                self.medians = [self.medians; str2double(temp(1:end-1))];
            end
            
            % skip another header row
            fgetl(fid);
            
            % get the data
            self.data = {};            
            for i=1:length(self.range)+1
                self.data{i} = [];
                temp = strsplit(fgetl(fid),'\t');
                while length(temp) > 3
                    self.data{i} = [self.data{i}; str2double(temp)];
                    tline = fgetl(fid);
                    if ischar(tline)
                        temp = strsplit(tline,'\t');
                    else
                        temp = '';
                    end
                end
            end

            fclose(fid);
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
