% Reads and writes Insight3 molecule lists
%
% Read a molecule list
% i3 = Insight3('D:/sample.bin');
% % to access the data in the molecule list you can use the following in your program
% i3.data;
% % or you can create a copy of the molecule list and use this new variable
% mList = i3.getData();
%
% Create a new molecule list with 1000 random x,y localizations for a 256 x 256 image
% i3 = Insight3();
% i3.setData( 256*rand(1000,2) );
% i3.write('D:/newMoleculeList.bin');
%
classdef Insight3 < handle    
    properties
        numMolecules; % int, number of molecules
        numFrames; % int, number of frames
        filename; % string, the path of the Insight3 file
        data; % 2D matrix, the molecule list
        columnNames; % list of strings describing what each column in the 'data' matrix represents
    end
    properties(Hidden)
        % assign the default parameters
        txtHeader = {'channel','x','y','xc','yc','height','area','width','phi','aspect','background','intensity','frame','trackLength','link','valid','z','zc'};
        binHeader = {'x','y','xc','yc','height','area','width','phi','aspect','background','intensity','channel','fitIterations','frame','trackLength','link','z','zc'};
        defaultX = 0.0;
        defaultXc = 0.0;
        defaultY = 0.0;
        defaultYc = 0.0;
        defaultHeight = 100.0; 
        defaultArea = 10000.0;
        defaultWidth = 300.0;
        defaultPhi = 0.0;
        defaultAspect = 1.0;
        defaultBackground = 0.0;
        defaultIntensity = 10000.0;
        defaultChannel = 1;
        defaultFitIter = 1;
        defaultValid = 1;
        defaultFrame = 1;
        defaultTrackLength = 1;
        defaultLink = -1;
        defaultZ = 0.0;
        defaultZc = 0.0;
        verbose; % boolean, indicates whether to display message to the command window
        forceOverwrite = false; % boolean, overwrite the file even if it exists
    end
    methods        
        function self = Insight3(filename, verbose)
        % function self = Insight3(filename, verbose)
        %
        % Only 'bin' and 'txt' filenames are supported
        %
        % Inputs
        % ------
        % filename : string
        %   the path of the Insight3 file
        %   you can set filename to the strings 'txt' or 'bin' if you do not
        %   care about reading/writing an actual file at when creating an Insight3 object
        %
        % verbose (optional argument) : boolean
        %   display messages to the Command Window, default is true
            
            if nargin == 0
               %error('ERROR Insight3() :: you must specify a .bin or a .txt filename.\nIf you do not want to specify a file then just pass ''%s'' or ''%s'' in as a string\nso that Matlab knows how to structure the columns of the molecule list.', 'bin', 'txt');
               self.verbose = true;
               filename = 'dummy.bin';               
            elseif nargin < 2
                self.verbose = true;
            else
                self.verbose = verbose;
            end
            
            self.numMolecules = 0;
            self.numFrames = 0;
            self.data = [];            
            self.setFilename(filename);
            
            % using the java Insight3io.jar file to read/write a molecule list from/to 
            % a file is significantly faster than using Matlab to perform the task
            path = fullfile(fileparts(mfilename('fullpath')), 'Insight3io.jar');
            if isempty(strfind(javaclasspath, path))
                javaaddpath(path);
            end
            
            % if this file already exists then we probably want to read it
            if exist(self.filename, 'file') == 2                
                self.read();
            end
        end

        function read(self)
        % function read()
        %
        % Reads the molecule list into a Matlab matrix (self.data)
        %
            self.data = [];
        
            % make sure this Insight3 file exists
            if exist(self.filename, 'file') ~= 2
                error('ERROR Insight3.read() :: Insight3 file does not exist\n\t%s', self.filename)
            end
            
            try
                self.data = Insight3IO.read(self.filename, self.verbose);
                self.numFrames = Insight3IO.nFrames;
                self.numMolecules = Insight3IO.nMolecules;
            catch
                self.readSLOW()
            end

        end

        function write(self, filename)
        % function write(filename)
        %
        % Writes the values in 'data' to the current 'filename'
        %
        % Inputs
        % ------
        % filename (optional) : string 
        %   the filename name of the file to create
            
            if nargin == 2
                self.setFilename(filename);
            end
            
            if isempty(self.data)
                if self.verbose
                    fprintf(2, 'WARNING Insight3.write() :: The molecule list is empty. Not creating %s\n', self.filename);
                end
                return;
            end
        
            % check to see if you are going to overwrite a file
            if ~self.forceOverwrite && exist(self.filename, 'file') == 2
                error('ERROR Insight3.write() :: Insight3 molecule list exists. Not going to overwrite\n\t%s\nYou can call forceFileOverwrite(true) to always overwite a file', self.filename);
            end
            
            [~,~,ex] = fileparts(self.filename);            
            if isempty(ex)
                ex = strcat('.', self.getFormat());
                self.filename = strcat(self.filename, ex);
            end
            
            try
                Insight3IO.write(self.filename, self.data, self.numFrames, self.verbose);
            catch
                self.writeSLOW()
            end
            
        end
        
        function success = convertFormat(self)
        % function success = convertFormat()
        %
        % Converts the molecule list to the other file format. 
        %
        % For example,
        %
        %   If the molecule list is currently in a .bin file format then 
        %   convert it to be in a .txt format
        %   WARNING if you convert bin to txt you will loose the information 
        %           in the 'FitIterations' column and the 'Valid' column
        %           will be populated with 1's
        %
        %   If the molecule list is currently in a .txt file format then 
        %   convert it to be in a .bin format
        %   WARNING if you convert txt to bin you will loose the information 
        %           in the 'Valid' column and the 'FitIterations' column
        %           will be populated with 1's
        %
        % Returns true if the conversion was successful
            
            success = ~isempty(self.data);
        
            if ~success
                error('ERROR Insight3.convertFormat() :: The molecule list is empty');
            end
        
            if strcmp(self.getFormat(), 'bin')
                if self.verbose
                    fprintf('Converting the %s molecule list to a %s molecule list\n', '.bin', '.txt');
                end
                self.columnNames = self.txtHeader;
                self.filename = strcat(self.filename(1:end-3), 'txt');
                self.data(:,:) = self.data(:,[12 1 2 3 4 5 6 7 8 9 10 11 14 15 16 13 17 18]);
                self.data(:,16) = 1; % set the 'Valid' column to be 1's
            else
                if self.verbose
                    fprintf('Converting the %s molecule list to a %s molecule list\n', '.txt', '.bin');
                end
                self.columnNames = self.binHeader;
                self.filename = strcat(self.filename(1:end-3), 'bin');
                self.data(:,:) = self.data(:,[2 3 4 5 6 7 8 9 10 11 12 1 16 13 14 15 17 18]);
                self.data(:,13) = 1; % set the 'FitIterations' column to be 1's
            end
            
        end
        
        function bool = isBinFormat(self)
        % function bool = isBinFormat()
        %
        % Returns true if the current format type is 'bin'
        % Returns false if the current formt type is 'txt'
            bool = strcmp(self.getFormat(), 'bin');
        end
        
        function extn = getFormat(self)
        % function extn = getFormat()
        %
        % Returns the file format that the molecule list is currently in.
        % Two return types are possible, 'bin' or 'txt'
            extn = self.filename(end-2:end);
        end
        
        function forceFileOverwrite(self, bool)
        % function forceFileOverwrite(bool)
        %
        % By default the molecule list you write will not overwrite a file if it
        % already exists. 
        %
        % Inputs
        % ------
        % bool : boolean (logical)
        %   'true' - always overwrite the dax file
        %   'false' - never overwrite the dax file
        %
            self.forceOverwrite = bool;
        end

        function nMolecules = getNumberOfMolecules(self)
        % function nMolecules = getNumberOfMolecules()
        %
        % Returns the number of molecules
            nMolecules = self.numMolecules;
        end

        function XY = getXY(self)
        % function XY = getXY()
        %
        % Returns the X,Y localization values
            if self.isBinFormat()
                XY = self.data(:,1:2);
            else
                XY = self.data(:,2:3);
            end
        end
                
        function XYc = getXYcorr(self)
        % function XYc = getXYcorr()
        %
        % Returns the drift-corrected X,Y localization values
            if self.isBinFormat()
                XYc = self.data(:,3:4);
            else
                XYc = self.data(:,4:5);
            end
        end

        function XYZ = getXYZ(self)
        % function XYZ = getXYZ()
        %
        % Returns the X,Y localization values and the Z value
            if self.isBinFormat()
                XYZ = self.data(:,[1 2 17]);
            else
                XYZ = self.data(:,[2 3 17]);
            end
        end

        function XYZc = getXYZcorr(self)
        % function XYZc = getXYZcorr()
        %
        % Returns the drift-corrected X,Y localization values and the drift-corrected Z value
            if self.isBinFormat()
                XYZc = self.data(:,[3 4 18]);
            else
                XYZc = self.data(:,[4 5 18]);
            end
        end
        
        function fname = getFilename(self)
        % function fname = getFilename()
        %
        % Returns the full path of the Insight3 file
            fname = self.filename;
        end

        function success = setFilename(self, filename)
        % function success = setFilename(filename)
        %
        % Sets the full path of the Insight3 file
        %
        % Inputs
        % ------
        % filename : string
        %   the name of the file to associate with this Insight3 object
        
            if length(filename) < 3
                error('ERROR Insight3.setFilename() :: Filename too short');
            end
            
            extn = filename(end-2:end);
            success = strcmp(extn, 'bin') || strcmp(extn, 'txt');
            if ~success
                error('ERROR Insight3.setFilename() :: Only .txt and .bin extensions are allowed');
            end
            
            % check if the format should change with this new filename
            if ~isempty(self.data) && ~strcmp(self.getFormat(), extn)
                self.filename = strcat(filename(1:end-3), self.getFormat());
                success = self.convertFormat();
            else
                self.filename = filename;
                if strcmp(extn, 'bin')
                    self.columnNames = self.binHeader;
                else
                    self.columnNames = self.txtHeader;
                end
            end
        end
        
        function idx = getColumnIndex(self, columnName)
        % function idx = getColumnIndex(columnName)
        %
        % The 'columnName' must be one of values in the 'columnNames' array
        % see getColumnNames()
        %
        % Inputs
        % ------
        % columnName : string
        %   The name of the column 
        %
        % Returns
        % -------
        % idx : integer,
        %   the column number in the molecule list, 
        %
            i = 1;
            while i < 19
                if strcmp(self.columnNames{i},columnName)
                    idx = i;
                    return
                end
                i = i + 1;
            end
            error('ERROR! Insight3.getColumnIndex() :: There is no ''%s'' column', columnName);
        end

        function header = getColumnNames(self)
        % function header = getColumnNames()
        %
        % Returns a 1D cell of strings describing what each column in the 'data' matrix is
            header = self.columnNames;
        end
        
        function data = getData(self)
        % function data = getData()
        %
        % Returns a copy of the molecule list
            data = self.data;
        end

        function setData(self, data)
        % function setData(data)
        %
        % Sets the molecule list to be equal to 'data'
        %
        % Inputs
        % ------
        % data : 2D matrix
        %   - can be a Nx2 matrix, where N is the number of molecules and the 2 columns are x, y values
        %   - can be a Nx3 matrix, where N is the number of molecules and the 3 columns are x, y, z values
        %   - can be a Nx18 matrix, where data contains all information
        %  The columns that you don't specify in 'data' get set to a default value
        %
            
            [self.numMolecules, cols] = size(data);
            
            % ensure that the number of columns in 'data' makes sense
            if ~ismember(cols , [2 3 18])
                error('ERROR Insight3.setData() :: The data can only have 2, 3 or 18 columns. Got %d columns.\nEither pass only x,y values, x,y,z values, or an entire molecule list', cols);
            end
    
            if strcmp(self.getFormat(), 'txt')           
                if cols == 18
                    self.numFrames = max(data(:,13));
                    self.data = data;
                else
                    self.numFrames = 1;
                    self.data = zeros(self.numMolecules,18);
                    self.data(:,1) = self.defaultChannel;
                    self.data(:,2) = data(:,1);
                    self.data(:,3) = data(:,2);
                    self.data(:,4) = data(:,1);
                    self.data(:,5) = data(:,2);
                    self.data(:,6) = self.defaultHeight;
                    self.data(:,7) = self.defaultArea;
                    self.data(:,8) = self.defaultWidth;
                    self.data(:,9) = self.defaultPhi;
                    self.data(:,10) = self.defaultAspect;
                    self.data(:,11) = self.defaultBackground;
                    self.data(:,12) = self.defaultIntensity;
                    self.data(:,13) = self.defaultFrame;
                    self.data(:,14) = self.defaultTrackLength;
                    self.data(:,15) = self.defaultLink;
                    self.data(:,16) = self.defaultValid;
                    if cols == 2                        
                        self.data(:,17) = self.defaultZ;
                        self.data(:,18) = self.defaultZc;
                    else
                        self.data(:,17) = data(:,3); 
                        self.data(:,18) = data(:,3);
                    end
                end
            else
                if cols == 18
                    self.numFrames = max(data(:,14));
                    self.data = data;
                else
                    self.numFrames = 1;
                    self.data = zeros(self.numMolecules,18);
                    self.data(:,1) = data(:,1);
                    self.data(:,2) = data(:,2);
                    self.data(:,3) = data(:,1);
                    self.data(:,4) = data(:,2);
                    self.data(:,5) = self.defaultHeight;
                    self.data(:,6) = self.defaultArea;
                    self.data(:,7) = self.defaultWidth;
                    self.data(:,8) = self.defaultPhi;
                    self.data(:,9) = self.defaultAspect;
                    self.data(:,10) = self.defaultBackground;
                    self.data(:,11) = self.defaultIntensity;
                    self.data(:,12) = self.defaultChannel;
                    self.data(:,13) = self.defaultFitIter;   
                    self.data(:,14) = self.defaultFrame;
                    self.data(:,15) = self.defaultTrackLength;
                    self.data(:,16) = self.defaultLink;
                    if cols == 2                        
                        self.data(:,17) = self.defaultZ;
                        self.data(:,18) = self.defaultZc;
                    else
                        self.data(:,17) = data(:,3); 
                        self.data(:,18) = data(:,3);
                    end
                end
            end
        end

        function setVerbose(self, verbose)
        % function setVerbose(verbose)
        %
        % Set whether to display messages to the command window
        %
        % Inputs
        % ------
        % verbose : boolean
        %  if 'true' then display warning/error messages to the screen
        %  
            self.verbose = verbose;
        end

        function createDefault(self, numMolecules)
        % function createDefault(numMolecules)
        %
        % Creates a new molecule list with 'numMolecules' molecules, each
        % initialized to a default value
        %
        % Inputs
        % ------
        % numMolecules : integer
        %   the number of molecules to create
        
            self.numMolecules = max(1, numMolecules);
            self.numFrames = 1;
            self.data = zeros(self.numMolecules,18);
            if strcmp(self.getFormat(), 'txt')           
                self.data(:,1) = self.defaultChannel;
                self.data(:,2) = self.defaultX;
                self.data(:,3) = self.defaultY;
                self.data(:,4) = self.defaultXc;
                self.data(:,5) = self.defaultYc;
                self.data(:,6) = self.defaultHeight;
                self.data(:,7) = self.defaultArea;
                self.data(:,8) = self.defaultWidth;
                self.data(:,9) = self.defaultPhi;
                self.data(:,10) = self.defaultAspect;
                self.data(:,11) = self.defaultBackground;
                self.data(:,12) = self.defaultIntensity;
                self.data(:,13) = self.defaultFrame;
                self.data(:,14) = self.defaultTrackLength;
                self.data(:,15) = self.defaultLink;
                self.data(:,16) = self.defaultValid;
                self.data(:,17) = self.defaultZ;
                self.data(:,18) = self.defaultZc;
            else
                self.data(:,1) = self.defaultX;
                self.data(:,2) = self.defaultY;
                self.data(:,3) = self.defaultXc;
                self.data(:,4) = self.defaultYc;
                self.data(:,5) = self.defaultHeight;
                self.data(:,6) = self.defaultArea;
                self.data(:,7) = self.defaultWidth;
                self.data(:,8) = self.defaultPhi;
                self.data(:,9) = self.defaultAspect;
                self.data(:,10) = self.defaultBackground;
                self.data(:,11) = self.defaultIntensity;
                self.data(:,12) = self.defaultChannel;
                self.data(:,13) = self.defaultFitIter;   
                self.data(:,14) = self.defaultFrame;
                self.data(:,15) = self.defaultTrackLength;
                self.data(:,16) = self.defaultLink;
                self.data(:,17) = self.defaultZ;
                self.data(:,18) = self.defaultZc;
            end
            
        end
        
        function [data, indices] = getChannels(self, channels)
        % function [data, indices] = getChannels(channels)
        %
        % Returns all molecules in the specified channel(s)
        %
        % Inputs
        % ------
        % channels : integer OR an array of channels
        %    eg. 0 <= channel <= 9  
        %    eg. [1 3] will return the molecules in channels 1 and 3
        % 
        % Returns
        % -------
        % data - 2D array
        %   the molecules list for the specified channels
        %
        % indices - 1D array
        %   the rows indices of the entire molecule list that each row in 'data' corresponds to
             
            % make sure all channel values are valid
            %i = 1;
            %while i <= length(channels) 
            %    if (channels(i) < 0) || (channels(i) > 9)
            %        channels(i) = []; % remove
            %    else
            %        i = i + 1;
            %    end
            %end
            
            if isempty(channels)
                data = [];
                indices = [];
            else
                idx = self.getColumnIndex('channel');
                members = ismember(self.data(:,idx), channels);
                temp = transpose(1:self.numMolecules);
                data = self.data(members,:);                
                indices = temp(members);
            end                
            
        end

        function [data, indices] = getFrames(self, frames)
        % function [data, indices] = getFrames(frames)
        %
        % Returns all molecules in the specified frames(s)
        %
        % Inputs
        % ------
        % frames : integer OR an array of frames
        %    eg. frames = 1:4:101 will return all the molecules in frames 1, 5, 9, 13, 17, ..., 101
        %        frames = [10 20 123] will return all the molecules in frames 10, 20 and 123 
        %
        % Returns
        % -------
        % data - 2D array
        %   the molecules list for the specified frame(s)
        %
        % indices - 1D array
        %   the rows indices of the entire molecule list that each row in 'data' corresponds to
             
            if isempty(frames)
                data = [];
                indices = [];
            else
                idx = self.getColumnIndex('frame');
                members = ismember(self.data(:,idx), frames);
                temp = transpose(1:self.numMolecules);
                data = self.data(members,:);                
                indices = temp(members);
            end
        end
        
        function photons = getPhotons(self, conversionFactor)
        % function photons = getPhotons(conversionFactor)
        %
        % Returns the number of photons for each molecule in the molecule list
        %
        % Inputs
        % ------
        % conversionFactor : double
        %    The photon conversion factor, eg. 0.41 for STORM and 0.14 for NSTORM
        % 
            idx = self.getColumnIndex('area');
            photons = conversionFactor * self.data(:,idx);
        end
        
        function filter(self, varargin)
        % function filter(varargin)
        %
        % Removes the molecules from the molecule list that are not within the 
        % specified photon-number, track-length and channel criteria.
        %
        % Equivalent to only saving the high-resolution molecule list as a .txt 
        % file in Insight3.        
        %
        % Allowed keys are 'minPhotons', 'maxPhotons', 'maxTrackLength', 'keepChannels'
        %
        % The only required input parameter is the photonFactor. All other values 
        % will be set to their default value if you don't specify them.
        %
        % Examples: 
        % i3.filter(0.41)
        % i3.filter(0.14, 'maxTrackLength', 10)
        %
        % The default values are:
        %   minPhotons = 500
        %   maxPhotons = 1e4
        %   maxTrackLength = 5
        %   keepChannels = 0:3 (ie. keep channels 0, 1, 2, 3)
        %
        % Inputs
        % ------
        % photonFactor : double
        %     The photon conversion factor, eg. 0.41 for STORM and 0.14 for N-STORM
        %     The value is valid in the range 0 < photon_factor < 1
        %
        % minPhotons (optional parameter) : double
        %     The minimum number of photons, default 500
        %
        % maxPhotons (optional parameter) : double
        %     The maximum number of photons, default 10000
        %
        % maxTrackLength (optional parameter) : integer
        %     The maximum track length, default 5
        %
        % keepChannels (optional parameter) : integer OR an array of integers
        %    Keep only the molecules that are with the specified channels, default 0:3
        %
            if nargin == 1
                error('ERROR! Insight3.filter() :: You must specify the photon conversion factor\nFor example, i3.filter(0.41)');
            end
            
            photonFactor = varargin{1};
            if ~(isa(photonFactor,'double') && photonFactor > 0 && photonFactor <= 1)
                error('ERROR! Insight3.filter() :: The first value must be the photon\nconversion factor and it must be a number between 0 and 1')
            end                
            
            % parse the input params
            opts   = {'minPhotons','maxPhotons','maxTrackLength','keepChannels'};
            values = {         500,         1e4,               5,           0:3};
            dict = containers.Map(opts,values);
            for i=2:2:nargin-1
                if isKey(dict, varargin{i})
                    dict(varargin{i}) = varargin{i+1};
                end
            end
            
            photons = self.getPhotons(photonFactor);
            self.data = self.data(photons >= dict('minPhotons') & photons <= dict('maxPhotons') ,:);
            
            idx = self.getColumnIndex('trackLength');
            self.data = self.data(self.data(:,idx) <= dict('maxTrackLength') ,:);
            
            [self.data, ~] = self.getChannels(dict('keepChannels'));
            
            idx = self.getColumnIndex('frame');
            self.numFrames = max(self.data(:,idx));
            
            self.numMolecules = size(self.data, 1);
        end
        
        function data = getColumn(self, columnName)
        % function data = getColumn(columnName)
        %
        % Returns the data in the specified columnName
        % see getCoumnsNames
        %
        % Inputs
        % ------
        % column_name : string
        %   The name of the column
        %
            data = self.data(:,self.getColumnIndex(columnName));
        end        
        
        function writeViSP(self, nm_per_pixel, dimension)
        % function writeViSP(nm_per_pixel, dim)
        %
        % Creates a 2D or 3D file for rendering the molecule list in ViSP
        % http://www.nature.com/nmeth/journal/v10/n8/full/nmeth.2566.html
        %
        % Currently, does not support creating the .2dlp or .3dlp files
        %
        % Inputs
        % ------
        % nm_per_pixel : double
        %    The number of nanometers in 1 pixel
        %
        % dimension : int
        %   Specify whether you want to create a 2D (dimension=2) or 3D 
        %   (dimension=3) ViSP file
        %
            if (dimension ~= 2) && (dimension ~= 3)
                dimension = 2; % use 2D by default
            end
            visp = zeros(self.numMolecules, dimension+2);
            visp(:,1) = self.getColumn('xc')*nm_per_pixel;
            visp(:,2) = self.getColumn('yc')*nm_per_pixel;
            visp(:,dimension+1) = self.getColumn('area');
            visp(:,dimension+2) = self.getColumn('frame');
            if dimension == 3
                visp(:,3) = self.getColumn('zc');
            end
            fout = [self.filename(1:end-3) num2str(dimension) 'd'];
            dlmwrite(fout, visp, 'delimiter', '\t');
            if self.verbose
                fprintf('Created %s\n', fout);
            end
        end

        function writeSharpViSu(self, nm_per_pixel, photon_factor)
        % function writeSharpViSu(nm_per_pixel, photon_factor)
        %
        % Save the molecule list in the SharpViSu format
        % http://bioinformatics.oxfordjournals.org/content/early/2016/03/26/bioinformatics.btw123
        %
        % Inputs
        % ------
        % nm_per_pixel : double
        %    The number of nanometers in 1 pixel
        %
        % photon_factor : double
        %     The photon conversion factor, eg. 0.41 for oSTORM and 0.14 for N-STORM
        %     The value is valid in the range 0 < photon_factor < 1
        %
            visu1 = zeros(self.numMolecules, 9);
            visu1(:,2) = self.getColumn('frame');
            visu1(:,4) = self.getColumn('xc') * nm_per_pixel;
            visu1(:,5) = self.getColumn('yc') * nm_per_pixel;
            visu1(:,6) = self.getColumn('zc');
            visu1(:,7) = self.getPhotons(photon_factor);
            visu1(:,8) = self.getColumn('width') ./ self.getColumn('aspect');
            visu1(:,9) = self.getColumn('width') .* self.getColumn('aspect');
            
            visu = sortrows(visu1, 2);
            
            % set the eventID's in column 3
            fr = visu(1,2);
            ev = 1;
            for i = 1:size(visu,1)
                if visu(i,2) == fr
                   visu(i,3) = ev;
                   ev = ev + 1;
                else
                   ev = 1;
                   fr = visu(i,2);
                   visu(i,3) = ev;
                   ev = ev + 1;
                end
            end

            fout = [self.filename(1:end-3) 'ascii'];
            fid = fopen(fout, 'w');
            fprintf(fid, '#stackID(UINT32), frameID(UINT32), eventID(UINT32), x0(float), y0(float), z0(float), photon_count0(float), sigmaX0(float), sigmaY0(float)\r\n');
            fclose(fid);
            dlmwrite(fout, visu, '-append', 'delimiter', ',', 'precision', 6);
            if self.verbose
                fprintf('Created %s\n', fout);
            end
        end

        function writeSRTesseler(self)
        % function writeSRTesseler()
        %
        % Save the molecule list in the SR-Tesseler format
        % http://www.nature.com/nmeth/journal/v12/n11/full/nmeth.3579.html
        %
            tess = zeros(self.numMolecules, 4);
            tess(:,1) = self.getColumn('xc');
            tess(:,2) = self.getColumn('yc');
            tess(:,3) = self.getColumn('area');
            tess(:,4) = self.getColumn('frame') - 1;
            
            fout = [self.filename(1:end-4) '_SR_Tesseler.txt'];
            fid = fopen(fout, 'w');
            fprintf(fid, '%d\t%d\r\n', self.numFrames, self.numMolecules);
            fclose(fid);
            dlmwrite(fout, tess, '-append', 'delimiter', '\t', 'precision', 6);
            if self.verbose
                fprintf('Created %s\n', fout);
            end
        end
        
        function writeRapidSTORM(self, nm_per_pixel)
        % function writeRapidSTORM(nm_per_pixel)
        %
        % Save the molecule list in the RapidSTORM format
        %
        % Inputs
        % ------
        % nm_per_pixel : double
        %    The number of nanometers in 1 pixel
        %
            rstorm = zeros(self.numMolecules, 8);
            rstorm(:,1) = self.getColumn('xc') * nm_per_pixel;
            rstorm(:,2) = self.getColumn('yc') * nm_per_pixel;
            rstorm(:,3) = self.getColumn('frame');
            rstorm(:,4) = self.getColumn('area');
            rstorm(:,5) = self.getColumn('width') ./ self.getColumn('aspect');
            rstorm(:,6) = self.getColumn('width') .* self.getColumn('aspect');
            rstorm(:,8) = self.getColumn('background');

            fout = [self.filename(1:end-4) '_RapidSTORM.txt'];
            fid = fopen(fout, 'w');
            fprintf(fid, '# <localizations insequence="true" repetitions="variable"><field identifier="Position-0-0" syntax="floating point with . for decimals and optional scientific e-notation" semantic="position in sample space in x dimension" unit="nanometer" min="0 m" max="%.4e m" /><field identifier="Position-1-0" syntax="floating point with . for decimals and optional scientific e-notation" semantic="position in sample space in y dimension" unit="nanometer" min="0 m" max="%.4e m" /><field identifier="ImageNumber-0-0" syntax="integer" semantic="frame number" unit="frame" min="0 fr" /><field identifier="Amplitude-0-0" syntax="floating point with . for decimals and optional scientific e-notation" semantic="emission strength" unit="A/D count" /><field identifier="PSFWidth-0-0" syntax="floating point with . for decimals and optional scientific e-notation" semantic="PSF FWHM in x dimension" unit="nanometer" /><field identifier="PSFWidth-1-0" syntax="floating point with . for decimals and optional scientific e-notation" semantic="PSF FWHM in y dimension" unit="nanometer" /><field identifier="FitResidues-0-0" syntax="floating point with . for decimals and optional scientific e-notation" semantic="fit residue chi square value" unit="dimensionless" /><field identifier="LocalBackground-0-0" syntax="floating point with . for decimals and optional scientific e-notation" semantic="local background" unit="A/D count" /></localizations>\r\n', max(rstorm(:,1))*1e-9, max(rstorm(:,2))*1e-9);
            fclose(fid);
            
            dlmwrite(fout, rstorm, '-append', 'delimiter', '\t', 'precision', 6);
            if self.verbose
                fprintf('Created %s\n', fout);
            end
        end
        
        function convertV4290(self, img_dimension)
        % function convertV4290(img_dimension)
        %
        % If you used a version of Insight3 > v4.29.0 to make the molecule 
        % list then you can use this function to transform the molecule 
        % list to an Insight3 version that is <= v4.29.0
        %
        % Also, if you used a version of Insight3 <= v4.29.0 to make the molecule 
        % list then you can use this function to transform the molecule 
        % list to an Insight3 version that is > v4.29.0
        %
        % After v4.29.0 the Insight3 canvas was rotated 90 deg cw and flipped 
        % vertically so that a DAX file shown in Insight3 looks the same as 
        % a DAX file opened in ImageJ.
        %
        % Inputs
        % ------
        % img_dimension : integer
        %   the dimension of the original image file in pixels, eg. 256 for
        %   a 256 x 256 image
        %
            temp = self.getColumn('x');
            self.data(:, self.getColumnIndex('x')) = img_dimension + 1 - self.data(:, self.getColumnIndex('y'));
            self.data(:, self.getColumnIndex('y')) = img_dimension + 1 - temp;
            temp = self.getColumn('xc');
            self.data(:, self.getColumnIndex('xc')) = img_dimension + 1 - self.data(:, self.getColumnIndex('yc'));
            self.data(:, self.getColumnIndex('yc')) = img_dimension + 1 - temp;
        end

        function show(self)
        % function show()
        %
        % Draw the drift-corrected localizations in a 2D histogram
        %
            figure;
            
            ix = self.getColumnIndex('xc');
            iy = self.getColumnIndex('yc');
            
            xmax = max(self.data(:,ix));
            ymax = max(self.data(:,iy));
            
            width = 1;
            while width < xmax
                width = width * 2;
            end
            
            height = 1;
            while height < ymax
                height = height * 2;
            end
            
            pixels = zeros(height, width);
            
            for i=1:self.numMolecules
                x = max(1, min(width, round(self.data(i,ix))));
                y = max(1, min(height, round(self.data(i,iy))));
                pixels(y,x) = pixels(y,x) + 1;
            end
            
            clip = min(max(pixels(:)), 100); % set a range limit in case some pixels have significantly more localizations than the other pixels
            imagesc(pixels, [0 clip]);
            title(strrep(strrep(self.getFilename(), '\', '\\'), '_', '\_')); % set the figure title
            axis('image'); % force the aspect ratio of the plot to be 1
            colormap(gray);
        end
        
        function appendData(self, data)
        % function appendData(data)
        %
        % Append the values in 'data' to the molecule list
        %
        % Inputs
        % ------
        % data : 2D matrix
        %   - can be a Nx2 matrix, where N is the number of molecules and the 2 columns are x, y values
        %   - can be a Nx3 matrix, where N is the number of molecules and the 3 columns are x, y, z values
        %   - can be a Nx18 matrix, where data contains all information
        %  The columns that you don't specify in 'data' get set to a default value
        %
            [numMol, cols] = size(data);
            self.numMolecules = self.numMolecules + numMol;
            
            % ensure that the number of columns in 'data' makes sense
            if ~ismember(cols , [2 3 18])
                error('ERROR Insight3.appendData() :: The data can only have 2, 3 or 18 columns. Got %d columns.\nEither pass only x,y values, x,y,z values, or an entire molecule list', cols);
            end
    
            if strcmp(self.getFormat(), 'txt')           
                if cols == 18
                    self.numFrames = max(self.numFrames, max(data(:,13)));
                    self.data = [self.data; data];
                else
                    self.numFrames = 1;
                    tempdata = zeros(numMol,18);
                    tempdata(:,1) = self.defaultChannel;
                    tempdata(:,2) = data(:,1);
                    tempdata(:,3) = data(:,2);
                    tempdata(:,4) = data(:,1);
                    tempdata(:,5) = data(:,2);
                    tempdata(:,6) = self.defaultHeight;
                    tempdata(:,7) = self.defaultArea;
                    tempdata(:,8) = self.defaultWidth;
                    tempdata(:,9) = self.defaultPhi;
                    tempdata(:,10) = self.defaultAspect;
                    tempdata(:,11) = self.defaultBackground;
                    tempdata(:,12) = self.defaultIntensity;
                    tempdata(:,13) = self.defaultFrame;
                    tempdata(:,14) = self.defaultTrackLength;
                    tempdata(:,15) = self.defaultLink;
                    tempdata(:,16) = self.defaultValid;
                    if cols == 2                        
                        tempdata(:,17) = self.defaultZ;
                        tempdata(:,18) = self.defaultZc;
                    else
                        tempdata(:,17) = data(:,3); 
                        tempdata(:,18) = data(:,3);
                    end
                    self.data = [self.data; tempdata];
                end
            else
                if cols == 18
                    self.numFrames = max(self.numFrames, max(data(:,14)));
                    self.data = [self.data; data];
                else
                    self.numFrames = 1;
                    tempdata = zeros(numMol,18);
                    tempdata(:,1) = data(:,1);
                    tempdata(:,2) = data(:,2);
                    tempdata(:,3) = data(:,1);
                    tempdata(:,4) = data(:,2);
                    tempdata(:,5) = self.defaultHeight;
                    tempdata(:,6) = self.defaultArea;
                    tempdata(:,7) = self.defaultWidth;
                    tempdata(:,8) = self.defaultPhi;
                    tempdata(:,9) = self.defaultAspect;
                    tempdata(:,10) = self.defaultBackground;
                    tempdata(:,11) = self.defaultIntensity;
                    tempdata(:,12) = self.defaultChannel;
                    tempdata(:,13) = self.defaultFitIter;   
                    tempdata(:,14) = self.defaultFrame;
                    tempdata(:,15) = self.defaultTrackLength;
                    tempdata(:,16) = self.defaultLink;
                    if cols == 2                        
                        tempdata(:,17) = self.defaultZ;
                        tempdata(:,18) = self.defaultZc;
                    else
                        tempdata(:,17) = data(:,3); 
                        tempdata(:,18) = data(:,3);
                    end
                    self.data = [self.data; tempdata];
                end
            end            
        end
        
        function sort(self, factor)
        % function sort(factor)
        %
        % Sort the molecule list by the integer representation of xc then by yc.
        %
        % Multiplies the xc and yc values by 'factor', rounds the result
        % down to the nearest integer and then sorts the molecule list by
        % the integer representation of xc then by yc.
        %
        % Inputs
        % ------
        % factor (optional) : double
        %   the factor to multiply the xc and yc values by before sorting
        %   default value is 1.0
        %
            if nargin == 1
                factor = 1.0;
            end
            
            % get the index of the xc and yc columns
            ix = self.getColumnIndex('xc');
            iy = self.getColumnIndex('yc');
            
            % multiply xc and yc by the factor value and 
            % convert to integer values rounding down
            temp = floor(self.data(:,[ix iy])*factor);
            
            % sort by xc then yc
            [~,ind] = sortrows(temp(:,[1 2]));
            self.data = self.data(ind,:);
        end
        
        function removeChannels(self, channels)
        % function removeChannel(channels)
        %
        % Deletes the molecules in the specified channel(s)
        %
        % Inputs
        % ------
        % channel : integer or list of integers
        %   The channel(s) to delete from the molecule list.
        %   Examples
        %   channels = 0 -> remove channel 0
        %   channels = [2 6] -> remove channels 0 and 6
        %
            ic = self.getColumnIndex('channel');
            for i = 1:length(channels)
                self.data = self.data(self.data(:,ic) ~= channels(i), :);
            end
            self.numMolecules = size(self.data, 1);
        end

    end
    
    
    methods(Hidden)
        
        % THIS IS THE ORIGINAL WAY FOR WRITING A MOLECULE LIST TO A FILE
        function success = writeSLOWSLOW(self)
        % function success = writeSLOWSLOW(self)
        %
        % Writes the molecule list to the file defined in self.filename
        % Returns true if the molecule list was successfully written
            
            success = false;

            if isempty(self.data)
                if self.verbose
                    fprintf(2, 'ERROR Insight3.write() :: The molecule list is empty. Nothing to write\n');
                end
                return
            end
        
            % check to see if you are going to overwrite a file
            if ~self.forceOverwrite && exist(self.filename, 'file') == 2
                if self.verbose
                    fprintf(2, 'ERROR Insight3.write() :: Insight3 molecule list exists. Not going to overwrite\n\t%s\n', self.filename);
                end
                return
            end
            
            [~,name,ex] = fileparts(self.filename);
            
            if isempty(ex)
                ex = strcat('.', self.getFormat());
                self.filename = strcat(self.filename, ex);
            end

            % write the data to the file
            if strcmp(self.getFormat(), 'txt')
                if self.verbose
                    fprintf('Writing %s%s, this can take a while... ',name,ex);
                end
                fileID = fopen(self.filename, 'w');
                fprintf(fileID, sprintf('Cas%d\tX\tY\tXc\tYc\tHeight\tArea\tWidth\tPhi\tAx\tBG\tI\tFrame\tLength\tLink\tValid\tZ\tZc\n', self.numMolecules));
                dlmwrite(self.filename, self.data, '-append', 'delimiter', '\t', 'precision', 7);
                fclose(fileID);
                if self.verbose
                    fprintf('DONE\n');
                end
            else
                % create the file reference
                fid = fopen(self.filename, 'wb', 'l'); % byte order is little endian 

                % write the file header
                fwrite(fid, 'M425');
                fwrite(fid, self.numFrames, 'int32');
                fwrite(fid, 6, 'int32');
                fwrite(fid, self.numMolecules, 'int32');

                % define the data types for each column in the bin file
                recordType = {'single' 'single' 'single' 'single' 'single' ...
                    'single' 'single' 'single' 'single' 'single' 'single' ...
                    'int32' 'int32' 'int32' 'int32' 'int32' 'single' 'single'};

                % the number of bytes for each record type
                recordLen = [4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4];

                % write the data to the file in a column-by-column manner
                if self.verbose
                    fprintf('Writing %s%s [                  ]', name, ex);
                end
                for i=1:numel(recordType)
                    if self.verbose
                        fprintf([repmat(sprintf('\b'), 1, 19), repmat('.',1,i), blanks(18-i), ']']);
                    end
                    
                    % go to the correct position in the file to write the first value
                    fseek(fid, 16 + sum(recordLen(1:i-1)), 'bof'); % 16 bytes are for the header

                    % write the first value (the value in row=1)
                    fwrite(fid, self.data(1,i), recordType{i});

                    % write the rest of the values in the column by skipping the specified number of bytes
                    fwrite(fid, self.data(2:end,i), recordType{i}, sum(recordLen) - recordLen(i));
                end

                % write 0 numFrames times so that Insight3 can open the file without crashing
                fwrite(fid, zeros(self.numFrames,1), 'int32');
                fclose(fid);
                if self.verbose
                    fprintf(' DONE\n');
                end
            end                    
            success = true;
        end

        % THIS IS VERSION 2 FOR WRITING A MOLECULE LIST TO A FILE
        function writeSLOW(self)
        % function writeSLOW()
        %
        % Writes the molecule list to the file defined in object.filename
        % Returns true if the molecule list was successfully written
            
            if isempty(self.data)
                error('ERROR Insight3.write() :: The molecule list is empty. Nothing to write');
            end
        
            % check to see if you are going to overwrite a file
            if ~self.forceOverwrite && exist(self.filename, 'file') == 2
                error('ERROR Insight3.write() :: Insight3 molecule list exists. Not going to overwrite\n\t%s\nYou can call i3.forceFileOverwrite(true) to always overwite a file', self.filename);
            end
            
            [~,name,ex] = fileparts(self.filename);
            
            if isempty(ex)
                ex = strcat('.', self.getFormat());
                self.filename = strcat(self.filename, ex);
            end

            % write the data to the file
            if strcmp(self.getFormat(), 'txt')
                if self.verbose
                    fprintf('Writing %s%s, this can take a while... ',name,ex);
                end
                fileID = fopen(self.filename, 'w');
                fprintf(fileID, sprintf('Cas%d\tX\tY\tXc\tYc\tHeight\tArea\tWidth\tPhi\tAx\tBG\tI\tFrame\tLength\tLink\tValid\tZ\tZc\n', self.numMolecules));
                dlmwrite(self.filename, self.data, '-append', 'delimiter', '\t', 'precision', 7);
                fclose(fileID);
                if self.verbose
                    fprintf('DONE\n');
                end
            else
                
                if self.verbose
                    fprintf('Writing %s%s [                  ]', name, ex);
                end

                % create the file, byte order is little endian
                fid = fopen(self.filename, 'wb', 'l');
                
                % write the file header
                fwrite(fid, 'M425');
                fwrite(fid, self.numFrames, 'int32');
                fwrite(fid, 6, 'int32');
                fwrite(fid, self.numMolecules, 'int32');
                
                % Write everything as a 'single' data type.
                % The 5 columns that are 'int32' will be overwritten with the proper data type.
                % Creating the file in this way seems to be ~3.5 times
                % faster than the previous method, see writeSLOWSLOW()
                fwrite(fid, transpose(self.data), 'single');
                
                if self.verbose
                    fprintf([repmat(sprintf('\b'), 1, 19), repmat('.',1,3), blanks(15), ']']);
                end

                % overwrite the 5 'int32' columns, columns 12 to 16
                for i=12:16                    
                    % go to the correct position in the file to overwrite the first value
                    fseek(fid, 16 + (i-1)*4, 'bof'); % skip bytes: 16 bytes are for the header and (i-1)*4 are for the previous columns.

                    % write the first value (the value in row=1)
                    fwrite(fid, self.data(1,i), 'int32');

                    % write the rest of the values in the column by skipping the specified number of bytes
                    remain = self.numMolecules - 1; % the number of values left to write
                    chunk_size = 50000; % write this many values at a time
                    m = 2; % we already wrote the 1st value
                    while remain > 0
                        k = min(chunk_size,remain);
                        fwrite(fid, self.data(m:m+k-1,i), 'int32', 68); % skip 68 bytes (= 4*18 - 4) between writes
                        m = m + k;
                        remain = remain - k;
                    end
                    if self.verbose
                        fprintf([repmat(sprintf('\b'), 1, 16-(i-12)*3), repmat('.',1,3), blanks(15-(i-11)*3), ']']);
                    end
                end
                fseek(fid, 0, 'eof');
                fwrite(fid, zeros(self.numFrames,1), 'int32');
                fclose(fid);
                if self.verbose
                    fprintf(' DONE\n');
                end
            end
        end

        % THIS IS THE SLOW WAY FOR READING A MOLECULE LIST FROM A FILE
        function readSLOW(self)
        % function readSLOW()
        %
        % Returns true if reading the molecule list was successful
            
            self.data = [];
        
            % make sure this Insight3 file exists
            if exist(self.filename, 'file') ~= 2
                error('ERROR Insight3.read() :: Insight3 file does not exist\n\t%s',self.filename)
            end
            
            [~,name,ex] = fileparts(self.filename);

            % want to read a .txt Insight3 file type
            if strcmp(self.getFormat(), 'txt')
                fid = fopen(self.filename, 'rt');
                line = fgets(fid);
                fclose(fid);
                if ~strcmp(line(1:3), 'Cas')
                    error('ERROR Insight3.read() :: This is not a valid .txt molecule list file\n\t%s', self.filename);
                end                    
                if self.verbose
                    fprintf('Reading %s%s, this can take a while... ', name,ex);
                end
                self.data = dlmread(self.filename, '\t', 1, 0);
                self.numMolecules = size(self.data, 1);
                self.numFrames = self.data(end,13);
                if self.verbose
                    fprintf('DONE\n');
                end                
            else
                try
                    % byte order is little endian
                    fid = fopen(self.filename, 'r', 'l');  

                    % read the file header
                    fread(fid, 4, '*char');
                    self.numFrames = fread(fid, 1, 'int32');
                    fread(fid, 1, 'int32');
                    self.numMolecules = fread(fid, 1, 'int32');
                    
                    % define the data types for each column in the bin file
                    recordType = {'*single' '*single' '*single' '*single' '*single' ...
                        '*single' '*single' '*single' '*single' '*single' '*single' ...
                        '*int32' '*int32' '*int32' '*int32' '*int32' '*single' '*single'};
                    
                    % each data type is 4 bytes long
                    recordLen = [4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4]; 

                    % read the file column-by-column
                    self.data = zeros(self.numMolecules, numel(recordType));
                    if self.verbose
                        fprintf('Reading %s%s [                  ]', name, ex);
                    end
                    for i=1:numel(recordType)
                        if self.verbose
                            fprintf([repmat(sprintf('\b'), 1, 19), repmat('.',1,i), blanks(18-i), ']']);
                        end

                        % seek to the first field of the first record
                        fseek(fid, 16+sum(recordLen(1:i-1)), 'bof'); % 16 bytes are the header offset

                        % read column with specified format, skipping required number of bytes
                        self.data(:,i) = fread(fid, self.numMolecules, recordType{i}, sum(recordLen)-recordLen(i));
                    end
                    if self.verbose
                        fprintf(' DONE\n');
                    end
                    fclose(fid);
                catch
                    fclose(fid);
                    error('ERROR Insight3.read() :: Failed to read the .bin molecule list\n\t%s', self.filename);
                end
            end
        end 
        
        function heapError(~, e)
            if isempty(strfind(e.message, 'java.lang.OutOfMemoryError: Java heap space'))
                error(e.message)                 
            else
                 msg = ['Hi AFIBer!, the size of your Java Heap Memory is too low to successfully',char(10),...
                        'read/write Insight3 files. Please increase your Java Heap Memory to be > 300 MB.',char(10),char(10),...
                        'Depending on your version of Matlab there are different ways to do this.',char(10),char(10),...
                        'Method 1 (Matlab 2013b)',char(10),...
                        'Go to Home -> Preferences -> MATLAB -> General -> Java Heap Memory',char(10),char(10),...
                        'Method 2 (Matlab 2013a)', char(10),...
                        'Go to Home -> Preferences -> General -> Java Heap Memory',char(10),char(10),...
                        'Method 3 (Matlab 2011a)', char(10),...
                        'Go to File -> Preferences -> General -> Java Heap Memory',char(10),char(10),...
                        'Method 4 (Matlab ?)', char(10),...
                        'google is your friend...',char(10),char(10),...
                        'You will have to restart Matlab after making this change.'];
                error(msg);
            end
        end

    end
    
end