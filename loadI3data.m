% function to check an input value presumed to be related to an Insight3
% localization list and return that list.
%
% Input types & corresponding outputs:
%   Empty or missing
%       Call format: LL = checkI3Input;
%                    LL = checkI3Input( [] );
%       Here the user is asked to select a .bin or .txt localization list
%       The list is loaded and returned as an Insight3 object
%   String indicating a full file path of a localization list in either 
%       .bin or .txt format:
%       Call format: LL = checkI3Input('path\file.bin');
%                    LL = checkI3Input('path\file.txt');
%       Here the Localization list is loaded and returned
%       Error checking of string inputs:
%           If the input is a character string that is lacking an extension 
%               that would otherwise be a valid and existing localization 
%               list .bin or .txt file, then it is loaded
%           If the input is a file name without a path or an extension that
%               is not .bin or .txt, then the user is asked to select a 
%               valid .bin or .txt list
%   A matrix of either 2, 3, or 18 columns
%       In this case the user is asked to save the data matrix as a .bin file
%       then the Insight3 object is returned
%   An Insight3 object 
%       Then the same object is returned

function LL = loadI3data(LLname)

if ~exist('LLname','var') || ...
        isempty(LLname) || ...
        ( ischar(LLname) && (isempty(strfind(LLname,'.bin')) && isempty(strfind(LLname,'.txt'))) )
    % Then the user included something in the LLname field, so check if it
    % is and incorrectly formatted file fpath name
    askToSelect = true;
    if exist('LLname','var') && ~isempty(LLname) && ischar(LLname)
        [fpath,~,ext]=fileparts(LLname);
        if ~isempty(fpath) && isempty(ext) 
            % name includes fpath, but no extension included
            % test to see if only an extension is missing ...
            if exist([LLname '.bin'],'file') 
                % then it's a .bin Loc List
                LLname = [LLname '.bin'];
                askToSelect = false;
            elseif exist([LLname '.txt'],'file') 
                % then it's a .txt Loc List
                LLname = [LLname '.txt'];
                askToSelect = false;
            end
%         else
%             % presume no fpath in the input name or some weird extention 
%             askToSelect = true;
        end
%     else 
%         askToSelect = true;
    end
    if askToSelect
        if exist('LLname','var') && ~isempty(LLname) && ischar(LLname) && exist(LLname,'dir')
            % input is a starting directory to select a .bin file
            [file, fpath] = uigetfile([LLname '/*.bin'],'Select Insight3 Molecule List');
        else
            % then query user to select .bin molecule list file saved by Insight3
            [file, fpath] = uigetfile({'*.bin';'*.txt'},'Select Insight3 Molecule List');
        end
        if ~exist(fullfile(fpath,file),'file')
            error(['Molecule list does not exist, input filename = ' file])
        end
    else
        [fpath,file,ext]=fileparts(LLname);
        file = [file ext];
    end
    LL = Insight3(fullfile(fpath,file));
    
elseif isobject(LLname) &&  isa(LLname,'Insight3')
    % then input is actually a molecule list
    LL = LLname;
    
elseif ischar(LLname) && (~isempty(strfind(LLname,'.bin')) || ~isempty(strfind(LLname,'.txt')))
    % then input is correctly formatted LL file location
%     idx = find(LLname=='\',1,'last');
%     fpath = LLname(1:idx);
%     file = LLname(idx+1:end);
    [fpath,file,ext]=fileparts(LLname);
    file = [file ext];
    if ~exist(fullfile(fpath,file),'file')
        error(['Molecule list does not exist, input filename = ' file])
    end
    LL = Insight3(LLname);
    
elseif isnumeric(LLname) 
    % input is a molecule list or part thereof
    % ask user to save the list and return it as an Insight3 object
    if sum(size(LLname,2)==[2,3,18]) % presume it's going to have the common column assignment
      [file,fpath,~] = uiputfile('*.bin','Save matrix to .bin List');
      LL = Insight3();
      LL.setData( LLname );
      LL.setFilename( fullfile(fpath,file) );
      LL.write();
    else
        error('input Matrix is of unrecognized format, it must have either 2, 3 or 18 columns')
    end
else
    % then the user input something for LLname with an unknown format
    error('input is of unrecognized format')
end