% Input an [x,y] list of localization positions and generate a mask image
% enclosing the area they comprise.  The image can optionally be resized to
% obtain a desired pixel size
%
% Inputs
%   xy - two-column matrix with x and y coordinates of localizations in
%       presumed units are pixels, since this makes the calculation faster
%   finalPix - (Optional) If this value is not included, then the output
%       mask will have pixel size of equivalent units to the input xy
%       values. If a value is included, then the final mask will have a
%       pixel size corresponding to this input value.  
%       For example, finalPix = 20 will give a mask having 20nm pixel size
%   pix2nm - (Optional) conversion from units of pixels to nanometers
%   showText - (optional) binary input to determine if the user should be
%       updated about calculation status
%
% Outputs
%   mask - binary image with 1 corresponding to the presence of a localization
%   area - total area of the mask having localizations, units are equal to
%       those of finalPix if provided, otherwise in units of xy if no 
%       finalPix is specified
%   pixLocs - cell array of dimensions equal to mask where each cell
%       contains the indicies of the input xy coordinates contained within
%       each pixel in the mask

function [mask,area,pixLocs] = Locs2Mask( xy, finalPix, pix2nm, showText )


% parse inputs
if min(size(xy)) ~= 2
    error('Input xy list must be a 2-column matrix')
end
if size(xy,2) > size(xy,1)
    xy = transpose(xy);
end
if ~exist('pix2nm','var') || isempty(pix2nm)
    pix2nm = 1;
end
if ~exist('finalPix','var') || isempty(finalPix)
    finalPix = 1;
end
if ~exist('showText','var') || isempty(showText)
    showText = true;
else
    showText = logical(showText);
end
% set coordinates to correspond with final mask size
xy = xy*pix2nm/finalPix;


%% find the limits where localizations extend through
xyLims = [floor(min(xy(:,2))),ceil(max(xy(:,2)));... %y row1
          floor(min(xy(:,1))),ceil(max(xy(:,1)))];   %x row2
maxdim = max(xyLims(:,2));
% initialize mask
mask = false(maxdim,maxdim);

% round input localizations
xyRd = round(xy); 
% sort for rapid pixel assignment
[xyS,ixy] = sortrows(xyRd,[1 2]); % xyRd(ixy) == xyS
% remove possible zeros
xyS = xyS( xyS(:,1)>0 ,:);
xyS = xyS( xyS(:,2)>0 ,:);
% obtain unique rounded values for pixel assignment
unqXY = unique( xyS,'rows' ); 

%% find localizations contained within each pixel
pixLocs = cell(maxdim,maxdim);
% number of occurrances of unique column X 
XxyS = histc(xyS(:,1),unqXY(:,1)); 
% number of times column X appears in unique values
XunqXY = find(XxyS); 
% number of times column X appears in rounded, sorted values
XxyS = XxyS(XunqXY);
% convert to indicies corresponding to xyS
XxyS = cumsum(XxyS);
XunqXY = [0;XunqXY];
XxyS = [0;XxyS];

idxPrev = 1;
n_xcols = length(XunqXY);
% tmp = []; m=0; tmp2 = {};
tmr = tic;
strCR = '';
%%
for i = 2:n_xcols
    % number of occurrances of row Y pixels in selected column X
    YxyS = histc( xyS(XxyS(i-1)+1:XxyS(i),2), unqXY(XunqXY(i-1)+1:XunqXY(i),2) ); 
    
    for j = 1:length(YxyS)
        idxNext = YxyS(j)+idxPrev-1;
        maskidx = [unqXY(XunqXY(i-1)+j,2),unqXY(XunqXY(i-1)+j,1)];
        pixLocs{ maskidx(1),maskidx(2) } = ...
            ixy(idxPrev:idxNext);
        mask( maskidx(1),maskidx(2) ) = true;
        idxPrev = idxNext+1;
    end
    
    %% update user
    if showText
        strout = sprintf('Identifying pixels containing localizations: %.1f %%%% \n',100*i/n_xcols);
        fprintf([strCR strout]);
        strCR = repmat('\b',1,length(strout)-1);
    end
    
end

if showText
    toc(tmr)
end
        
% calculate the area containing localizations
if nargout > 1
    area = bwarea(mask);% * (finalPix^2);% * (pix2nm^2);
end
end