% a function using a localization list and conventional pixel size +
% desired S.R. pixel size and will calculate two image masks:
%   mask.conv that is the conventional-mask corresponding to a conventional
%   image
%   mask.SR that is a mask created from a super-resolved image
%
% Inputs
%   LL = Insight3 object containing localization list data
%   params = One of the following input formats are acceptable:
%       1) a 2-element array containing the camera pixel size in nanometers
%           and the super-resolved pixel size in nanometers
%       2) a structure containg fields
%           'original_pixel_size'
%           'localization_precision' 
%          Each field having a single value in units of nanometers
%       3) Output from FindClustersStruct.m
%
% Output
%   mask = structure having the following fields
%       'conv' = binary mask with pixel size equal to either the larger of
%           the two 'param' values or the 'original_pixel_size' field
%       'convArea' = image area in micrometer^2 
%       'SR' = binary mask with pixel size equal to either the smaller of
%           the two 'param' values or the 'localization_precision' field
%       'SRArea' = Super Resolved image area in micrometer^2


function mask = calcMaskAreas( LL, params, showText )

if nargin < 2
    error('Must provide two inputs to calcMaskAreas function')
elseif nargin == 2
    showText = false;
else 
    showText = logical( showText );
end

LL = loadI3data( LL );

reqFields = {'original_pixel_size','localization_precision'};
if ~isa(params,'FindClustersStruct')
    if isnumeric(params)
        if numel(params) == 2
            pix2nm = max(params);
            finalPix = min(params);
            params = struct(...
                reqFields{1}, pix2nm,...
                reqFields{2}, finalPix...
                );
        else
            error('The second input to calcMaskAreas should have been a 2-element array: [ pix2nm, finalPix ]')
        end
    elseif isstruct(params)
        errStr = 'The second input to calcMaskAreas does not have fields required fields: ';
        err = false;
        for fld = 1:length(reqFields)
            if ~isfield(params,reqFields{fld})
                errStr = [errStr reqFields{fld} ','];
                err = true;
            end
        end
        if err
            error(errStr(1:end-1))
        end
    else
        error('The second input to calcMaskAreas is not of FindClustersStruct data type');
    end
end

% Calculate conventional image mask from localizations and get its area in um^2
[mask.conv,mask.convArea] = Locs2Mask(LL.getXYcorr,[],[],showText);
mask.convArea = mask.convArea*(params.original_pixel_size/1000)^2; %pix^2 * (um/pix)^2

% Calculate Super Resolution (S.R.) image mask and get its area in um^2
[mask.SR,mask.SRArea] = Locs2Mask(LL.getXYcorr, ...
                                params.localization_precision, ...
                                params.original_pixel_size,...
                                showText);
mask.SRArea = mask.SRArea*(params.localization_precision/1000)^2; %pix^2 * (um/pix)^2

end % of function