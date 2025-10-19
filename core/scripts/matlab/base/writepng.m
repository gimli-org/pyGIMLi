function writepng(data, map, filename, varargin)
%WRITEPNG Write a PNG file to disk.
%   WRITEPNG(I,[],FILENAME) writes the grayscale image I
%   to the file specified by the string FILENAME.
%
%   WRITEPNG(RGB,[],FILENAME) writes the truecolor image
%   represented by the M-by-N-by-3 array RGB.
%
%   WRITEPNG(X,MAP,FILENAME) writes the indexed image X with
%   colormap MAP.  The resulting file will contain the equivalent
%   truecolor image.
%
%   WRITEPNG(...,PARAM,VAL,...) sets the specified parameters.
%
%   See also IMREAD, IMWRITE, IMFINFO.

%   Copyright 1984-2002 The MathWorks, Inc. 
%   $Revision: 1.1.1.1 $  $Date: 2005/06/01 10:29:16 $

if (ndims(data) > 3)
    error(sprintf('%d-D data not supported for PNG files', ndims(data)));
end

% Color type values (as in PNG library defs)
PNG_COLOR_TYPE_GRAY = 0;
PNG_COLOR_TYPE_RGB = 2;
PNG_COLOR_TYPE_PALETTE = 3;
PNG_COLOR_TYPE_GRAY_ALPHA = 4;
PNG_COLOR_TYPE_RGB_ALPHA = 6;

% Set default parameters
bitdepth = [];
sigbits = [];
interlace = 'none';
transparency = [];
alpha = [];
background = [];
gamma = [];
chromaticities = [];
xres = [];
yres = [];
resunit = [];
textchunks = cell(0,2);

% Process param/value pairs
propStrings = ['interlacetype  '
    'transparency   '
    'bitdepth       '
    'significantbits'
    'alpha          '
    'background     '
    'gamma          '
    'chromaticities '
    'xresolution    '
    'yresolution    '
    'resolutionunit '
    'title          '
    'author         '
    'description    '
    'copyright      '
    'creationtime   '
    'software       '
    'disclaimer     '
    'warning        '
    'source         '
    'comment        '];

for k = 1:2:length(varargin)
    prop = lower(varargin{k});
    if (~isstr(prop))
        error('Parameter name must be a string');
    end
    idx = strmatch(prop, propStrings);
    if (isempty(idx))
        keyword = varargin{k};
        textItem = varargin{k+1};
        keyword = CheckKeyword(keyword);
        textItem = CheckTextItem(textItem);
        textchunks{end+1,1} = keyword;
        textchunks{end,2} = textItem;
    
    elseif (length(idx) > 1)
        error(sprintf('Ambiguous parameter name "%s"', prop));
        
    else
        prop = deblank(propStrings(idx,:));
        switch prop
            case 'bitdepth'
                bitdepth = varargin{k+1};
                
            case 'significantbits'
                sigbits = varargin{k+1};
            
            case 'interlacetype'
                interlace = varargin{k+1};
                
            case 'transparency'
                transparency = varargin{k+1};
                
            case 'alpha'
                alpha = varargin{k+1};
                
            case 'background'
                background = varargin{k+1};
                
            case 'gamma'
                gamma = varargin{k+1};
                
            case 'chromaticities'
                chromaticities = varargin{k+1};
                
            case 'xresolution'
                xres = varargin{k+1};
                
            case 'yresolution'
                yres = varargin{k+1};
                
            case 'resolutionunit'
                resunit = varargin{k+1};
                
            case 'title'
                keyword = 'Title';
                textItem = CheckTextItem(varargin{k+1});
                textchunks{end+1,1} = keyword;
                textchunks{end,2} = textItem;
                
            case 'author'
                keyword = 'Author';
                textItem = CheckTextItem(varargin{k+1});
                textchunks{end+1,1} = keyword;
                textchunks{end,2} = textItem;
                
            case 'description'
                keyword = 'Description';
                textItem = CheckTextItem(varargin{k+1});
                textchunks{end+1,1} = keyword;
                textchunks{end,2} = textItem;
                
            case 'copyright'
                keyword = 'Copyright';
                textItem = CheckTextItem(varargin{k+1});
                textchunks{end+1,1} = keyword;
                textchunks{end,2} = textItem;
                
            case 'creationtime'
                keyword = 'Creation Time';
                if (ischar(varargin{k+1}))
                    textItem = datestr(datenum(varargin{k+1}), 0);
                else
                    textItem = datestr(varargin{k+1}, 0);
                end
                textchunks{end+1,1} = keyword;
                textchunks{end,2} = textItem;
                
            case 'software'
                keyword = 'Software';
                textItem = CheckTextItem(varargin{k+1});
                textchunks{end+1,1} = keyword;
                textchunks{end,2} = textItem;
                
            case 'disclaimer'
                keyword = 'Disclaimer';
                textItem = CheckTextItem(varargin{k+1});
                textchunks{end+1,1} = keyword;
                textchunks{end,2} = textItem;
                
            case 'warning'
                keyword = 'Warning';
                textItem = CheckTextItem(varargin{k+1});
                textchunks{end+1,1} = keyword;
                textchunks{end,2} = textItem;
                
            case 'source'
                keyword = 'Source';
                textItem = CheckTextItem(varargin{k+1});
                textchunks{end+1,1} = keyword;
                textchunks{end,2} = textItem;
                
            case 'comment'
                keyword = 'Comment';
                textItem = CheckTextItem(varargin{k+1});
                textchunks{end+1,1} = keyword;
                textchunks{end,2} = textItem;
                
        end
    end
    
end

if ((ndims(data) > 3) | (~ismember(size(data,3), [1 3])))
    error('Invalid input image.');
end

if (~ismember({class(data)}, {'double', 'logical', 'uint8', 'uint16'}))
    error('Unsupported input data class.');
end

if (~isempty(alpha) & ((size(alpha,1) ~= size(data,1)) | ...
                    (size(alpha,2) ~= size(data,2))))
    error('ALPHA must have the same number of rows and columns as the image data.');
end

%
% Identify color type
%
isTruecolor = (size(data,3) == 3);
paletteUsed = ~isempty(map) & ~isTruecolor;
colorUsed = paletteUsed | isTruecolor;
alphaUsed = ~isempty(alpha);
colortype = paletteUsed + 2*colorUsed + 4*alphaUsed;
if (colortype == 7)
    error('Cannot specify alpha channel with an indexed image.');
end

%
% Set default bitdepth if not specified
%
if (isempty(bitdepth))
    switch class(data)
        case 'logical'
            bitdepth = 1;

        case {'uint8', 'double'}
            bitdepth = 8;

        case 'uint16'
            bitdepth = 16;
    end
end

%
% Validate bitdepth
%
switch colortype
    case PNG_COLOR_TYPE_GRAY
        if (~ismember(bitdepth, [1 2 4 8 16]))
            error('Invalid bitdepth for grayscale image; must be 1, 2, 4, 8, or 16.');
        end
        
    case PNG_COLOR_TYPE_RGB
        if (~ismember(bitdepth, [8 16]))
            error('Invalid bitdepth for RGB image; must be 8 or 16.');
        end
        
    case PNG_COLOR_TYPE_PALETTE
        if (~ismember(bitdepth, [1 2 4 8]))
            error('Invalid bitdepth for indexed image; must be 1, 2, 4, or 8.');
        end
        
    case PNG_COLOR_TYPE_GRAY_ALPHA
        if (~ismember(bitdepth, [8 16]))
            error('Invalid bitdepth for grayscale image with alpha; must be 8 or 16.');
        end
        
    case PNG_COLOR_TYPE_RGB_ALPHA
        if (~ismember(bitdepth, [8 16]))
            error('Invalid bitdepth for RGB image with alpha; must be 8 or 16.');
        end
end

%
% Scale image if necessary to match requested bitdepth
%
switch class(data)
    case 'double'
        if (colortype == PNG_COLOR_TYPE_PALETTE)
            data = data - 1;
            data = uint8(data);
        
        else
            % Grayscale or RGB; clamp data to [0,1] dynamic range before
            % scaling, rounding, and casting.
            data = max(min(data,1),0);
            switch bitdepth
                case 8
                    data = uint8(255*data + 0.5);
                    
                case 16
                    data = uint16(65535*data + 0.5);
                    
                case 4
                    data = uint8(15*data + 0.5);
                    
                case 2
                    data = uint8(3*data + 0.5);
                    
                case 1
                    data = uint8(data ~= 0);
            end
        end
        
    case 'uint8'
        if (colortype == PNG_COLOR_TYPE_PALETTE)
            % Nothing to do
            
        else
            switch bitdepth
                case 16
                    data = uint16(data);
                    data = bitor(bitshift(data,8),data);
                    
                case 8
                    % Nothing to do
                    
                case 4
                    data = bitshift(data,-4);
                    
                case 2
                    data = bitshift(data,-6);
                    
                case 1
                    % Nothing to do
            end
        end
        
    case 'uint16'
        if (colortype == PNG_COLOR_TYPE_PALETTE)
            error('PNG does not allow 16-bit indexed images.');

        else
            switch bitdepth
                case 16
                    % Nothing to do
                    
                case 8
                    data = uint8(bitshift(data,-8));
                    
                case 4
                    data = uint8(bitshift(data,-12));
                    
                case 2
                    data = uint8(bitshift(data,-14));
                    
                case 1
                    data = uint8(data ~= 0);
            end
        end
end

if (ismember(colortype, [PNG_COLOR_TYPE_GRAY_ALPHA, ...
                        PNG_COLOR_TYPE_RGB_ALPHA]))
    %
    % Scale alpha data if necessary to match data class
    %
    switch bitdepth
        case 8
            switch class(alpha)
                case 'double'
                    alpha = max(min(alpha,1),0);
                    alpha = uint8(255 * alpha + 0.5);
                    
                case 'uint16'
                    alpha = uint8(bitshift(alpha, -8));
                    
                case 'uint8'
                    % nothing to do
                    
                otherwise
                    error('Invalid class for alpha');
            end
            
        case 16
            switch class(alpha)
                case 'double'
                    alpha = max(min(alpha,1),0);
                    alpha = uint16(65535 * alpha + 0.5);
                    
                case 'uint16'
                    % nothing to do
                    
                case 'uint8'
                    alpha = uint16(alpha);
                    alpha = bitor(bitshift(alpha, 8), alpha);
                    
                otherwise
                    error('Invalid class for alpha');
            end
    end
end

% Be friendly about specifying resolutions
if (~isempty(xres) & isempty(yres))
    yres = xres;

elseif (~isempty(yres) & isempty(xres))
    xres = yres;
end

if (~isempty(xres) & isempty(resunit))
    resunit = 'unknown';
end

if (isempty(xres) & isempty(yres) & ~isempty(resunit))
    error('X and Y resolutions required when specifying resolution unit.');
end
        
png('write', data, map, filename, colortype, bitdepth, ...
                sigbits, alpha, interlace, ...
                transparency, background, gamma, ...
                chromaticities, xres, yres, ... 
                resunit, textchunks);


function out = CheckKeyword(in)
%CheckKeyword
%   out = CheckKeyWord(in) checks the validity of the input text chunk keyword.

if (isempty(in))
    error('Text chunk keywords must not be empty')
end
if ((in(1) == 32) | (in(end) == 32))
    error('PNG does not allow leading or trailing spaces in text chunk keywords.');
end
if (prod(size(in)) > 80)
    error('Keyword too long; PNG spec limits keyword size to 80 characters.');
end
if (any(~ismember(in,[32:126 161:255])))
    error('Nonprintable characters found in text chunk keyword.');
end

out = in;


function out = CheckTextItem(in)
%CheckTextItem
%   out = CheckTextItem(in) strips out control characters from text; PNG spec
%   discourages them.  It also replaces [13 10] by 10; then it replaces 13 
%   by 10.  The PNG spec says newlines must be represented by a single 10.

if (~ischar(in))
    error('Text chunk must be a string.');
end

out = in;
out = strrep(out, char([13 10]), char(10));
out = strrep(out, char(13), char(10));
badChars = find((out < 32) & (out ~= 10));
if (~isempty(badChars))
    warning('Stripping control characters from text chunk.');
    out(badChars) = [];
end
