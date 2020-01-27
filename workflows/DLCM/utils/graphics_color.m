function rgb = graphics_color(name)
%GRAPHICS_COLOR Color selection.
%   RGB = GRAPHICS_COLOR(NAME) selects a suitable color for use in
%   scientific plots and presentation, following best advice on
%   clarity and visibility [1]. NAME is a color, see table below, and
%   RGB is a length 3 vector in [0,1].
%
%   Name            RGB (0--255)      CMYK (%)
%   ----------------------------------------------
%   Black           (0,0,0)           (0,0,0,100)
%   Orange          (230,159,0)       (0,50,100,0)
%   Sky blue        (86,180,233)      (80,0,0,0)
%   Bluish green    (0,158,115)       (97,0,75)
%   Yellow          (240,228,66)      (10,5,90,0)
%   Blue            (0,114,178)       (100,50,0,0)
%   Vermillion      (213,94,0)        (0,80,100,0)
%   Reddish purple  (204,121,167)     (10,70,0,0)
%
%   Reference:
%     [1] B. Wong: "Points of view: Color blindness", Nature
%     Methods 8(6):441 (2011), doi:10.1038/nmeth.1618

% S. Engblom 2017-08-29

switch lower(name)
 case 'black', rgb = [0,0,0];
 case 'orange', rgb = [230,159,0];
 case 'sky blue', rgb = [86,180,233];
 case 'bluish green', rgb = [0,158,115];
 case 'yellow', rgb = [240,228,66];
 case 'blue', rgb = [0,114,178];
 case 'vermillion', rgb = [213,94,0];
 case 'reddish purple', rgb = [204,121,167];
 otherwise, error('Color not found.');
end
rgb = rgb/255;
