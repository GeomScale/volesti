function ts = validate_color(s)
%
% checks the validity of given color
% If the color is supplied as string, it checks whethere the name matches
% the list of available colors. If provided numerically, it must be 3D
% vector defining RGB colors according to Matlab syntax.
%

narginchk(1, 1);

% prepare output
ts = true;

% check RGB values
if isnumeric(s) || islogical(s)
    validate_realvector(s);
    if numel(s)~=3
        error('RGB vector must consist of only 3 elements.');
    end
    if any(s>1) || any(s<0)
       error('The values of RGB vector must lie inside [0, 1] interval.')
    end
    
    return;
end

% assumes that argument is char
if ~ischar(s)
    error('Color name must be given as a string or numerically as RGB vector.');
end  

% list of supported colors
list = {...
    'b';
    'g';
    'r';
    'c';
    'm';
    'y';
    'k';
    'w';
    'aliceblue';
    'antiquewhite';
    'aqua';
    'aquamarine';
    'azure';
    'beige';
    'bisque';
    'black';
    'blanchedalmond';
    'blue';
    'blueviolet';
    'brown';
    'burlywood';
    'cadetblue';
    'chartreuse';
    'chocolate';
    'coral';
    'cornflowerblue';
    'cornsilk';
    'cyan';
    'darkblue';
    'darkcyan';
    'darkgoldenrod';
    'darkgray';
    'darkgreen';
    'darkkhaki';
    'darkmagenta';
    'darkolivegreen';
    'darkorange';
    'darkorchid';
    'darkred';
    'darksalmon';
    'darkseagreen';
    'darkslateblue';
    'darkslategray';
    'darkturquoise';
    'darkviolet';
    'deeppink';
    'deepskyblue';
    'dimgray';
    'dodgerblue';
    'firebrick';
    'floralwhite';
    'forestgreen';
    'fuschia';
    'gainsboro';
    'ghostwhite';
    'gold';
    'goldenrod';
    'gray';
    'green';
    'greenyellow';
    'honeydew';
    'hotpink';
    'indianred';
    'ivory';
    'khaki';
    'lavender';
    'lavenderblush';
    'lawngreen';
    'lemonchiffon';
    'lightblue';
    'lightcoral';
    'lightcyan';
    'lightgoldenrod';
    'lightgoldenrodyellow';
    'lightgray';
    'lightgreen';
    'lightpink';
    'lightsalmon';
    'lightseagreen';
    'lightskyblue';
    'lightslateblue';
    'lightslategray';
    'lightsteelblue';
    'lightyellow';
    'lime';
    'limegreen';
    'linen';
    'magenta';
    'maroon';
    'mediumaquamarine';
    'mediumblue';
    'mediumorchid';
    'mediumpurple';
    'mediumseagreen';
    'mediumslateblue';
    'mediumspringgreen';
    'mediumturquoise';
    'mediumvioletred';
    'midnightblue';
    'mintcream';
    'mistyrose';
    'moccasin';
    'navajowhite';
    'navy';
    'oldlace';
    'olive';
    'olivedrab';
    'orange';
    'orangered';
    'orchid';
    'palegoldenrod';
    'palegreen';
    'paleturquoise';
    'palevioletred';
    'papayawhip';
    'peachpuff';
    'peru';
    'pink';
    'plum';
    'powderblue';
    'purple';
    'red';
    'rosybrown';
    'royalblue';
    'saddlebrown';
    'salmon';
    'sandybrown';
    'seagreen';
    'seashell';
    'sienna';
    'silver';
    'skyblue';
    'slateblue';
    'slategray';
    'snow';
    'springgreen';
    'steelblue';
    'tan';
    'teal';
    'thistle';
    'tomato';
    'turquoise';
    'violet';
    'violetred';
    'wheat';
    'white';
    'whitesmoke';
    'yellow';
    'yellowgreen'};

% check if the string is inside the list
index = find(strcmpi(s,list), 1);

if isempty(index)
    error('Given color name is not in the list of supported colors.');
end
