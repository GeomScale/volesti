function c = charToColor(s, color_map, index, colororder)
%
% convert color string to RGB vector using an optional color map
%
% color_map must be either a string or a compatible double array

global MPTOPTIONS
if isempty(MPTOPTIONS)
	MPTOPTIONS=mptopt;
end

narginchk(1, 4);

default_index = 1;
default_colororder = 'fixed';
if nargin<2
	% default color map
	color_map = 'mpt';
	index = default_index;
	colororder = default_colororder;
elseif nargin<3
	index = default_index;
	colororder = default_colororder;
elseif nargin<4
	colororder = default_colororder;
end

validate_indexset(index);
if isnumeric(color_map) && size(color_map, 2)==3
	% numerical colormap provided
	colors = color_map;
elseif ischar(color_map)
	if isequal(lower(color_map), 'mpt')
		% use MPT's defualt color map
		colors = MPTOPTIONS.colormap;
	elseif isequal(lower(color_map), 'lightgray')
		% custom color map for publications
		colors = feval('gray')*0.5+0.4;
	else
		% for 'hsv', 'bone', etc.
		colors = feval(color_map);
	end
else
	error('Wrong type of the color map. Must be either double or char.');
end

if isequal(colororder, 'random')
	% permute the colors randomly
	colors = colors(randperm(size(colors, 1)), :);
end

if isempty(s)
    s = 'x';
end
if ~ischar(s)
    error('Color name must be a string.');
end

% make string lowercase
s = lower(s);

switch s
    case 'b', c = [ 0 0 1 ];
    case 'g', c = [ 0 1 0 ];
    case 'r', c = [ 1 0 0 ];
    case 'c', c = [ 0 1 1 ];
    case 'm', c = [ 1 0 1 ];
    case 'y', c = [ 1 1 0 ];
    case 'k', c = [ 0 0 0 ];
    case 'w', c = [ 1 1 1 ];
    case 'aliceblue', c = [240 	248 	255]/255;
    case 'antiquewhite', c = [250 	235 	215]/255;
    case 'aqua', c = [0 	255 	255]/255;
    case 'aquamarine', c =	[127 	255 	212]/255;
    case 'azure', c = [240 	255 	255]/255;
    case 'beige', c = [245 	245 	220]/255;
    case 'bisque', c = [255 	228 	196]/255;
    case 'black', c = [0 0 	0];
    case 'blanchedalmond', c = [255 	235 	205]/255;
    case 'blue', c = [0 	0 	255]/255;
    case 'blueviolet', c = [138 	43 	226]/255;
    case 'brown',   c = [165 	42 	42]/255;
    case 'burlywood', c = [222 	184 	135]/255;
    case 'cadetblue', c = [95 	158 	160]/255;
    case 'chartreuse', c =  [127 	255 	0]/255;
    case 'chocolate', c = [210 	105 	30]/255;
    case 'coral', c = [255 	127 	80]/255;
    case 'cornflowerblue', c = [100 	149 	237]/255;
    case 'cornsilk', c = [255 	248 	220]/255;
    case 'cyan', c = [0 	255 	255]/255;
    case 'darkblue', c = [0 	0 	139]/255;
    case 'darkcyan', c = [0 	139 	139];
    case 'darkgoldenrod', c = [184 	134 	11]/255;
    case 'darkgray', c = [169 	169 	169]/255;
    case 'darkgreen', c = [0 	100 	0]/255;
    case 'darkkhaki', c = [189 	183 	107]/255;
    case 'darkmagenta', c = [139 	0 	139]/255;
    case 'darkolivegreen', c = [85 	107 	47]/255;
    case 'darkorange', c = [255 	140 	0]/255;
    case 'darkorchid', c = [153 	50 	204]/255;
    case 'darkred', c = [139 	0 	0]/255;
    case 'darksalmon', c = [233 	150 	122]/255;
    case 'darkseagreen', c = [143 	188 	143]/255;
    case 'darkslateblue', c = [72 	61 	139]/255;        
    case 'darkslategray', c = [47 	79 	79]/255;
    case 'darkturquoise', c = [0 	206 	209]/255;
    case 'darkviolet',  c = [148 	0 	211]/255;
    case 'deeppink', c = [255 	20 	147]/255;
    case 'deepskyblue', c = [0 	191 	255]/255;
    case 'dimgray', c = [105 	105 	105]/255;
    case 'dodgerblue', c = [30 	144 	255]/255;
    case 'firebrick', c = [178 	34 	34]/255;
    case 'floralwhite', c = [255 	250 	240]/255;
    case 'forestgreen', c = [34 	139 	34]/255;
    case 'fuschia', c = [255 	0 	255]/255;
    case 'gainsboro', c = [220 	220 	220]/255;
    case 'ghostwhite', c = [255 	250 	250]/255;
    case 'gold', c = [255 	215 	0]/255;
    case 'goldenrod', c = [218 	165 	32]/255;
    case 'gray', c = [128 	128 	128]/255;
    case 'green', c = [0 	128 	0]/255;
    case 'greenyellow', c = [173 	255 	47]/255;
    case 'honeydew', c = [240 	255 	240]/255;
    case 'hotpink', c = [255 	105 	180]/255;
    case 'indianred', c = [205 	92 	92]/255;
    case 'ivory', c = [255 	255 	240]/255;
    case 'khaki', c = [240 	230 	140]/255;
    case 'lavender', c = [230 	230 	250]/255;
    case 'lavenderblush', c = [255 	240 	245]/255;
    case 'lawngreen', c = [124 	252 	0]/255;
    case 'lemonchiffon', c = [255 	250 	205]/255;
    case 'lightblue', c = [173 	216 	230]/255;
    case 'lightcoral', c = [240 	128 	128]/255;
    case 'lightcyan', c = [224 	255 	255]/255;
    case 'lightgoldenrod', c = [238 	221 	130]/255;
    case 'lightgoldenrodyellow', c = [250 	250 	210]/255;
    case 'lightgray', c = [211 	211 211]/255;
    case 'lightgreen', c = [144 	238 	144]/255;
    case 'lightpink', c = [255 	182 	193]/255;
    case 'lightsalmon', c = [255 	160 	122]/255;
    case 'lightseagreen', c = [32 	178 	170]/255;
    case 'lightskyblue', c = [135 	206 	250]/255;
    case 'lightslateblue', c = [132 	112 	255]/255;
    case 'lightslategray', c = [119 	136 	153]/255;
    case 'lightsteelblue', c = [176 	196 	222]/255;
    case 'lightyellow', c = [255 	255 	224]/255;
    case 'lime', c = [0 	255 	0]/255;
    case 'limegreen', c = [50 	205 	50]/255;
    case 'linen', c = [250 	240 	230]/255;
    case 'magenta', c = [255 	0 	255]/255;
    case 'maroon', c = [128 	0 	0]/255;
    case 'mediumaquamarine', c = [102 	205 	170]/255;
    case 'mediumblue', c = [0 	0 	205]/255;
    case 'mediumorchid', c = [186 	85 	211]/255;
    case 'mediumpurple', c = [147 	112 	219]/255;
    case 'mediumseagreen', c = [60 	179 	113]/255;
    case 'mediumslateblue', c = [123 	104 	238]/255;
    case 'mediumspringgreen', c = [0 	250 	154]/255;
    case 'mediumturquoise', c = [72 	209 	204]/255;
    case 'mediumvioletred', c = [199 	21 	133]/255;
    case 'midnightblue', c = [25 	25 	112]/255;
    case 'mintcream', c = [245 	255 	250]/255;
    case 'mistyrose', c = [255 	228 	225]/255;
    case 'moccasin', c = [255 	228 	181]/255;
    case 'navajowhite', c = [255 	222 	173]/255;
    case 'navy', c = [0 	0 	128]/255;
    case 'oldlace', c = [253 	245 	230]/255;
    case 'olive', c = [128 	128 	0]/255;
    case 'olivedrab', c = [107 	142 	35]/255;
    case 'orange', c = [255 	165 	0]/255;
    case 'orangered', c = [255 	69 	0]/255;
    case 'orchid', c = [218 	112 	214]/255;
    case 'palegoldenrod', c = [238 	232 	170]/255;
    case 'palegreen', c = [152 	251 	152]/255;
    case 'paleturquoise', c = [175 	238 	238]/255;
    case 'palevioletred', c = [219 	112 	147]/255;
    case 'papayawhip', c = [255 	239 	213]/255;
    case 'peachpuff', c = [255 	218 	185]/255;
    case 'peru', c = [205 	133 	63]/255;
    case 'pink', c = [255 	192 	203]/255;
    case 'plum', c = [221 	160 	221]/255;
    case 'powderblue', c = [176 	224 	230]/255;
    case 'purple', c = [128 	0 	128]/255;
    case 'red', c = [255 	0 	0]/255;
    case 'rosybrown', c = [188 	143 	143]/255;
    case 'royalblue', c = [65 	105 	225]/255;
    case 'saddlebrown', c = [139 	69 	19]/255;
    case 'salmon', c = [250 	128 	114]/255;
    case 'sandybrown', c = [244 	164 	96]/255;
    case 'seagreen', c = [46 	139 	87]/255;
    case 'seashell', c = [255 	245 	238]/255;
    case 'sienna',  c = [160 	82 	45]/255;
    case 'silver', c = [192 	192 	192]/255;
    case 'skyblue', c = [135 	206 	235]/255;
    case 'slateblue', c = [106 	90 	205]/255;
    case 'slategray', c = [112 	128 	144]/255;
    case 'snow', c = [255 	250 	250]/255;
    case 'springgreen', c = [0 	255 	127]/255;
    case 'steelblue', c = [70 	130 	180]/255;
    case 'tan', c = [210 	180 	140]/255;
    case 'teal', c = [0 	128 	128]/255;
    case 'thistle', c = [216 	191 	216]/255;
    case 'tomato', c = [255 	99 	71]/255;
    case 'turquoise', c = [64 	224 	208]/255;
    case 'violet', c = [238 	130 	238]/255;
    case 'violetred', c = [208 	32 	144]/255;
    case 'wheat', c = [245 	222 	179]/255;
    case 'white', c = [255 	255 	255]/255;
    case 'whitesmoke', c = [245 	245 	245]/255;
    case 'yellow', c = [255 	255 	0]/255;
    case 'yellowgreen', c = [154 	205 	50]/255;
    otherwise,
		% adjust the index to the size of the number of colors
		c = colors(mod(index-1, size(colors, 1))+1, :);
end

end
