function hysdel(filename,simname,options)
% hysdel(filename,simname,options)
% Compiles the HYSDEL list <filename.hys>, 
% generates the M-file <filename.m>
% 
% INPUT: 
% filename: the hysdel source 
% simname : if not empty, generates the HYSDEL simulator (see manual)
% options : a string of space separated command line switches to append to 
%           HYSDEL call (see HYSDEL manual, it includes: -p -a -5 
%           --no-symbol-table --no-row-info --no-params-checks 
%           --matlab-symbolic -v[0-3])
%
% OUTPUT:
% filename.m and simname.m on disk 
%
% (C) 2000--2002 F.D. Torrisi,
% Automatic Control Laboratory, ETH Zentrum, CH-8092 Zurich, Switzerland
% torrisi@aut.ee.ethz.ch
%
% see license.txt for the terms and conditions.

% Removes extension

global MPTOPTIONS
if isempty(MPTOPTIONS)
	MPTOPTIONS=mptopt;
end

hysdel_cmd = MPTOPTIONS.modules.hysdel2.hysdel_executable;

if (nargin < 3)
	options = [];
	if (nargin < 2)
		simname = [];
		if (nargin < 1)
        		% print version and exit
			cmd_line = ' -V';
			eval(['!' hysdel_cmd cmd_line ]); 
			return 
		end
	end		
end
	

j=findstr(filename,'.');
if ~isempty(j),
	if strcmpi(filename(j(end):end),'.hys')
		filename(j(end):end)=[];
	end
end

j=findstr(simname,'.');
if ~isempty(j),
	if strcmpi(simname(j(end):end),'.m')
        	simname(j(end):end)=[];
	end
end


cmd_line = [' -i' filename '.hys -m' filename ];
if ~isempty(simname),
	cmd_line = [cmd_line ' -s' simname ' '];
end

% compile
eval(['!' hysdel_cmd cmd_line ' ' options]);
    
% execute HYSDEL output to load the system S 
% eval(filename); % this will not work if your system has symbolic parameters
                  % and therefore has been disabled
% write 'filename.mat'
% eval(['save ' filename]);



