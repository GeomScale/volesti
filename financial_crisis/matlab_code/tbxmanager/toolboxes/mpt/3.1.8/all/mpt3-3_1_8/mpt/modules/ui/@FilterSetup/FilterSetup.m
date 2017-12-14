classdef FilterSetup
    % Class representing filter's setup

	properties
		replaces
		removes
		incompatibleWith
		callback
		dependsOn
		hidden
		transient % transient filters are not saved
		parser
		filter_name
	end
	
    methods
        function obj = FilterSetup()
            % Constructor
			
			% TODO: optionally disable removal of a filter
			obj.replaces = containers.Map;
			obj.removes = containers.Map;
			obj.incompatibleWith = containers.Map;
			obj.callback = containers.Map;
			obj.dependsOn = containers.Map;
			obj.parser = inputParser;
			obj.hidden = false;
			obj.transient = false;
			
			% determine the filter name
			s = dbstack;
			obj.filter_name = s(2).name(8:end); % kick out the 'filter_' prefix
		end

		function out = actsOn(obj, which)
			% Returns true if the filter modifies a given method
			
			out = obj.callback.isKey(which);
		end
		
		function obj = addField(obj, varargin)
			% Adds a new field to the filter
		
			obj.parser.addParamValue(varargin{:});
		end
		
		function R = parse(obj, varargin)
			% Parses input data
			
			obj.parser.parse(varargin{:});
			R = obj.parser.Results;
			if length(fieldnames(R))==1
				f = fieldnames(R);
				R = R.(f{1});
			end
		end
		
		function out = validate(obj, varargin)
			% Returns true if input data satisfy all validators
			
			if nargin==2 && ~isstruct(varargin{1})
				% single parameter, we need to convert it back to a
				% structure before parsing
				v.(obj.parser.Parameters{1}) = varargin{1};
				try
					obj.parser.parse(v);
				catch
					fprintf('\nValidation failed for filter "%s":\n', obj.filter_name);
					rethrow(lasterror);
				end
			else
				obj.parser.parse(varargin{:});
			end
			out = true;
		end
		
	end
end
