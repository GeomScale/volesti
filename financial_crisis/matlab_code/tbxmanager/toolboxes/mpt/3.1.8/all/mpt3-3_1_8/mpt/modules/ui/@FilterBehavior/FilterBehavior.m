classdef FilterBehavior < MPTUIHandle
    % Class representing filters
    %
    % Any class derived from FilterBahavior will inherit the addFilter()
    % and applyFilters() methods
	
    methods
        function obj = FilterBehavior()
            % Constructor
			obj.Internal.filter_map = containers.Map;
		end

		function obj = with(obj, filter)
			% Shortcut call to addFilter()
			
			obj.addFilter(filter);
		end
		
		function obj = without(obj, filter)
			% Shortcut call to removeFilter()
			
			obj.removeFilter(filter);
		end

		function keys = listFilters(obj, varargin)
			% Lists all filters associated with this object
			
			% TODO: implement recursion
			ip = inputParser;
			ip.addParamValue('actsOn', '', @ischar);
			ip.addParamValue('incompatibleWith', '', @ischar);
			ip.addParamValue('disabled', true, @islogical);
			ip.addParamValue('enabled', true, @islogical);
			ip.addParamValue('recursive', false, @islogical);
			ip.parse(varargin{:});
			options = ip.Results;

			keys = obj.Internal.filter_map.keys;
			
			if ~isempty(options.actsOn)
				keep = false(1, length(keys));
				for i = 1:length(keys)
					keep(i) = obj.getFilterSetup(keys{i}).actsOn(options.actsOn);
				end
				keys = keys(keep);
			end
			
			enabled = false(1, length(keys));
			for i = 1:length(keys)
				enabled(i) = obj.isFilterEnabled(keys{i});
			end
			if options.disabled && options.enabled
				% keep all filters
				keep = true(1, length(keys));
			elseif options.disabled && ~options.enabled
				% only keep disabled filters
				keep = enabled==0;
			elseif ~options.disabled && options.enabled
				% only keep enabled filters
				keep = enabled==1;
			else
				% keep no filters
				keep = false(1, length(keys));
			end
			keys = keys(keep);
			
			if nargout==0
				for i = 1:length(keys)
					if obj.isFilterEnabled(keys{i})
						status = 'enabled';
					else
						status = 'disabled';
					end
					fprintf('%s (%s)\n', keys{i}, status);
				end
				clear keys
			end
		end
		
		function out = hasFilter(obj, name, varargin)
			% Checks whether the object contains filter "name"
			
			% TODO: return a logical array if 'name' is a cell array

			out = any(obj.Internal.filter_map.isKey(name));
			if ~out && ~isempty(varargin)
				% Should we dive recursively? Note that parsing inputs
				% before determining the status of "out" would be an
				% overkill.
				ip = inputParser;
				ip.addParamValue('recursive', false, @islogical);
				ip.parse(varargin{:});
				options = ip.Results;

				if options.recursive
					% dive recursively
					p = properties(obj);
					for i = 1:length(p)
						v = obj.(p{i});
						for j = 1:numel(v)
							if isa(v(j), 'FilterBehavior')
								out = out || v(j).hasFilter(name, 'recursive', true);
							end
						end
					end
				end
			end
		end

	end

	methods(Hidden)
		% private APIs

		function obj = disableFilter(obj, names)
			% Deactivate given filters

			if ~iscell(names), names = {names}; end
			for i = 1:length(names)
				if obj.hasFilter(names{i})
					f = obj.Internal.filter_map(names{i});
					f.active = false;
					obj.Internal.filter_map(names{i}) = f;
				end
			end
		end

		function obj = enableFilter(obj, names)
			% Activate given filters

			% TODO: when enabling a filter acting on 'generate' we should
			% check whether there isn't an another active filter
			if ~iscell(names), names = {names}; end
			for i = 1:length(names)
				if obj.hasFilter(names{i})
					f = obj.Internal.filter_map(names{i});
					f.active = true;
					obj.Internal.filter_map(names{i}) = f;
				end
			end
		end

		function out = isFilterEnabled(obj, name)
			% Returns activity status of the filter

			out = obj.hasFilter(name) && obj.Internal.filter_map(name).active;
		end
		
		function out = applyFilters(obj, from, concat_type)
			% Applies all filters
			%
			% concat_type can be either '+', 'map', or '[]'
			
			if nargin<3
				concat_type = '+';
			end
			
			filters = obj.listFilters('actsOn', from);
			if isempty(filters)
				out = obj.filterDefaultOutput(from);
			elseif isequal(concat_type, 'map')
				out = containers.Map;
				for i = 1:length(filters)
					out(filters{i}) = obj.callFilterFrom(filters{i}, from, obj.(filters{i}));
				end
			elseif isequal(concat_type, '[]')
				out = obj.callFilterFrom(filters{1}, from, obj.(filters{1}));
				for i = 2:length(filters)
					out = [out, obj.callFilterFrom(filters{i}, from, obj.(filters{i}))];
				end
			else
				out = obj.callFilterFrom(filters{1}, from, obj.(filters{1}));
				for i = 2:length(filters)
					out = out + obj.callFilterFrom(filters{i}, from, obj.(filters{i}));
				end
			end
		end

		function obj = addFilter(obj, filter)
			% Add a filter
			
			% TODO: support multiple filters of the same kind
			if obj.hasFilter(filter)
				warning('Filter "%s" is already present.', filter);
				return
			end
			
			% mark the object as modified
			obj.markAsModified();
			
			% call the filter to determine the setup
			setup = feval(obj.filterFunction(filter), obj);
			if ~isa(setup, 'FilterSetup')
				error('Filter "%s" must return an instance of the FilterSetup class.', filter);
			end
			data = setup.parse();
			
			% only a single filter acting on the 'generate' action can be
			% assigned
			%
			% TODO: this should be part of compatibility checking
			% TODO: we should also print the name of offending filter(s)
			if setup.callback.isKey('generate') && ...
					~isempty(obj.listFilters('actsOn', 'generate', 'disabled', false))
				error('Only a single filter acting on the "generate" action can be used.');
			end
			
			% check dependencies
			keys = setup.dependsOn.keys;
			for i = 1:length(keys)
				recursion = isequal(setup.dependsOn(keys{i}), 'recursive');
				if ~obj.hasFilter(keys{i}, 'recursive', recursion)
					error('Filter "%s" depends on "%s", which is not present.', ...
						filter, keys{i});
				end
			end
			
			% check compatibility
			keys = setup.incompatibleWith.keys;
			for i = 1:length(keys)
				recursion = isequal(setup.incompatibleWith(keys{i}), 'recursive');
				if obj.hasFilter(keys{i}, 'recursive', recursion)
					error('Filter "%s", is incompatible with "%s".', filter, keys{i});
				end
			end
			
			% removes any filters replaced by the current on
			obj.removeFilter(setup.removes.keys);
			
			% deactivate filters if necessary
			obj.disableFilter(setup.replaces.keys);
			
			% create a new dynamic property and assign data
			H = addprop(obj, filter);
			obj.(filter) = data;
			if setup.callback.isKey('set')
				% use provided setter method
				H.SetMethod = setup.callback('set');
			end
			if setup.callback.isKey('get')
				% use provided getter method
				H.GetMethod = setup.callback('get');
			end
			
			% listen to changes of the filter's value
			H.SetObservable = true;
			obj.addlistener(filter, 'PostSet', @obj.markAsModified);
			
			% hide the filter if requested
			H.Hidden = setup.hidden;
			
			% structure in which we keep the filter settings
			f.handle = H;
			f.setup = setup;
			f.active = true;
			obj.Internal.filter_map(filter) = f;
			
			try
				% immediately execute the filter if needed
				obj.callFilterFrom(filter, 'addFilter');
			catch
				% filter failed, remove it
				obj.removeFilter(filter);
				rethrow(lasterror);
			end
			
		end
		
		function obj = removeFilter(obj, filter)
			% Remove a filter
			
			% mark the object as modified
			obj.markAsModified();

			if iscell(filter)
				% support removing multiple filters
				for i = 1:length(filter)
					obj.removeFilter(filter{i});
				end
				
			elseif obj.hasFilter(filter)
				% remove this filter
				
				obj.callFilterFrom(filter, 'removeFilter');
				delete(obj.Internal.filter_map(filter).handle);
				obj.Internal.filter_map.remove(filter);
			end
		end

		function new = copyFiltersFrom(new, obj)
			% removes all filters from 'new' and replaces them by filters
			% defined for 'obj'
			
			new.Internal.filter_map = containers.Map;
			filters = obj.Internal.filter_map.keys;
			% TODO: add filters in the order in which they were added
			for i = 1:length(filters)
				% do not copy transient filters
				filter = obj.Internal.filter_map(filters{i});
				if filter.setup.transient
					% do not copy transient filters
					continue
				end
				new.addFilter(filters{i});
				% re-assign values
				new.(filters{i}) = obj.(filters{i});
			end
		end
		
		function obj = saveAllFilters(obj)
			% saves arguments of all filters to
			% obj.Internal.save.filters
			
			% store arguments of filters
			filters = obj.Internal.filter_map.keys;
			arguments = [];
			for i = 1:length(filters)
				filter = obj.Internal.filter_map(filters{i});
				if filter.setup.transient
					% do not save transient filters
					continue
				end
				arguments.(filters{i}) = obj.(filters{i});
			end
			obj.Internal.save.filters = arguments;
			% remove all filters from the object before saving (note that
			% we need to work with a copy of the object in saveobj).
			% removing filters is necessary, or otherwise "clear classes"
			% would no longer work correctly after subsequent load.
			obj.removeFilter(obj.Internal.filter_map.keys);
		end
		
		function obj = loadSavedFilters(obj)
			% re-initializes filters using arguments stored in
			% obj.Internal.save.filters

			if isfield(obj.Internal, 'save') && ...
					isfield(obj.Internal.save, 'filters') && ...
					isstruct(obj.Internal.save.filters)
				arguments = obj.Internal.save.filters;
				filters = fieldnames(arguments);
				% clear filters
				obj.Internal.filter_map = containers.Map;
				% now add the filters back and re-assign arguments
				%
				% TODO: we should respect the order in which filters were added
				for i = 1:length(filters)
					obj.addFilter(filters{i});
					obj.(filters{i}) = arguments.(filters{i});
				end
			end
		end
	end
	
	methods(Access=protected)

		function out = callFilterFrom(obj, name, from, varargin)
			% Executes a given filter

			filters = obj.Internal.filter_map;
			if filters(name).active && ...
					filters(name).setup.actsOn(from)
				setup = filters(name).setup;
				setup.validate(varargin{:});
				callback = setup.callback(from);
				out = callback(obj, varargin{:});
			else
				out = obj.filterDefaultOutput(from);
			end
		end

		function out = isFilterCompatible(obj, name)
			% Checks whether the filter is compatible with an another
			% filter
			
			% TODO: finish, must support recursive call
			out = true;
		end
		
		function out = getFilterSetup(obj, name)
			% Returns setup of a given filter
			
			if obj.hasFilter(name)
				f = obj.Internal.filter_map(name);
				out = f.setup;
			else
				out = [];
			end
		end

		function out = hasFilterOnAction(obj, action, varargin)
			% Checks whether the object has some filters which responds to
			% a given action
			
			ip = inputParser;
			ip.addParamValue('recursive', false, @islogical);
			ip.addParamValue('includeDisabled', false, @islogical);
			ip.parse(varargin{:});
			options = ip.Results;
			
			out = false;
			filters = obj.listFilters();
			for i = 1:length(filters)
				setup = obj.getFilterSetup(filters{i});
				hasFilter = setup.callback.isKey(action);
				if ~options.includeDisabled
					% only consider enabled filters
					hasFilter = hasFilter && obj.isFilterEnabled(filters{i});
				end
				out = out || hasFilter;
				if ~out && options.recursive
					% dive recursively
					p = properties(obj);
					for i = 1:length(p)
						v = obj.(p{i});
						for j = 1:numel(v)
							if isa(v(j), 'FilterBehavior')
								out = out || v(j).hasFilterOnAction(action, 'recursive', true);
							end
						end
					end
				end
				
			end
		end

    end
    
    methods(Static, Hidden)
		
		function out = filterFunction(name)
			out = ['filter_' name];
		end
       
        function [out, exit] = filterDefaultOutput(from, desired)

			switch from
				case {'constraints', 'instantiate', 'uninstantiate', ...
						'generate', 'construct', 'getVariables', ...
						'setup', 'addFilter', 'removeFilter' };
					out = [];
				case 'objective'
					out = 0;
				otherwise
					error('Unrecognized caller "%s".', from);
			end

			exit = nargin==2 && ~isequal(from, desired);
		end
    end

end
