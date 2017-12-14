classdef ComponentBehavior < MPTUIHandle
    % Class representing a storage of components
   
    methods(Hidden)
        function obj = ComponentBehavior()
            % Constructor
            
            obj.Internal.component.type = '';
            obj.Internal.component.map = containers.Map;
        end
        
        function out = getType(obj)
            % Returns component's type
            
            out = obj.Internal.component.type;
        end

        function obj = setType(obj, t)
            % Sets component's type
            
            obj.Internal.component.type = t;
        end
        
        
        function out = isType(obj, t)
            % Returns true if the component's type matches the input
            
            out = all(isequal(obj.getType(), t));
        end
        
        function obj = addComponent(obj, name, data)
            % Adds an arbitrary object to the composition under a given
            % name
            
            % check that the name is allowed
            forbidden_names = {'Internal'};
            if ~isempty(strmatch(name, forbidden_names, 'exact'))
                error('Name "%s" is not allowed.', name);
            end
            
            props = properties(obj);
            if ~isempty(strmatch(name, props, 'exact'))
                % update existing property
                obj.(name)(end+1) = data;
				m = obj.Internal.component.map(name);
				m{end+1} = data;
				obj.Internal.component.map(name) = m;
            else
                % create a new dynamic property
                addprop(obj, name);
                obj.(name) = data;
				obj.Internal.component.map(name) = {data};
			end
        end
        
        
        function out = getComponents(obj, kind, out)
            % Returns all components
            
            if nargin<3
                out = [];
            end
            if nargin<2
                kind = 'all';
			end
			map = obj.Internal.component.map;
			keys = map.keys;
            for i = 1:length(keys)
				components = map(keys{i});
				for j = 1:length(components)
					component = components{j};
					if isa(component, 'ComponentBehavior')
						% dive recursively
						out = component.getComponents(kind, out);
					end
					if isa(component, 'SystemSignal')
						% TODO: test any other class with the isKind()
						% method (none at the moment besides SystemSignal)
						if isequal(kind, 'all') || component.isKind(kind)
							if numel(out)==0
								out = component;
							else
								out(end+1) = component;
							end
						end
					end
				end
            end
        end
        
        function out = listComponents(obj, kind)
            % Returns names of all components of a given kind
            
            if nargin < 2
                kind = 'all';
            end
            out = {};
            vars = obj.getComponents(kind);
            for i = 1:length(vars)
                out{end+1} = vars(i).name;
            end
        end

        function out = forEachComponent(obj, kind, fun, varargin)
            % Applies a given function to all components (recursively)
            
			if isequal(kind, 'all')
				components = obj.getComponents();
			else
				components = obj.getComponents(kind);
			end
            out = cell(1, length(components));
            for i = 1 :length(components)
                out{i} = feval(fun, components(i), varargin{:});
            end
		end
		
	end
	
	methods(Hidden)
		% private APIs

		function new = copyComponentsFrom(new, obj)
			% copies components of 'obj' to 'new'
			
			new.Internal.component.type = '';
			new.Internal.component.map = containers.Map;
			map = obj.Internal.component.map;
			keys = map.keys;
			for i = 1:length(keys)
				c = map(keys{i});
				for j = 1:length(c)
					new.addComponent(keys{i}, c{j}.copy());
				end
			end
		end

		function obj = saveAllComponents(obj)
			% saves arguments of all components to
			% obj.Internal.save.components
			
			components = obj.Internal.component.map.keys();
			arguments = [];
			for i = 1:length(components)
				arguments.(components{i}) = obj.Internal.component.map(components{i});
			end
			obj.Internal.save.components = arguments;
			% TODO: should we remove all components before saving,
			% similarly to how we remove filters in
			% FilterBehavior/saveAllFilters() ??
		end
		
		function obj = loadSavedComponents(obj)
			% re-initializes components using arguments stored in
			% obj.Internal.save.filters
			
			if isfield(obj.Internal, 'save') && ...
					isfield(obj.Internal.save, 'components') & ...
					isstruct(obj.Internal.save.components)
				arguments = obj.Internal.save.components;
				components = fieldnames(arguments);
				% clear components
				obj.Internal.component.type = '';
				obj.Internal.component.map = containers.Map;
				for i = 1:length(components)
					c = arguments.(components{i});
					for j = 1:length(c)
						obj.addComponent(components{i}, c{j});
					end
				end
			end
		end

	end
end
