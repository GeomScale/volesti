classdef MPTUIHandle < dynamicprops & matlab.mixin.Copyable
    % Class representing a storage of components and filters
   
    properties(Hidden,SetAccess=protected)
        Internal
	end
	properties(Hidden,Constant)
		Version = 1.1
		% History:
		%   1.1 (June 5, 2013): "internal_properties" changed to "Internal"
	end

    methods
        function obj = MPTUIHandle()
		end
	end

	methods(Static)
		function new = loadobj(obj)
			% load method
			
			if ~isobject(obj)
				% objects saved in MPT 3.0.6 or older cannot be loaded due
				% to rename of "internal_properties" to "Internal"
				fprintf('Cannot load objects saved in MPT 3.0.6 and older.\n');
				new = obj;
				return
			else
				new = obj;
			end

			if isfield(obj.Internal, 'save')
				if isfield(obj.Internal.save, 'filters')
					% restore filters
					new.loadSavedFilters();
				end
				if isfield(obj.Internal.save, 'components')
					% restore components
					new.loadSavedComponents();
				end
				if isfield(obj.Internal.save, 'value')
					% restore signals' values (for SystemSignal objects
					% only)
					%
					% TODO: what about sdpvars introduced in filters?
					new.loadSdpvarValue();
				end
				
				new.Internal = rmfield(obj.Internal, 'save');
			end
		end
	end
	
	methods(Access=protected)

		function out = wasModified(obj)
			% Returns boolean flag indicating whether properties of the
			% object were changed by the user
			
			if ~isfield(obj.Internal, 'modified')
				obj.Internal.modified = true;
			end
			out = obj.Internal.modified;
		end
		
		function obj = markAsUnmodified(obj, source, event)
			% Marks the object as unmodified
			
			obj.Internal.modified = false;
		end

		function obj = markAsModified(obj, source, event)
			% Marks the object as unmodified
			
			obj.Internal.modified = true;
		end
		
		function new = copyElement(obj)
			% Copy constructor
			
			new = copyElement@matlab.mixin.Copyable(obj);

			if isfield(obj.Internal, 'component')
				% clone-copy components
				new.copyComponentsFrom(obj);
			end
			
			if isfield(obj.Internal, 'filter_map')
				% clone-copy filters
				new.copyFiltersFrom(obj);
			end
		end

	end

end
