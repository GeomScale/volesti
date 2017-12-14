function options = mpt_subModules
%% additional properties, depending on the directory structure

% list of modules represented by folders in this directory
% trying to detect modules by reading directory structure such that in
% future there should not be problems when adding new modules

d = fileparts(which('mpt_subModules.m'));

modules_list = dir(d);
% other modules to be added
%             list2 = {'analysis';
%                 'auxiliary';
%                 'control';
%                 'geometry';
%                 'graphics';
%                 'gui';
%                 'simulink'};

% for each module evaluate a function that specifies given
% options
for i=1:length(modules_list)
    if ~any(strcmpi(modules_list(i).name,{'.','..'})) && modules_list(i).isdir
        f=['mpt_',modules_list(i).name,'_options'];
        if ~exist(f,'file')
           error('Please, provide option specifications for "%s" m-file for "%s" module.',f,modules_list(i).name);
        end
        options.(modules_list(i).name) = feval(f);
    end
end

end