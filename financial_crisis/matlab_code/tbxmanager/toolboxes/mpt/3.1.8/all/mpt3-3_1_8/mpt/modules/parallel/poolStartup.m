function poolStartup
%
% Initialization of global settings for local worker when "matlabpool open"
% command opens up parallel computing channels. 
% Parallel Computing toolbox is needed to allow this feature.
%

try
    mpt_init;
end

end
