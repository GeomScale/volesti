function S = mpt_hysdel2_options
%
% Option settings for the "hysdel2" module.
%

% first check for hysdel.exe or hysdel.bin
if ispc
	S.hysdel_executable = which('hysdel.exe');
else
	S.hysdel_executable = which('hysdel.bin');
end

if isempty(S.hysdel_executable)
	% check current directory
	hysdel_path = fileparts(which(mfilename));
	if ispc
		hysdel_cmd = 'hysdel.exe';
	else
		hysdel_cmd = 'hysdel';
	end
	S.hysdel_executable = [hysdel_path filesep hysdel_cmd];
end

if ~exist(S.hysdel_executable, 'file')
	fprintf('Warning: hysdel compiler could not be found.');
end
