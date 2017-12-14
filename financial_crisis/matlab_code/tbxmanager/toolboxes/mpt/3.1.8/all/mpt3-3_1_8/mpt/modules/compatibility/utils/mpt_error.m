function mpt_error

global mptOptions MPTOPTIONS
if isempty(MPTOPTIONS)
	MPTOPTIONS=mptopt;
end

mptOptions=MPTOPTIONS.modules.compatibility;
