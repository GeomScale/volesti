function [info, entry_i] = syminfo(S, vname, vtype, vkind, vindex)
% info = syminfo(S, vname, vtype, vkind, vindex)
% extract infos from the symbol table
% look up the symbol table for entries matching the query 
% item.name == name AND item.type == type AND item.kind == kind
% empty field ('') matches all entries
% 
% INPUT:
% S     : the system name
% vname : var name as in the Hysdel source
% vtype : var type ('b','r')
% vkind : var kind ('x','u','d','z')
% vindex: var index
%
% OUTPUT
% info  : cell array of syminfo record satisfying the query
%         and index of the record
% 
% REMARK
% the contents of the syminfo records are detailed in the HYSDEL manual in the 
% section describing the MLD-structure
% 
% (C) 2001 by F.D. Torrisi
% Automatic Control Laboratory, ETH Zentrum, CH-8092 Zurich, Switzerland
% torrisi@aut.ee.ethz.ch
%
% added index to output, Tobias Geyer, 12. Jan. 2003
%
% see license.txt for the terms and conditions.

if nargin == 0
   disp('Version 0.01')
   info='0.01';
   return
end  

narginchk(1, 5);

% ==============================================================================
% Define optional arguments
% ==============================================================================

if nargin <= 4
   vindex = '';
end
if nargin <= 3
	vkind = '';
end
if nargin <= 2
   vtype = '';
end   
if nargin <= 1
   vname = '';
end

N_s = size(S.symtable,2);
rep = 0;
info = []; entry_i = [];
for i = 1:N_s
   if (isempty(vname) | strcmp(S.symtable{i}.name,vname)) & ...
      (isempty(vtype) | strcmp(S.symtable{i}.type,vtype)) & ...
      (isempty(vkind) | strcmp(S.symtable{i}.kind,vkind)) & ...
      (isempty(vindex) | S.symtable{i}.index == vindex) ,
      rep = rep + 1;
      
      info{rep} = S.symtable{i};
      entry_i(rep) = i;
   end
end
