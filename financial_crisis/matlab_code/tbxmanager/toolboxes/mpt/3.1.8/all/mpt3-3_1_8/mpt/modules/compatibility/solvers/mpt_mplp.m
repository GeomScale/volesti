function [Pn,Fi,Gi,activeConstraints,Phard,details]=mpt_mplp(Matrices,Options)
%MPT_MPLP Explicitly solves the given linear program (LP)
%
% [Pn,Fi,Gi,activeConstraints,Phard,details]=mpt_mplp(Matrices,Options)
%
% NOTE: This function is obsolete. Use the @Opt class to define and solve
% mpLP problems. 
%
% This method simply converts the parametric problem into an instance of
% the Opt class, which will solve it using the default parametric solver
% (either PLCP or MPLP).
%
% See "help mpt_mplp_26" for a detailed description of input/output
% arguments.

% Copyright is with the following author(s):
%
% (C) 2012 Michal Kvasnica, Slovak University of Technology in Bratislava
%     michal.kvasnica@stuba.sk

% ---------------------------------------------------------------------------
% Legal note:
%          This program is free software; you can redistribute it and/or
%          modify it under the terms of the GNU General Public
%          License as published by the Free Software Foundation; either
%          version 2.1 of the License, or (at your option) any later version.
%
%          This program is distributed in the hope that it will be useful,
%          but WITHOUT ANY WARRANTY; without even the implied warranty of
%          MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%          General Public License for more details.
% 
%          You should have received a copy of the GNU General Public
%          License along with this library; if not, write to the 
%          Free Software Foundation, Inc., 
%          59 Temple Place, Suite 330, 
%          Boston, MA  02111-1307  USA
%
% ---------------------------------------------------------------------------

% Since input and output arguments are identical and because the Opt class
% automatically determines whether we want to solve pLP or pQP, we simply
% re-use the code of the mpt_mpqp gateway.

[Pn,Fi,Gi,activeConstraints,Phard,details] = mpt_mpqp(Matrices);
