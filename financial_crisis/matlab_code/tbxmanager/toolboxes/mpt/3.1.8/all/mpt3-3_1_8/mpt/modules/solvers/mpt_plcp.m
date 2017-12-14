function sol = mpt_plcp(opt)
%
%  MPT_PLCP: Parametric linear complementarity solver (PLCP) (without errorchecks) 
%  ================================================================================
%  
%  
%  SYNTAX
%  ------
%     
%      R = mpt_plcp(S)
%    
%  
%  DESCRIPTION
%  -----------
%     Implementation of the lexicographic PLCP solver. The PLCP is given as: 
%                                 w - Mz = q + Qtheta         (1) 
%                                              w >= 0         (2) 
%                                              z >= 0         (3) 
%                                       T                         
%                               w(theta) z(theta) = 0         (4) 
%                                                                 
%                              A     theta <= b               (5) 
%                               theta          theta              
%                                                             (6) 
%     where the matrices M, Q, A_theta, and vectors q, b_theta  are the problem
%  data, then z, w  are the decision variables and theta  are the parameters.
%  Structure of the algorithm: 
%    
%     1. Initialisation phase: For any feasible starting point theta_f solve
%     non-parametric LCP to get feasible basis. The basis is used to determine the
%     initial region P_0 with the corresponding optimal pair w and z.  
%     2. Exploration phase: For each facet of the initial region P_0 compute its
%     neighbors P_i by performing lex-pivot steps in the PLCP subproblem. Compute
%     the graph based on the found neighbors. Store basis of every region to a
%     corresponding hash-table. 
%     3. Termination phase: Verify the graph of the solution, evaluate the
%     statistics. 
%    The algorithm uses variable step approach for exploration of the parameter
%  space by default. The fixed step approach can be turned on via option
%  fixed_step. Parameter exploration is based on a graph search: depth-first search
%  (BFS) and breadth-first search (BFS) methods are considered, the default is BFS.
%  When the computer runs in a parallel mode, the exploration phase run in a
%  parallel for-loop automatically, which can potentially increases the
%  computational speed. Any change in default options must be done using mptopt
%  class. Following options can be modified: 
%    
%     - bfs - Logical value that determines if to use BFS for exploration of the
%     parameter space (default=1). 
%     - dfs - Logical value that determines if to use DFS for exploration of the
%     parameter space (default=0). 
%     - debug - Integer value that determines the debugging level (default=0). 
%     - maxlayers - For BFS-type of exploration, this value limits the number of
%     "layers" starting from the  initial region P_0. The first layer is defined as
%     the collection of regions that are adjacent to P_0, the second layer is given
%     as the new found neighbors to all regions in the first layer etc. The default
%     value is Inf. 
%     - maxregions - The maximal number of regions to be produced by PLCP
%     algorithm. This option is useful when solving large PLCP problems when there
%     are memory problems. The default value is Inf. 
%     - QRfactor - Logical value that select the type of factorization to use for
%     pivoting. If true, then recursive QR factorization is used instead of direct
%     sparse LU factorization by default. The default value is 0. 
%     - checkoverlaps - Logical value that launches overlap checking if true. Very
%     time consuming operation, therefore the default value is 0. 
%     - rescue - If the variable step approach fails to find all the neighbors,
%     this logical statement indicates if the use the fixes-step approach as a
%     backup. If true, then the fixed-step approach is run whenever variable step
%     approach does not  find overlaps. The disadvantage of fixed step is, however,
%     overlaps can be produced, specially for degenerate cases. The default value
%     is 0. 
%     - maxsteps - If the fixed step approach is turned on, this value limits the
%     number of steps to be performed to find  neighbor. The step size is given by
%     the region_tol. The default value is 200. 
%     - maxpivots - The maximum limit on the lex-pivot steps to be performed when
%     finding a neighbor. Typically, it suffices 1-2 pivots to find a neighbor. If
%     the problem is very degenerate, or badly posed, or due to numerical problems
%     involved in factorization, the pivot steps are repeated up to 100-times,
%     which is the default value. 
%     - cacheregions - This flag causes that regions that have been discovered are
%     remembered and used for faster exploration of the parameter space. The
%     default value is 1. 
%    Note that to properly preprocess data to PLCP, use Opt class whenever
%  possible. This will avoid unnecessary numerical problems caused by improper
%  formulation of the problem.
%  
%  INPUT
%  -----
%     
%        
%          S                  Structure of the Opt class.              
%                             Class: struct                            
%          S.M                Data matrix for linear-complementarity   
%                             problem w-Mx=q.                          
%                             Class: double                            
%          S.q                Right hand side vector for               
%                             linear-complementarity problem w-Mx=q.   
%                             Class: double                            
%          S.Ath              Linear part of the inequality            
%                             constraints A_thetatheta <= b_theta.     
%                             Class: double                            
%                             Default: []                              
%          S.bth              Right hand side of the inequality        
%                             constraints A_thetatheta <= b_theta.     
%                             Class: double                            
%                             Default: []                              
%          S.recover          Affine map for MPLP/MPQP problems after  
%                             transformation to LCP.                   
%                             Class: struct                            
%          S.recover.uX       Matrix of the affine map x = uX(         
%                              w                                       
%                              z  ) + uTh(                             
%                              theta                                   
%                                1    ) . The map is from the          
%                             optimization variables involed in LCP    
%                             w(theta),z(theta) |->> x  and in the     
%                             original LP/QP.                          
%                             Class: double                            
%                             Default: []                              
%          S.recover.uTh      Matrix of the affine map x = uX(         
%                              w                                       
%                              z  ) + uTh(                             
%                              theta                                   
%                                1    ) . The map is from the          
%                             optimization variables involed in LCP    
%                             w(theta),z(theta) |->> x  and in the     
%                             original LP/QP.                          
%                             Class: double                            
%                             Default: []                              
%          S.recover.lambdaX  Matrix of the affine map x = lambdaX(    
%                              w                                       
%                              z  ) + lambdaTh(                        
%                              theta                                   
%                                1    ) . The map is from the          
%                             optimization variables involed in LCP    
%                             w(theta),z(theta) |->> lambda  and the   
%                             Lagrangian multipliers in the original   
%                             LP/QP.                                   
%                             Class: double                            
%                             Default: []                              
%          S.recover.lambdaTh Matrix of the affine map x = lambdaX(    
%                              w                                       
%                              z  ) + lambdaTh(                        
%                              theta                                   
%                                1    ) . The map is from the          
%                             optimization variables involed in LCP    
%                             w(theta),z(theta) |->> lambda  and the   
%                             Lagrangian multipliers in the original   
%                             LP/QP.                                   
%                             Class: double                            
%                             Default: []                              
%                             Default: []                              
%          S.Internal         Internal data that came from             
%                             transformation from LP/QP to LCP.        
%                             Class: struct                            
%          S.Internal.H       Quadratic part of the objective          
%                             function.                                
%                             Class: double                            
%                             Default: []                              
%          S.Internal.f       Linear part of the objective function.   
%                             Class: double                            
%          S.Internal.pF      Linear part of the objective function    
%                             for parameters.                          
%                             Class: double                            
%                             Default: []                              
%                             Default: []                              
%                               
%  
%  
%  OUTPUT
%  ------
%     
%        
%          R                       Result structure                         
%                                  Class: struct                            
%          R.xopt                  Partition of the polyhedra with the      
%                                  associated function values for z  and w  
%                                  variables.                               
%                                  Class: PolyUnion                         
%          R.exitflag              An integer value that informs if the     
%                                  result was feasible (1), or otherwise    
%                                  (different from 1).                      
%                                  Class: double                            
%          R.how                   A string that informs if the result was  
%                                  feasible ('ok'), or if any problem       
%                                  appeared through optimization.           
%                                  Class: char                              
%          R.solveTime             How long did the computation take in     
%                                  seconds.                                 
%                                  Class: double                            
%          R.stats                 Statistical information about the        
%                                  computation: the total number of pivots  
%                                  performed, the total number of facets    
%                                  traversed.                               
%                                  Class: struct                            
%          R.stats.pivs            The total number of pivots performed.    
%                                  Class: double                            
%          R.stats.facetsTraversed The total number of facets that have     
%                                  been traversed.                          
%                                  Class: double                            
%                                    
%  
%  
%  SEE ALSO
%  --------
%     mpt_solvemp,  lcp
%  

%  AUTHOR(s)
%  ---------
%     
%    
%   (c) 2010-2013  Colin Neil Jones: EPF Lausanne
%   mailto:colin.jones@epfl.ch 
%     
%    
%   (c) 2010-2013  Martin Herceg: ETH Zurich
%   mailto:herceg@control.ee.ethz.ch 
%  
%  

%  LICENSE
%  -------
%    
%    This program is free software; you can redistribute it and/or modify it under
%  the terms of the GNU General Public License as published by the Free Software
%  Foundation; either version 2.1 of the License, or (at your option) any later
%  version.
%    This program is distributed in the hope that it will be useful, but WITHOUT
%  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
%  FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
%    You should have received a copy of the GNU General Public License along with
%  this library; if not, write to the  Free Software Foundation, Inc.,  59 Temple
%  Place, Suite 330,  Boston, MA 02111-1307 USA
%  -------------------------------------------------------------------------------
%    
%      This document was translated from LaTeX by HeVeA (1).
%  ---------------------------------------
%    
%    
%   (1) http://hevea.inria.fr/index.html
 
 
global MPTOPTIONS BNDH INPUT_LCP
if isempty(MPTOPTIONS)
    MPTOPTIONS = mptopt;
end

% halfspace representation of the feasible set
BNDH = zeros(0, opt.d+1);

% copy of the input problem (used in is_boundary_facet())
INPUT_LCP = opt.copy();

% pick DFS or BFS but not both
if MPTOPTIONS.modules.solvers.plcp.bfs
	if MPTOPTIONS.modules.solvers.plcp.dfs
		disp('Both BFS and DFS options are turned on, using only BFS.');
		MPTOPTIONS.modules.solvers.plcp.dfs = false;
	end
elseif MPTOPTIONS.modules.solvers.plcp.dfs
	if MPTOPTIONS.modules.solvers.plcp.bfs
		disp('Both BFS and DFS options are turned on, using only DFS.');
		MPTOPTIONS.modules.solvers.plcp.bfs = false;
	end
elseif ~MPTOPTIONS.modules.solvers.plcp.bfs && ~MPTOPTIONS.modules.solvers.plcp.dfs
    % default
    MPTOPTIONS.modules.solvers.plcp.dfs = true;
end
% correct maxsteps in case is less than 2
if MPTOPTIONS.modules.solvers.plcp.maxsteps<2
    MPTOPTIONS.modules.solvers.plcp.maxsteps=2;
end
% correct maximum number of pivots
if MPTOPTIONS.modules.solvers.plcp.maxpivots<1
    MPTOPTIONS.modules.solvers.plcp.maxpivots = 1;
end
if MPTOPTIONS.modules.solvers.plcp.maxpivots>get(0,'RecursionLimit')
    error('The "maxpivots" setting for PLCP solver is bigger than "RecursionLimit". Please, change the "RecursionLimit" via "set(0,''RecursionLimit'',N)"');    
end


% use to paralellize the code
% check if the matlab pool is running
if ~isempty(which('parpool'))
    % R2014b and newer
    poolobj = gcp('nocreate'); % If no pool, do not create new one.
    if isempty(poolobj)
        poolsize = 0;
    else
        poolsize = poolobj.NumWorkers
    end
    ISPARALLEL = poolsize>0;
elseif ~isempty(which('matlabpool'))
    ISPARALLEL = matlabpool('size')>0;
else
    ISPARALLEL = false;
end

if ISPARALLEL
    disp('Running in parallel');
end

% if there are some equality constraints present, throw error
if opt.me>0
   error('PLCP solver does not solve problems with equality constraints. Try a different solver.');
end

% put LCP arguments together
lc.M  = opt.M;
lc.Q  = opt.Q;
lc.q  = opt.q;
lc.Ht = [opt.Ath opt.bth];

% lc = varargin2struct(varargin{:});

lc.d = size(lc.Q,2);
lc.n = size(lc.M,1);
lc.A = [eye(lc.n) -lc.M];
lc.recover = opt.recover;
lc.obj = opt.Internal;
if ~isempty(lc.recover)
    lc.P = speye(size(lc.recover.uX,1));
    % reorder variables is this came from Yalmip
    if ~isempty(opt.varOrder)
        lc.P = lc.P(opt.varOrder.requested_variables,:);
    end
end

% if isfield(lc, 'bndA') && ~isempty(lc.bndA),
%   lc.Ht = [lc.bndA lc.bndb];
% elseif ~isfield(lc, 'Ht')
%   lc.Ht = zeros(0, lc.d+1);
% end

% Initialize hash table
HASH = containers.Map;

%% region exploration

tStart = clock; tic;
iter = 0;
layer = 2;

% first layer is always the first region
layer_list{1} = 1; 

% pick the initial feasible parameter theta
fo = lcp_initialfeas(lc);


% Compute an initial basis
%[w,z,Io,result,piv]=lcp_mex(lc.M,lc.q+lc.Q*fo);
[z, w, Io, result, piv] = lcp(lc.M, lc.q+lc.Q*fo, MPTOPTIONS.modules.solvers.lcp);

if (result ~= MPTOPTIONS.OK)
    % pick again the initial feasible parameter theta
    fo = lcp_initialfeas(lc);
    [z, w, Io, result, piv] = lcp(lc.M, lc.q+lc.Q*fo, MPTOPTIONS.modules.solvers.lcp);
    if (result~=MPTOPTIONS.OK)        
        % retry with a matlab version
        [z,w,Io,result,piv] = lcpm(lc.M,lc.q+lc.Q*fo);
        if (result ~= MPTOPTIONS.OK)
            sol.xopt = PolyUnion;
            sol.xopt.setInternal('theta0', fo);
            sol.exitflag = MPTOPTIONS.INFEASIBLE;
            sol.how = 'infeasible';
            sol.stats.pivs = piv;
            sol.stats.facetsTraversed = 0;
            sol.stats.time = etime(clock, tStart);
            return
            %error('Initial parameter theta is not feasible.');
        end
    end
end

% sort the indices before adding to hash table
Io = sort(Io(:))';

R = lcp_getRegion(lc, Io, HASH, [], [], piv);

% The region returned by the first shot might be empty. This happens
% usually around origin for LP. We retry to find another region by
% perturbing initial theta by a small value.
if builtin('isempty',R) || R.isEmptySet || ~R.isFullDim || ~R.isBounded
    % try using the interior point approach
    P = Polyhedron('H', [-lc.Q -lc.M lc.q; zeros(lc.n,lc.d) -eye(lc.n) zeros(lc.n,1)]);
    ip = P.interiorPoint;
    if isempty(ip.x)
        sol.xopt = PolyUnion;
        sol.xopt.setInternal('theta0', fo);
        sol.exitflag = MPTOPTIONS.INFEASIBLE;
        sol.how = 'infeasible';
        sol.stats.pivs = piv;
        sol.stats.facetsTraversed = 0;
        sol.stats.time = etime(clock, tStart);
        return;
        %error('Initial parameter theta is not feasible.');
    end
    
    % first d variables are theta
    fo = ip.x(1:lc.d);

    [z, w, Io, result, piv] = lcp(lc.M, lc.q+lc.Q*fo, MPTOPTIONS.modules.solvers.lcp);
    if (result ~= MPTOPTIONS.OK)
        % slightly perturb and try again
        fo = fo + 1e-2*randn(lc.d,1);
        [z, w, Io, result, piv] = lcp(lc.M, lc.q+lc.Q*fo, MPTOPTIONS.modules.solvers.lcp);
        if result~=MPTOPTIONS.OK
            % retry with a matlab version
            [z,w,Io,result,piv] = lcpm(lc.M,lc.q+lc.Q*fo);
            if (result ~= MPTOPTIONS.OK)
                sol.xopt = PolyUnion;
                sol.xopt.setInternal('theta0', fo);
                sol.exitflag = MPTOPTIONS.INFEASIBLE;
                sol.how = 'infeasible';
                sol.stats.pivs = piv;
                sol.stats.facetsTraversed = 0;
                sol.stats.time = etime(clock, tStart);
                return;
                %error('Initial parameter theta is not feasible.');
            end
        end
    end
    
    % sort the indices before adding to hash table
    Io = sort(Io(:))';
    
    R = lcp_getRegion(lc, Io, HASH, [], [], piv);

    if builtin('isempty',R) || R.isEmptySet || ~R.isFullDim || ~R.isBounded
        % if fo is exactly zero, the region might be empty, therefore we
        % perturb a bit and try again
        fo = fo + 1e-2*randn(lc.d,1);
        [z, w, Io, result, piv] = lcp(lc.M, lc.q+lc.Q*fo, MPTOPTIONS.modules.solvers.lcp);
        if result~=MPTOPTIONS.OK
            % retry with a matlab version
            [z,w,Io,result,piv] = lcpm(lc.M,lc.q+lc.Q*fo);
            if (result ~= MPTOPTIONS.OK)
                sol.xopt = PolyUnion;
                sol.xopt.setInternal('theta0', fo);
                sol.exitflag = MPTOPTIONS.INFEASIBLE;
                sol.how = 'infeasible';
                sol.stats.pivs = piv;
                sol.stats.facetsTraversed = 0;
                sol.stats.time = etime(clock, tStart);
                return;                
                %error('Initial parameter theta is not feasible.');
            end
        end
        % sort the indices before adding to hash table
        Io = sort(Io(:))';
        
        R = lcp_getRegion(lc, Io, HASH, [], [], piv);
        if builtin('isempty',R) || R.isEmptySet || ~R.isFullDim || ~R.isBounded
            % unfortunately, now we must capitulate here
            sol.xopt = PolyUnion;
            sol.xopt.setInternal('theta0', fo);
            sol.exitflag = MPTOPTIONS.INFEASIBLE;
            sol.how = 'infeasible';
            sol.stats.pivs = piv;
            sol.stats.facetsTraversed = 0;
            sol.stats.time = etime(clock, tStart);
            return;            
            %error('Initial region is empty or not bounded.');
        end
    end
end

% number of estimated regions
est_regions = 1;

% array of unexplored regions
UNEX = R;

% store index set I and data determining the region R
hdat.basis = R.Internal.I;
hdat.redundant_rows = R.Internal.redundant_rows;
HASH = map_put(HASH, Io, hdat);

% initialize counters
regions = Polyhedron.emptySet(lc.d);
stats.pivs = 0;
stats.facetsTraversed = 0;
flag = 1;
how = 'ok';

if MPTOPTIONS.verbose >= 1
    hw = waitbar(0,'1','Name','Solving PLCP, please wait ...',...
        'CreateCancelBtn','setappdata(gcbf,''interrupt'',1)');
    setappdata(hw,'interrupt',0);
end

while ~builtin('isempty',UNEX)
    
  if toc > MPTOPTIONS.report_period
	  fprintf('regions: %4i, unexplored: %i \n', ...
		  length(regions), length(UNEX));
	  tic;
  end
  
  
  % Pop an unexplored region off the stack
  R    = UNEX(1);
  UNEX = UNEX(2:end);

  if MPTOPTIONS.modules.solvers.plcp.checkoverlaps %|| isfield(R.Internal,'th')
      % check if it is not overlapping existing regions
      [ts,iov] = check_overlaps(R,regions);
      
      if ts
          if (MPTOPTIONS.verbose >= 2 ) || (MPTOPTIONS.modules.solvers.plcp.debug >= 1)
              disp('Overlapping region discovered.');
          end
          % store index set and solution
          Ia = R.Internal.I;
          %xa = R.Internal.x;
          pivx = R.Internal.piv;
          redx = R.Internal.redundant_rows;
          
          % do the cut
          Rd = R\regions(iov);
          if ~builtin('isempty',Rd)              
             Rd(Rd.isEmptySet) = []; 
             Rd(~Rd.isFullDim) = [];
             Rd(~Rd.isBounded) = [];
			 Rd.copyFunctionsFrom(R);
          end
          
          for i=1:length(Rd)
              % correct for case when Rd==R (due to bad numerics) to avoid loop
              if ~(Rd(i)==R)
                  % regions might have redundant facets
                  Rd(i).minHRep;
                  
                  Rd(i).setInternal('I',Ia);
                  %Rd(i).setInternal('x',xa);
                  Rd(i).setInternal('piv',pivx);
                  Rd(i).setInternal('redundant_rows',redx);
                  
                  % put at the beginning in order to explore it next after this cut
                  UNEX = [Rd(i); UNEX];
                  
                  % update layers
                  if iter>max(layer_list{layer-1})
                      layer = layer + 1;
                  end
              end
          end
          continue
      end
  end
  

  % check if it is not empty
  if (MPTOPTIONS.verbose >= 2 ) || (MPTOPTIONS.modules.solvers.plcp.debug >= 1)
      if R.isEmptySet || ~R.isBounded
          disp('Empty or unbounded region discovered ... skipping.');
          continue;
      end
  end

     
  % Store solution
  iter = iter + 1; % current iteration
  actual_region = map_getIndex(HASH, R.Internal.I);
  regions(actual_region) = R; % store regions in array as they were explored
  er = ([regions.Dim]==0); % empty regions
  regions(er) = Polyhedron.emptySet(lc.d); % set the dimensio of empty regions to lc.d
  adj_list{actual_region} = cell(size(R.H,1),1); % adjacency list
  layer_list{layer+1} = []; % layers list 
  regions_found = 0;
         
  if MPTOPTIONS.modules.solvers.plcp.debug >= 2
      drawnow
      axis tight;
      fprintf('Press any key to continue\n');
      pause
      clf;
      hold on; grid on;
      % actual region in red
      R.plot('color','r','linewidth',1);
      
      % plot initial region in light blue
      plot(regions(1),'color','b','alpha',0.3);
      % add numbers to regions
      ip1 = regions(1).chebyCenter; 
      if lc.d==1          
          text(ip1.x,0,num2str(1));
      elseif lc.d==2
          text(ip1.x(1),ip1.x(2),num2str(1));
      end
      
      % plot remaining regions
      if length(regions)>1
          nr = 2:length(regions);
        plot(regions(nr),'alpha',0);
        for ii=nr
           ip = regions(ii).chebyCenter;
           if ~regions(ii).isEmptySet
               if lc.d==1
                   text(ip.x,0,num2str(map_getIndex(HASH, regions(ii).Internal.I)));
               elseif lc.d==2
                   text(ip.x(1),ip.x(2),num2str(map_getIndex(HASH, regions(ii).Internal.I)));
               end
           end
        end        
      end
  end
    
   
  % for each facet of the current regions, find its neighbors
  Radj = cell(1,size(R.H,1));
  if ISPARALLEL
      parfor i = 1:size(R.H,1)
                    
          if MPTOPTIONS.modules.solvers.plcp.fixedstep %FIXEDSTEP
              %always perform fixed step over the facet
              Radj{i} = lcp_ComputeAdj_fixedstep(R,i,lc, HASH, regions, UNEX, ISPARALLEL);
          else
              % compute the step size via pivoting from R
              Radj{i} = lcp_ComputeAdj(R,i,lc, HASH, regions, UNEX, ISPARALLEL);
              
              if MPTOPTIONS.modules.solvers.plcp.rescue
                  % if the variable step approach returns empty, but feasible
                  % region, try again with fixed step
                  if builtin('isempty',Radj{i}) || all(Radj{i}.isEmptySet)
                      Radj{i} = lcp_ComputeAdj_fixedstep(R,i,lc, HASH, regions, UNEX, ISPARALLEL);
                  end
              end
          end
          
      end

  else
      for i = 1:size(R.H,1),
          
          if MPTOPTIONS.modules.solvers.plcp.debug >= 2
              % plot actual facet
              exci= setdiff(1:size(R.H,1),i)';
              Pfacet=Polyhedron('He',R.H(i,:),'H',R.H(exci,:));
              Pfx = Pfacet.chebyCenter;
              if ~isempty(Pfx.x)
                  if lc.d==1
                      text(Pfx.x,0,['f',num2str(i)]);
                  elseif lc.d==2
                      text(Pfx.x(1),Pfx.x(2),['f',num2str(i)]);
                  end
              end
              
              fprintf('Press any key to continue\n');
              pause
          end
          
          if MPTOPTIONS.modules.solvers.plcp.fixedstep
              %always perform fixed step over the facet
              Radj{i} = lcp_ComputeAdj_fixedstep(R,i,lc, HASH, regions, UNEX, ISPARALLEL);
          else
              % compute the step size via pivoting from R
              Radj{i} = lcp_ComputeAdj(R,i,lc, HASH, regions, UNEX, ISPARALLEL);
              
              if MPTOPTIONS.modules.solvers.plcp.rescue
                  % if the variable step approach returns empty, but feasible
                  % region, try again with fixed step
                  if builtin('isempty',Radj{i}) || all(Radj{i}.isEmptySet)
                      if (MPTOPTIONS.verbose>=2) || (MPTOPTIONS.modules.solvers.plcp.debug>=1)
                          disp('Neighbor is empty, correcting with fixed step over the facet.');
                      end
                      Radj{i} = lcp_ComputeAdj_fixedstep(R,i,lc, HASH, regions, UNEX, ISPARALLEL);
                  end
              end
          end
          
      end
  end
  
  for i=1:length(Radj)
      % For each new found region, put the description given by index I to hash table
      % Note that different index sets can give the same region!
      for j=1:length(Radj{i})
                         
                    
          % do not add it to the list of unexplored regions if is
          % empty, lower-dimensional, unbounded or it overlaps existing
          % regions                    
          if ~Radj{i}(j).isEmptySet && Radj{i}(j).isFullDim && Radj{i}(j).isBounded

              % put the discovered basis into the hash table regardless of
              % the overlap
              hdat.basis = Radj{i}(j).Internal.I;
              hdat.redundant_rows = Radj{i}(j).Internal.redundant_rows;
              [HASH, S] = map_put(HASH, Radj{i}(j).Internal.I,hdat);
              
              % count discovered regions
              if (S==0)
                  regions_found = regions_found + 1;
              end
              
              % build graph represented by adjacency list              
              region_index = map_getIndex(HASH, Radj{i}(j).Internal.I);
                           
              % for each facet i we found j neighbors
              adj_list{actual_region}{i}(j) = region_index;
              
              % DFS
              if MPTOPTIONS.modules.solvers.plcp.dfs
                  if (S == 0)
                      UNEX = [Radj{i}(j);UNEX];
                  end;
              end
              
              % BFS
              if MPTOPTIONS.modules.solvers.plcp.bfs
                  if (S == 0)
                      UNEX = [UNEX;Radj{i}(j)];
                      % check also which layer do we explore right now
                      if iter>max(layer_list{layer-1})
                          layer = layer + 1;
                      end
                      layer_list{layer} = [layer_list{layer}; region_index];
                      
                  end
              end
              
              
              % do some plots for debugging purposes
              if MPTOPTIONS.modules.solvers.plcp.debug >= 2
                  % plot discovered neighbors
                  Radj{i}(j).plot('color','g','alpha',0.2);
                  % show legend
                  legend(-1,'actual region','initial region')
              end
          end
      end
    
    % if there are zero indices (empty neighbor), we need to remove them
    to_delete = adj_list{actual_region}{i}<=0;
    if ~isempty(to_delete)
        adj_list{actual_region}{i}(to_delete)=[];
    end
  end
  
  if MPTOPTIONS.verbose >= 1
      % update number of estimated regions
      est_regions = est_regions + regions_found;
      % update waitbar
      if MPTOPTIONS.modules.solvers.plcp.bfs
          waitbar(iter/est_regions,hw,sprintf('Region %d from %d estimated. Level %d.',iter,est_regions,layer));
      else
          waitbar(iter/est_regions,hw,sprintf('Region %d from %d estimated.',iter,est_regions));
      end
      if getappdata(hw,'interrupt')
          break;
      end
  end

  % break in case of maximum layers
  if layer>=MPTOPTIONS.modules.solvers.plcp.maxlayers && MPTOPTIONS.modules.solvers.plcp.bfs
      disp('Maximum number of allowed layers achieved. Interrupting ...');
      flag = -3;
      how = 'maximum number of layers achieved';
      break;
  end

  % break in case of maximum regions
  if iter>=MPTOPTIONS.modules.solvers.plcp.maxregions
      disp('Maximum number of regions achieved. Interrupting ...');
      flag = -2;
      how = 'maximum number of regions achieved';
      break;
  end


  stats.pivs = stats.pivs + R.Internal.piv;
  stats.facetsTraversed = stats.facetsTraversed + size(R.H,1);
  
end
% end of region exploration

%% evaluate and post process solution

time = etime(clock, tStart);

fprintf('mpt_plcp: %d regions\n',length(regions));
if MPTOPTIONS.verbose >= 1
    % close waitbar
    delete(hw);

    % show statistics
    fprintf('Total pivots / facets  : %3i / %3i\n', stats.pivs, stats.facetsTraversed);
    fprintf('Average pivots / facet : %.2f\n', stats.pivs / stats.facetsTraversed);
    fprintf('Computation time : %.2f seconds\n', time);
end

% Build the output functions
if (MPTOPTIONS.verbose >= 1) || (MPTOPTIONS.modules.solvers.plcp.debug >=1)
    fprintf('Post-processing...');
end
% delete last field which is empty
layer_list(end) = [];

adj_list_verified = false;
if MPTOPTIONS.modules.solvers.plcp.adjcheck
    % check graph, if neighbors are correct (if not, remove those neighbors
    % who do not link to each other)
    adj_list = verify_graph(regions,adj_list);
    adj_list_verified = true;
end

% compute the set of feasible parameters: first try the boundary facet
% identified during region exploration
hull = Polyhedron('H', BNDH).normalize().minHRep();

if hull.isEmptySet && any(~regions.isEmptySet)
    if MPTOPTIONS.verbose>=1
        fprintf('Feasible set is empty despite non-empty regions, recomputing via the adjacency list...\n');
    end
    hull = feasible_set(regions, adj_list, opt);
    
    % quick sanity check: does the feasible set contain the chebycenters of all
    % regions?
    if ~sub_check_feasible_set(hull, regions)
        
        if adj_list_verified
            % the adjacency list was verified previously, yet the feasible set
            % is wrong. Declare it as an empty set such that it is recomputed
            % by projection lated
            hull = Polyhedron.emptySet(opt.d);
        else
            % verify the adjacency list
            if MPTOPTIONS.verbose>=0
                fprintf('Fixing the adjacency list...\n');
            end
            adj_list = verify_graph(regions,adj_list);
            adj_list_verified = true;
            if MPTOPTIONS.verbose>=0
                fprintf('...done.\n');
            end
            
            % and reconstruct the feasible set based on the fixed list
            hull = feasible_set(regions, adj_list, opt);
            
            % check correctness of the feasible set
            if ~sub_check_feasible_set(hull, regions)
                % the feasible set is still wrong
                hull = Polyhedron.emptySet(opt.d);
            end
        end
    end
end

if hull.isEmptySet && any(~regions.isEmptySet)
    % feasible set is empty which means that the adjacency
    % list is not correct and there could be holes in the solution
    flag = -4;
    how = 'wrong adjacency list';
    if MPTOPTIONS.verbose>=0
        fprintf('Feasible set is wrong, returning outer approximation.\n');
    end
    hull = PolyUnion(regions).outerApprox();
end

% normalize regions to return normalized solution
regions.normalize;

if flag==1
    % solution terminated correctly
    sol.xopt = PolyUnion('Set',regions,'Convex',true,'Overlaps',false,...
    'Bounded',true,'Connected',true,'FullDim',true,'Domain', hull);
else
    % preterminated, can be nonconvex
    sol.xopt = PolyUnion('Set',regions,'Overlaps',false,...
    'Bounded',true,'Connected',true,'FullDim',true, 'Domain', hull);
end

% store feasible set 
sol.xopt.setInternal('convexHull', hull);

% store other important properties under Internal
sol.xopt.setInternal('hash_table_basis', HASH);
sol.xopt.setInternal('adj_list', adj_list);
sol.xopt.setInternal('adj_list_verified', adj_list_verified);
if MPTOPTIONS.modules.solvers.plcp.bfs
    sol.xopt.setInternal('layer_list', layer_list);
else
    sol.xopt.setInternal('layer_list', {[]});
end
sol.xopt.setInternal('bfs', MPTOPTIONS.modules.solvers.plcp.bfs);
sol.xopt.setInternal('dfs', MPTOPTIONS.modules.solvers.plcp.dfs);
sol.xopt.setInternal('theta0', fo);
% map to original variables
sol.xopt.setInternal('recover',opt.recover);

    
% sol.regions = Polyset(regions);
% sol.x       = PolySet;
% for i=1:length(regions)
%     sol.x.add(PolyAffineFunc(regions(i), 'funcDat', regions(i).Internal.x));
% end

% output from PLCP solver
sol.exitflag = flag;
sol.how = how;
stats.solveTime = time;
sol.stats = stats;


if (MPTOPTIONS.verbose >= 1) || (MPTOPTIONS.modules.solvers.plcp.debug >=1)
    fprintf('done.\n')
end

end



%% subfunctions
function answer = sub_check_feasible_set(K, regions)
% checks whether the feasible set K contains the chebycenter of each region

answer = true;
ch = regions.chebyCenter;
for i = 1:length(ch)
    if ~isempty(ch(i).x) && ~K.contains(ch(i).x)
        % the feasible set does not contain the chebyceter of the i-th
        % region, thus it is wrong
        answer = false;
        return
    end
end

end

%---
function fo = lcp_initialfeas(lc)
%
% Compute any feasible value of theta s.t. (M,q+Q th) is feasible
%

% load global options
global MPTOPTIONS
    
    % find an approximate solution by relaxing complementarity constraints
    % -Q*th + A*x = q
    %  Ath*th <= bth
    %   x >= 0
    
    S.Ae = [-lc.Q lc.A];
    S.be = lc.q;
    S.A = [lc.Ht(:,1:end-1) zeros(size(lc.Ht,1),2*lc.n); zeros(2*lc.n,lc.d) -eye(2*lc.n)];
    S.b = [lc.Ht(:,end); zeros(2*lc.n,1)];
    % instead of solving pure feasibility problem (which turns to be
    % numerically problematic), try to minimize all variables to some 
    % small reasonable value
    S.f = ones(2*lc.n+lc.d,1)+0.01*rand(2*lc.n+lc.d,1);
	S.lb = []; S.ub = []; S.quicklp = true;
    r = mpt_solve(S);
    
    if r.exitflag~=MPTOPTIONS.OK
        % retry again
        S.f = ones(2*lc.n+lc.d,1)+0.01*rand(2*lc.n+lc.d,1);
        r = mpt_solve(S);
        %error('mpt_plcp: Parameter set is empty');
    end
    
    if r.exitflag~=MPTOPTIONS.OK
        % infeasible
        fo = zeros(lc.d,1);
    else
        % feasible theta are first lc.d variables
        fo = r.xopt(1:lc.d);
    end

end

%%
function R = lcp_getRegion(lc, I, HASH, regions, UNEX, piv, computeHull)
% 
% For given index set I of the basis lc, return a corresponding Polyhedron.
% If there are only 6 arguments, perform also a redundancy elimitation
% such that returned Polyhedron is in its minimal representation.
%
% HASH and regions arguments are used to speed up the region computation
% if the region was discovered and stored earlier
%

global MPTOPTIONS

if nargin < 7
    computeHull = true;
end
if nargin < 6
    piv = 0;
end
 
% check if region already exist
s = map_getIndex(HASH, I);
% region exist, get it from the list and return quickly
if ~isempty(s)
    nR = numel(regions);
    if s<=nR
        R = regions(s);
    else
        R = UNEX(s-nR);
    end
    if ~R.isEmptySet
        % region exist, but it is not stored in the arrays of regions,
        % we need to evaluate it from the polyhedral description
        return
    end
end

% region description
A = lc.A(:,I)\[-lc.Q lc.q];

% if the matrix cannot be inverted, return empty region
if any(any(isnan(A),2))
    R = [];
    return
end
% intersect current representation and bounds for theta
Hn = [A; lc.Ht];

% create polyhedron {x | A*x<=b}
% 
% Note that the polyhedron constructor retains zero rows
R = Polyhedron(Hn(:, 1:end-1), Hn(:, end));

% detect zero rows
zero_rows = sum(abs(R.A),2)<MPTOPTIONS.zero_tol;
% scaling matrix
Hscaled = Hn;
% detect zero columns
Hnz = Hn(~zero_rows,:);
zero_cols = sum(abs(Hnz),1)<MPTOPTIONS.zero_tol;
% do not scale if there is any zero column
if ~any(zero_cols)
    [An,D1,D2] = mpt_scale_matrix(Hnz);
    % D1 must be nonnegative
    if (all(diag(D1)>MPTOPTIONS.zero_tol))
        % respect zero rows in the original system
        Hscaled(~zero_rows,:) = D1*Hn(~zero_rows,:);
        % create a scaled polyhedron to avoid numerical problems in redudancy
        % elimination and detection of empty polyhedra
        R = Polyhedron(Hscaled(:, 1:end-1), Hscaled(:, end));
    end
end

% exit quickly if the polyhedron is not fully dimensional
if ~R.isFullDim()
    R = [];
    return;
end

% store additional data inside Internal properties
R.setInternal('I',I);
x = zeros(lc.n*2,lc.d+1); % Solution
A(:,1:lc.d) = -A(:,1:lc.d);
x(I,:) = A;
%R.setInternal('x',x); % Solution
R.setInternal('piv',piv); % pivots

% set function handles for w and z variables
Lz = AffFunction(x(lc.n+1:2*lc.n,1:end-1),x(lc.n+1:2*lc.n,end));
Lw = AffFunction(x(1:lc.n,1:end-1),x(1:lc.n,end));
R.addFunction(Lz,'z');
R.addFunction(Lw,'w');

% if LCP was created from LP/QP, compute also primal, dual variables and
% the value function
if ~isempty(lc.recover)
    TT = lc.recover.uX*x + lc.recover.uTh;
    
    % compute primal variables
    Lprimal = lc.P*TT;
    R.addFunction(AffFunction(Lprimal(:,1:end-1),Lprimal(:,end)),'primal');
    
    % compute dual variables for inequalities
    %Ldual = lc.recover.lambdaX*x + lc.recover.lambdaTh;
    Lineq = lc.recover.lambda.ineqlin.lambdaX*x + lc.recover.lambda.ineqlin.lambdaTh;
    R.addFunction(AffFunction(Lineq(:,1:end-1),Lineq(:,end)),'dual-ineqlin');
    
    % dual variables for equalities
    Leq = lc.recover.lambda.eqlin.lambdaX*x + lc.recover.lambda.eqlin.lambdaTh;
    R.addFunction(AffFunction(Leq(:,1:end-1),Leq(:,end)),'dual-eqlin');
    
    % dual variables for lower bounds
    Llb = lc.recover.lambda.lower.lambdaX*x + lc.recover.lambda.lower.lambdaTh;
    R.addFunction(AffFunction(Llb(:,1:end-1),Llb(:,end)),'dual-lower');
    
    % dual variables for upper bounds
    Lub = lc.recover.lambda.upper.lambdaX*x + lc.recover.lambda.upper.lambdaTh;
    R.addFunction(AffFunction(Lub(:,1:end-1),Lub(:,end)),'dual-upper');
    
    % compute the objective value
    Y = TT(:,1:end-1);
    T = TT(:,end);
    
    if ~isempty(lc.obj.H)
        qt = 0.5*Y'*lc.obj.H*Y + lc.obj.pF'*Y + lc.obj.Y;
        lt = T'*lc.obj.H*Y + T'*lc.obj.pF + lc.obj.f'*Y + lc.obj.C;
        at = 0.5*T'*lc.obj.H*T + lc.obj.f'*T + lc.obj.c;
        R.addFunction(QuadFunction(qt,lt,at),'obj');
    else
        qt = lc.obj.pF'*Y + lc.obj.Y;
        lt = T'*lc.obj.pF + lc.obj.f'*Y + lc.obj.C;
        at = lc.obj.f'*T + lc.obj.c;
        if any(any(qt))
            R.addFunction(QuadFunction(qt,lt,at),'obj');
        else
            R.addFunction(AffFunction(lt,at),'obj');
        end
    end
end


if computeHull
    [R, hull] = R.minHRep();
    
	% store indices of redundant rows
	R.setInternal('redundant_rows', hull.I);
end

end


function [ts, io] = check_overlaps(R,regions)
% check if region R overlaps the array of regions and return true or false
%
% returns alson the region indices which overlap 

global MPTOPTIONS

ts = false;
io=[];
for i=1:length(regions)   
    P = intersect(R,regions(i));
    % rather test with rel_tol
    %if ~P.isEmptySet
    if P.chebyCenter.r>MPTOPTIONS.rel_tol 
        ts = true;
        io=[io,i];
    end    
end

end

function Radj = lcp_ComputeAdj_fixedstep(R,i,lc, HASH, regions, UNEX, ISPARALLEL)
%
%  computes an adjacent region by considering a small step beyond facet
%  with a fixed size
%
% inputs: R - current region
%         i - index of a facet to step over
%         lc - structure with LCP data (M, q, Q, Ath, bth)
%         HASH - hash table with optimal bases defining polyhedra
%         regions - an array of already discovered regions
%         UNEX - an array of unexplored regions
%
% outputs:  Radj - adjacent regions


global MPTOPTIONS

% extract information from the facet equation (normal to the facet)
A = R.A; b = R.b;
gamma = A(i,:);

% normalize facet equation s.t. ||gamma||_2 = 1
ng = norm(gamma);
g = gamma / ng;
h = b(i) / ng;

% get a feasible point on the facet
xc = chebyCenter(R,i);
thf = xc.x;

% if we cannot find feasible point on the facet, leave
if xc.exitflag~=MPTOPTIONS.OK
    Radj = [];
    return; 
end

% compute new theta as thnew = thf + g'*step_size
% if the ste size is too small then too many steps are needed to find a new
% feasible basis, otherwise, if too big, then we might skip over a region.
thn = thf+g'*MPTOPTIONS.region_tol/2;

% use feasible basis from region R to reduce the number of pivots to detect
% neighboring basis I
[Radj, piv] = lcp_getRegionFromTheta(thn, lc, R.Internal.I, HASH, regions, UNEX);
R.setInternal('piv', R.Internal.piv + piv);
if ~ISPARALLEL && isempty(Radj)
    % possible boundary
    update_boundary(R, i);
end

% if Radj is considered as empty region (but feasible), we need to adjust
% the step size until not empty or infeasible region is found
Rn= [];
k = 0;
if ~builtin('isempty',Radj) 
    while all(Radj.Internal.I==R.Internal.I) || check_overlaps(Radj,R)
       % if Radj is the same as R (or possibly overlaps it)
       thn = thn + g'*MPTOPTIONS.region_tol/2;
       [Radj, piv] = lcp_getRegionFromTheta(thn, lc, R.Internal.I, HASH, regions, UNEX);
       % add pivots
       R.setInternal('piv', R.Internal.piv + piv);
       if builtin('isempty',Radj)
           % did not find feasible region, possible boundary, return
           if ~ISPARALLEL
               update_boundary(R, i);
           end
           return;
       end
       k = k+1;
       if k>MPTOPTIONS.modules.solvers.plcp.maxsteps
           %safety break
           break
       end
    end
    % if the recomputed Radj is not empty, go to the end
    if isEmptySet(Radj)
        % possible boundary
        if ~ISPARALLEL
            update_boundary(R, i);
        end
        
        % new step
        tn = shoot(Radj,g,thn);
        if tn.exitflag==MPTOPTIONS.OK
            thn = tn.x;
            
            % new region
            [Rd, piv] = lcp_getRegionFromTheta(thn, lc, R.Internal.I, HASH, regions, UNEX);
            % add pivots
            R.setInternal('piv', R.Internal.piv + piv);
            
            if builtin('isempty',Rd)
                Rd = Radj;
            end
            % collect regions into an array
            Rn = [Radj, Rd];
            
            % mark if the region is not the same
            Rni = intersect(Rn(end),Rn(end-1));
            ts = ~Rni.isEmptySet;
            %ts = isfeasible(Rni);
            k = 0;
            
            while isEmptySet(Rd) || ts
                if ts
                    % if the regions are equal, remove the last added one
                    thn = thn + g'*MPTOPTIONS.region_tol/2;
                    Rn(end) = [];
                else
                    % new step
                    tn = shoot(Rd,g,thn);
                    if tn.exitflag~=MPTOPTIONS.OK
                        break;
                    end
                    thn = tn.x;
                end
                % new region
                [Rd, piv] = lcp_getRegionFromTheta(thn, lc, R.Internal.I, HASH, regions, UNEX);
                % add pivots
                R.setInternal('piv',R.Internal.piv + piv);
                if builtin('isempty',Rd)
                    break;
                end
                % merge
                Rn = [Rn, Rd];
                
                % check if the regions are not the same
                % this is more strict comparison than overlap because testing
                % overlaps works with region_tol -which does not apply here
                Rni = intersect(Rn(end),Rn(end-1));
                %ts = isfeasible(Rni);
                ts = ~Rni.isEmptySet;

                k = k+1;
                % safety break
                if k>MPTOPTIONS.modules.solvers.plcp.maxsteps
                    break;
                end
            end
        end
    end
end

if numel(Rn)>2
    % take the last nonempty discovered region
    Radj = Rn(end);
end

% We need to check if the discovered neighbor Radj does not overlap R. If
% yes, we add the sign-reversed inequality corresponding of the current
% facet to Radj in order to cut the overlapping part.
% 
if ~builtin('isempty',Radj) && check_overlaps(Radj, R)
   P = Polyhedron('H',[Radj.H; -g, -h] );
   
   % make P irredundant and store data
   P.minHRep();
   
   I=Radj.Internal.I;
   %x=Radj.Internal.x;
   piv = Radj.Internal.piv;
   redundant_rows = Radj.Internal.redundant_rows;
   th = Radj.Internal.th;
   if P.isEmptySet
       % remainder after cutting is empty - Radj is empty
       Radj = [];
   else
       % restore Radj with the cut       
       Radj = P;
       % include new index set
       Radj.setInternal('I', I);       
       %Radj.setInternal('x', x);
       Radj.setInternal('piv',piv);
       Radj.setInternal('redundant_rows',redundant_rows);
       Radj.setInternal('th',th);
   end
end

end

function [R, piv] = lcp_getRegionFromTheta(thn, lc, I, HASH, regions, UNEX)
%
% for fixed value of theta find a corresponding region by solving LCP
%

global MPTOPTIONS

% if theta is out of feasible bounds, return empty region
if any(lc.Ht*[thn;-1] > MPTOPTIONS.abs_tol)
    R = [];
    piv = 0;
    if MPTOPTIONS.verbose >= 2 || MPTOPTIONS.modules.solvers.plcp.debug >= 1
        %fprintf('Warning: No facet constraints are active!... continuing\n')
        fprintf('Boundary facet (detected by fixed step).\n');
    end
    return;
end

% start from feasible basis I
try
    [z,w,basis,flag,piv] = lcpm(lc.M,lc.q+lc.Q*thn, I);
    if flag~=MPTOPTIONS.OK
        % solve non-parametric LCP with new parameter
        [z,w,basis,flag,piv] = lcp(lc.M,lc.q+lc.Q*thn, MPTOPTIONS.modules.solvers.lcp);
    end
catch
    % if matlab version of lex-lcp fails, retry with mex version
    [z,w,basis,flag,piv] = lcp(lc.M,lc.q+lc.Q*thn, MPTOPTIONS.modules.solvers.lcp);
end

% return region only if feasible
if flag==MPTOPTIONS.OK
    % sort and transpose
    basis = sort(basis(:))';
    
    % if the region cannot be computed from basis, mark as empty
    if rank(lc.A(:,basis),MPTOPTIONS.abs_tol)~=lc.n
        R = [];
    else                
        % construct region from basis
        R = lcp_getRegion(lc, basis, HASH, regions, UNEX, piv);
        % store theta based on which the region was determined
        if ~builtin('isempty',R)
            R.setInternal('th', thn);       
        end
    end
else
    R = [];
end



end

function Radj = lcp_ComputeAdj(R,i,lc, HASH, regions, UNEX, ISPARALLEL)
%
% computes adjacent region along the given facet via pivoting from region R
%
% inputs: R - current region
%         i - index of a facet to step over
%         lc - structure with LCP data (M, q, Q, Ath, bth)
%         HASH - hash table with optimal bases defining polyhedra        
%         regions - an array of already discovered regions
%         UNEX - an array of unexplored regions
%
% outputs:  Radj - adjacent regions

global MPTOPTIONS

% extract information from the face equation (normal to the facet)
H = R.H;
gamma = H(i,1:end-1)';
b = H(i,end);

% normalize facet equation gamma'*th = b, s.t. ||gamma||_2 = 1; (optional)
b = b / norm(gamma); gamma = gamma / norm(gamma);

% Embed theta in the facet (eliminate the facet equation).
% Eliminate theta that corresponds to maximum gamma because
% there will be no problem when dividing with gmax.
% If there are multiple maxima, we take the first one (done automatically
% by function max).
[~,p] = max(abs(gamma));
gmax = gamma(p);

% Elimination proceeds as follows:
% g1*th1 + g2*th2 + ... + gmax*thp + ... + gd*thd = b
%
% extract thp:
% thp = 1/gmax * (b - g1*th1 -g2*th2 - ... - gd*thd) (except gmax*thp)
% thp = - 1/gmax * (g1*th1 +g2*th2 + ... + gd*thd) + b/gmax (except gmax*thp)
%
% put thp inside Q*th:
% Q*th = [q11*th1 + q12*th2 + ... + q1p*thp + ... + q1d*thd; ...]
%      = [(q11-q12*g1/gmax)*th1 + ... + (q1d-q12*gd/gmax)*thd + q12*b/gmax; ...] (except thp)
%
% From the above equation we extract matrices such that Qnew*thf + qnew 
% thf is theta except thp, Qnew corresponds to columns of Q that multiply thf
% qnew is extracted and added to remaining term q
% Qnew = Q(:,np) - Q(:,p)*gamma(np)'/gmax
% qnew = Q(:,p)*b/gmax + q

np    = true(lc.d,1);
np(p) = false; %indices of all except gp

% new LCP problem, including variable alpha for step size
% w - M*z  = qnew + Qnew*thf + Q*gamma*alpha  + eps
% [I -M -Q*gamma]*[w;z;alpha] = qnew + Qnew*thf + eps

% obtain matrices Qnew, qnew
if lc.d==1
    % for dimension 1 we reduce to non-parametric LCP
    lcnew.Q = zeros(lc.n,lc.d);
else    
    lcnew.Q = lc.Q(:,np) -lc.Q(:,p)*gamma(np)'/gmax;
end
lcnew.q = lc.Q(:,p)*b/gmax + lc.q;

% Augment A with the parametric column
% variables are [w;z;alpha] where alpha is the step size along this facet
lcnew.A = [lc.A -lc.Q*gamma];

% problem dimension remains the same
lcnew.n = lc.n;
% parameters thf are in dimension d-1
lcnew.d = lc.d-1;
% store old LCP form 
lcnew.lcold = lc;
% store indices of thp
lcnew.p = p;
lcnew.np = np;

% Compute H-description of the facet in the reduced dimension (d-1)
I = R.Internal.I; 

% use QR factorization when required
if MPTOPTIONS.modules.solvers.plcp.QRfactor
    %factorize Anew with QR and store factors to f
    [f.Q,f.R] = qr(lcnew.A(:,I));
    %store index set
    f.I = I;
    % solve upper triangular system based on Q, R factors
    Hth=QRsolve(lcnew.A,[-lcnew.Q,lcnew.q],I,f);
else
    f=[];
    Hth = lcnew.A(:,I)\[-lcnew.Q lcnew.q];
end

n = sqrt(sum(Hth.*Hth,2));
Hth(n < MPTOPTIONS.zero_tol,:) = []; % eliminate zero rows

% Bring gamma into the basis (initial entering variable)
enter = 2*lc.n+1;

[Adj, adj_flag] = lcp_pivot(I, enter, lcnew, -1, Hth, [], 0, f, HASH, regions, UNEX);
% check whether the facet is really boundary
if ~ISPARALLEL && adj_flag>0
    update_boundary(R, i);
end

Radj = cell(1,length(Adj));
for j=1:length(Adj)
    % must sort
    Adj(j).I = sort(Adj(j).I);
    %piv = piv + Adj(j).piv;
    
    % get the region for given basis I
    Radj{j} = lcp_getRegion(lc,Adj(j).I,HASH,regions,UNEX,Adj(j).piv);
end
% kick out empty regions
%Radj(cellfun(@isEmptySet,Radj)) = [];
Radj = [Radj{:}];

if (MPTOPTIONS.modules.solvers.plcp.debug >=1)
    if length(Adj) > 1
        disp('Adjacent region does not have facet-to-facet property or tiny regions have been discovered!');
    end
end


end


%% MAIN PIVOTING FUNCTION
function [Adj, adj_flag] = lcp_pivot(I, enter, lc, gammaVar, Hth, Adj, piv, f, HASH, regions, UNEX)
%
% for given basis LC (determined by index set I), do a pivot step with 
% "enter" as an entering variable
%
%  inputs: I - index set for actual basis
%          enter - entering variable
%          lc - structure with LCP variables (M, q, Q, A, Ath, bth, d, n)
%          gammaVar - index of the variable "alpha" that is used to find
%          the feasible step size along the face
%          f- Q, R factors of A(:,I) corresponding to index set I
%         HASH - hash table with optimal bases defining polyhedra
%         regions - an array of already discovered regions
%
%          Hth - H representation of the facet (in dimension d-1)
%          Adj - structure with adjacent basis
%          piv - pivot counter
%
%  outputs: Adj: updated structure with adjacent basis
%           adj_flag: 0 if an adjacent region was found
%                     1 if the facet appears to be boundary

% pivoting new LCP problem, including variable alpha for step size
% w - M*z  = qnew + Qnew*thf + Q*gamma*alpha  + eps
% [I -M -Q*gamma]*[w;z;alpha] = qnew + Qnew*thf + eps
% thf belongs to polyhedron given by Hth

% Basic solution 
% xb = inv(A(:,B))*(qnew + Qnew*thf + eps -A(:,enter)*xe)

global MPTOPTIONS
%persistent pp
%if isempty(opt)
%    opt = Opt; 
%end

adj_flag = 0;

% limit the maximum number of pivots (if singularity appears or cycles)
if piv>MPTOPTIONS.modules.solvers.plcp.maxpivots
    return
end

% Basis
% [I -M -Q*gamma]*[w;z;alpha] = qnew + Qnew*thf + eps
% A = [I -M -gamma]

if MPTOPTIONS.modules.solvers.plcp.QRfactor
    % recursive factorization based on previous index set
    [CR, f] = QRsolve(lc.A, [lc.Q lc.q lc.A(:,enter)], I, f);
else    
    % direct LU factorization
    [L,U,p] = lu(lc.A(:,I),'vector');
    if any(abs(diag(U))<MPTOPTIONS.zero_tol)
        % if the basis of A forms a rank deficient matrix, quit pivoting
        return;
    end
    % use factorized solution to compute CR
    % CR = lc.A(:,I) \ [lc.Q lc.q lc.A(:,enter)];
    br = [lc.Q lc.q lc.A(:,enter)];
    CRl = linsolve(L,br(p,:),struct('LT',true));
    CR = linsolve(U,CRl,struct('UT',true));    
    %CRn = lc.A(:,I) \ [lc.Q lc.q lc.A(:,enter)];
end

% remove small entries
CR(abs(CR) < MPTOPTIONS.zero_tol) = 0;


% Entering variable sensitivity
v  = CR(:,end);
CR(:,end) = [];

% Compute rows that are strictly active
% changing this tolerance higher causes more lex-pivots
% if the value is too small, it might not find all adjacent regions
% but if the value is too big, generates overlapping regions
Z = sqrt(sum(CR.*CR,2)) < MPTOPTIONS.rel_tol;

% Compute the possible leaving vars
P = v > MPTOPTIONS.zero_tol;

if all(~P) && gammaVar>0
    if v(gammaVar)>-MPTOPTIONS.zero_tol
        if (MPTOPTIONS.verbose >= 2 ) || (MPTOPTIONS.modules.solvers.plcp.debug >= 1)
            fprintf('Boundary facet.\n');
        end
        adj_flag = 1;
        return;
    end
end

if all(~Z)
    if MPTOPTIONS.verbose >= 2 || MPTOPTIONS.modules.solvers.plcp.debug >= 1
        fprintf('Warning: No facet constraints are active!... continuing\n')
%        fprintf('Boundary facet.\n');
    end
    adj_flag = 2;
  return
end



%  if all(~P)
%      if MPTOPTIONS.verbose >= 2 || MPTOPTIONS.modules.solvers.plcp.debug >= 1
%          fprintf('Warning: Region is unbounded, skipping\n')
%      end
%    return
%  end

% Z intersect P
PZ = Z & P;


% Can only be optimal if gamma has already entered
if gammaVar ~= -1
    % Gamma can only increase to a positive value if v(gammavar) < 0
    % and if there is space for the entering variable to increase a s.p.
    % amount
    
    % terminal condition for pivoting
    if v(gammaVar) < -MPTOPTIONS.zero_tol % !!! must be compared against zero_tol or 0
        ts = false;

       if any(PZ)
            % test whether new basis forms a region that is adjacent to
            % current facet
            
            % region formed from new basis
            Inew = I;
            Inew(gammaVar) = enter;
            % don't do redundancy elimination here, because it is numerically
            % sensitive and costs time
            Radj = lcp_getRegion(lc.lcold, Inew, HASH, regions, UNEX, piv, false);  % region in full parameter dimension d              
            ts = false;
            % if new region is not empty, test if it is empty with "region_tol"
            if ~builtin('isempty',Radj)
                
                %if Radj.isEmptySet && Radj.isfeasible
                if Radj.isEmptySet || ~Radj.isFullDim || ~Radj.isBounded
                    ts = false;
                    % store regions indices
                    Adj(end+1).I = Radj.Internal.I;
                    %Adj(end).x = Radj.Internal.x;
                    Adj(end).piv = piv;
                else
                    ts = true;
                end
                               
%                 % sets initial statements on adjacency as false
%                 ts = false(size(Radj.A,1),1);
%                 if ~Radj.isEmptySet && size(Radj.A,1)>1
%                     % function isempty tests H-representation with
%                     % "region_tol" and if it the size of the region is less
%                     % than this tolerance, it is ignored
% 
% 
%                     if MPTOPTIONS.modules.solvers.plcp.debug >1
%                         % plot adjacent region
%                         plot(Radj,'color','m','LineWidth',2);
%                         disp('Press any key to continue...');
%                         pause
%                     end
%                     
%                     % for each facet of the new region test if it is adjacent
%                     for i=1:size(Radj.A,1)
%                         % obtain matrices in reduced dimension d-1 (eliminating thp)
%                         A = Radj.A;
%                         b = Radj.b;
%                         
%                         % facet equation f*th = r
%                         f = A(i,:);
%                         r=b(i);
%                         % remaining inequalities A*th <= b
%                         A(i,:) = [];
%                         b(i) = [];
%                         
%                         % f1*th1 + f2*th2 + ... + fp*thp + ...+fd*thd = r
%                         % thp = 1/fp*( r - f1*th1 - f2*th2 - ... - fd*thd) (except fp*thp)
%                         % thp = -1/fp * (f1*th1 + f2*th2 + ... fd*thd) + 1/fp*r
%                         %
%                         % put inside A*th
%                         % A*th = [a11*th1 + a12*th2 + ... + a1p*thp + ... + a1d*thd; ...]
%                         %      = [(a11-a1p*f1/fp)*th1 + ... + (a1d-a1p*g1/fp)*thd + a1p*r/fp; ...] (except thp)
%                         
%                         if norm(f(lc.p),1)<MPTOPTIONS.abs_tol
%                             % If fp is almost zero, we cannot divide by
%                             % this value. It means that thp can be
%                             % anything. We choose thp=0 which means
%                             % eliminating corresponding column.  
%                             Anew = A;
%                             Anew(:,lc.p) = [];
%                             bnew = b;
%                         else
%                             Anew = A(:,lc.np)-A(:,lc.p)*f(lc.np)/f(lc.p);
%                             bnew = A(:,lc.p)*r/f(lc.p) + b;
%                         end
%                         
%                         % solve primal LP for adjacency
%                         %       max t
%                         % s.t.: A*thf + t <=b  % actual facet space in dimension d-1
%                         %       Anew*thf + t <=bnew % new facet space
%                         Ao=Hth(:,1:end-1);
%                         bo=Hth(:,end);
%                         
%                         lp.A = [Ao ones(size(Ao,1),1); Anew ones(size(Anew,1),1)];
%                         lp.b = [bo;bnew];
%                         lp.f = zeros(lc.d+1,1);
%                         lp.f(end) = -1;
%                         rs=mpt_solve(lp);
%                         
%                         % if the objective value in solving the primal
%                         % problem is around zero, perform duality test
%                         if norm(rs.obj,Inf)<=MPTOPTIONS.abs_tol || rs.exitflag~=MPTOPTIONS.OK
%                             % solve lex dual LP
%                             %      min [[tmax;b;bnew] I ]'*y
%                             %  s.t.: [0...0 1]'    [0]
%                             %        [  A   1]     [ ]
%                             %        [ An   1]*y = [1]
%                             %                  y >= 0
%                             dp.Ae = [zeros(1,size(Ao,2)) 1;
%                                 Ao    ones(size(Ao,1),1);
%                                 Anew  ones(size(Anew,1),1)]';
%                             dp.be = zeros(size(dp.Ae,1),1);
%                             dp.be(end) = 1; % increasing this number may help too, for instance 1000 works fine
%                             dp.A = [];
%                             dp.b = [];
%                             dp.lb = zeros(size(dp.Ae,2),1);
%                             ff = [[MPTOPTIONS.infbound; bo; bnew] eye(length(bo)+length(bnew)+1)];
%                             for j=1:size(ff,2)
%                                 dp.f = ff(:,j); % pick j-th column
%                                 % solve dp
%                                 rd = mpt_solve(dp);
%                                 if rd.exitflag ==MPTOPTIONS.OK
%                                     if rd.obj > MPTOPTIONS.lex_tol
%                                         % LP is lex-positive, region is adjacent
%                                         ts(i) = true;
%                                         break;
%                                     elseif rd.obj < -MPTOPTIONS.lex_tol
%                                         % LP is not lex-positive, region is not
%                                         % adjacent
%                                         break;
%                                     elseif j==size(ff,2)
%                                         % this should not happen, we need
%                                         % to check again with 0 tolerance
%                                         for jn=1:size(ff,2)
%                                             dp.f = ff(:,jn); % pick j-th column
%                                             % solve dp
%                                             rd = mpt_solve(dp);
%                                             if rd.exitflag ==MPTOPTIONS.OK
%                                                 if rd.obj > 0
%                                                     % LP is lex-positive, region is adjacent
%                                                     ts(i) = true;
%                                                     break;
%                                                 elseif rd.obj < 0
%                                                     % LP is not lex-positive, region is not
%                                                     % adjacent
%                                                     break;
%                                                 elseif jn>=size(ff,2)
%                                                     disp('Lexicographic approach to adjacency problem failed due to LP solver problem. Continuing anyway ...');
%                                                 end
%                                             end                                        
%                                         end                                        
%                                     end
%                                 else
%                                     % report failore
%                                     disp('Error in solver when solving dual-lex adjacency LP. Continuing anyway ...');
%                                 end
%                                 
%                             end
%                         else
%                             ts(i) = (-rs.obj>MPTOPTIONS.abs_tol);
%                         end
%                         % break if any value is true (which means that Radj is adjacent region and we do not have to do more tests)                       
%                         if ts(i)
%                             break;
%                         end
%                     end
%                 end
            end

       end
        % terminating conditions
        if all(~PZ) || any(ts)
            I(gammaVar) = enter;
            Adj(end+1).I = I;
            if all(~PZ)
                Adj(end).result = 'optimal';
            else
                Adj(end).result = 'lex-optimal';
            end
            Adj(end).piv = piv;
            return;
        end
    end
end

if all(~PZ)
  % If there is only one possible pivot, then we don't have to do any more testing
  if nnz(P) == 1
    PZ = P;
  else
    % If all of the parameter-dependant columns are zero, then there's only one pivot
    zCR = abs(CR(P,1:end-1)) > MPTOPTIONS.zero_tol;
    if ~any(zCR(:))
      PZ = P;
    end
  end  
end

if all(~PZ)

    %   fprintf('possible complex pivot...\n');
    
    %   if 0
    %     % T := [CR(P,:) iB(P,:)] ./ v
    %     iB = inv(lc.A(:,I));
    %     T = [CR(P,:) iB(P,:)] ./ repmat(v(P),1,size(CR,2)+size(iB,2));
    %     CRT = T(:,1:lc.d);
    %
    %     % Group all possible leaving rows based on parallel
    %     % parameter-dependent normals
    %     nCRT = CRT ./ repmat(sqrt(sum(CRT.*CRT,2)),1,size(CRT,2));
    %     [B,uI,uJ] = unique_rows(nCRT,MPTOPTIONS.abs_tol);
    %     % Note: uI is the index into P => convert to index into CR
    %
    %     fP = find(P);
    %     % Check if there is a full-dimensional subset of the facet for each row
    %     % of B
    %
    %     LEAVE = [];
    %
    %     for i=1:size(B,1)
    %       r = uI(i); % Row of T
    %
    %       if size(B,1) > 1 % Only have to do complex test if there are choices
    %
    %         % Require that CRT(r,:)*[th;1] <= CRT(uI / r,:)*[th;1]
    %         t = repmat(CRT(r,:),length(uI)-1,1) - CRT([uI(1:i-1);uI(i+1:end)],:);
    %
    %         % Test if there's a full-dimensional region where the constraint t*[th;1] <= 0 is
    %         % satisfied
    %
    %         H = [t(:,1:end-1) -t(:,end);CR];
    %         d = sqrt(sum(H(:,1:end-1).*H(:,1:end-1),2));
    %         H(d < MPTOPTIONS.abs_tol, :) = []; % Remove zero rows
    %         d(d < MPTOPTIONS.abs_tol) = [];
    %
    %         % Check if H is full-dimensional
    %         %pp.clear;
    %         sol = mpt_solve(struct('A', [H(:,1:end-1) d;zeros(1,size(H,2)-1) 1],...
    %             'b',[H(:,end);1],...
    %             'f',[zeros(1,size(H,2)-1) -1]'));
    %
    %         %sol = pp.solve;
    % %         sol = opt.setData('A', [H(:,1:end-1) d;zeros(1,size(H,2)-1) 1], 'b', [H(:,end);1],...
    % %           'f', [zeros(1,size(H,2)-1) -1])).solve;
    %         if sol.exitflag ~= MPTOPTIONS.OK
    %             if MPTOPTIONS.verbose >= 2 || MPTOPTIONS.modules.solvers.plcp.debug >= 1
    %                 warning('Solver error... continuing...');
    %             end
    %           continue
    %         end
    %
    %         if -sol.obj < MPTOPTIONS.abs_tol * 1e3
    %           % Lower-dimensional pivot
    %           continue
    %         end
    %       end
    %
    %       % This is a valid pivot.
    %
    %       % Of those rows parallel to B(i,:), choose the one that's lexmin
    %       test = find(uJ == i);
    %       lmin = test(lexmin([nCRT(test,:) T(test,lc.d+1:end)]));
    %
    %       LEAVE = [LEAVE;fP(lmin)];
    %     end
    %
    %     for i=1:length(LEAVE)
    %       leave = LEAVE(i);
    %
    %       % This is a valid pivot
    %       Iout = I;
    %       if(gammaVar == -1),
    %         if MPTOPTIONS.verbose >= 2 || MPTOPTIONS.modules.solvers.plcp.debug >=1
    %           warning('gammaVar should not be empty here... continuing')
    %         end
    %         continue
    %       end
    %
    %       comp_leave = comp(Iout(leave),lc.n);
    %       Iout(leave)   = enter;
    %
    %       fprintf('COMPLEX PIVOT\n');
    %       Adj = lcp_pivot(Iout, comp_leave, lc, gammaVar, Hth, Adj, piv+1);
    %     end
    %   end
       
    if MPTOPTIONS.modules.solvers.plcp.debug>=1
        fprintf('Taking complex pivot.\n');
    end
    
    % number of possible leaving variables
    nnzP = nnz(P);
    % indices of possible leaving variables
    fP   = find(P);
    if MPTOPTIONS.modules.solvers.plcp.QRfactor
        iB = QRsolve(lc.A, eye(lc.n), I, f);
    else
        % beta = inv(A(:,I))
        iB   = lc.A(:,I) \ eye(lc.n);
    end
    
    % test each possible leaving variable for lexico-positivity
    for l = 1:nnzP
        
        % Gamma is defined as v(l)*beta(P,:) - v(P)*beta(l,:)
        G = v(fP(l)) * iB(P,:) - v(P) * iB(fP(l),:);
          
        % need to check if G*[q I] is lex positive while G*Q=0
        GQ = G*lc.Q;
        iszero = sqrt(sum(GQ.*GQ,2)) < MPTOPTIONS.zero_tol;
        if ~isMatrixLexPos([G(iszero,:)*lc.q G(iszero,:)])
            % matrix is not lex-positive, continue with next possible
            % variable
            continue
        end
        
        %  matrix is lex-positive, get the region
        if any(abs(GQ(:)) > MPTOPTIONS.zero_tol)
            %       fprintf('NON-ZERO GQ!\n');
    
            % get the region
            CR(:,1:end-1) = -CR(:,1:end-1);
            %       H = [[-GQ G*lc.q];Hth;CR];
            Hth = [[-GQ G*lc.q];CR];
            
            % prevent NaN or Inf
            if any(any(~isfinite(Hth)))
                continue;
            end
                        
            reg=Polyhedron(Hth(:, 1:end-1), Hth(:, end));
            if reg.isEmptySet || ~reg.isFullDim || ~reg.isBounded
                continue;
            end

%             % Remove zero rows
%             d = sqrt(sum(H(:,1:end-1).*H(:,1:end-1),2));
%             H(d < MPTOPTIONS.zero_tol, :) = []; 
%             d(d < MPTOPTIONS.zero_tol) = [];
%                          
%             % Check if H is full-dimensional by formulating a chebyshev
%             % center problem with the upper bound the radius 1
%             %    max r
%             % s.t. a_i*x + ||a_i||_2*r <= b for i=1,...,m
%             %          r <= 1
%             lp.A = [H(:,1:end-1) d;zeros(1,size(H,2)-1) 1];
%             lp.b = [H(:,end);1];
%             lp.f = [zeros(1,size(H,2)-1) -1];
%             sol = mpt_solve(lp);
% 
%             % if the direct approach failes, go thru polyhedron
%             if sol.exitflag ~= MPTOPTIONS.OK
%                 % we need to find other way to test full-dimensionality
%                 
%                 %       P = Polyhedron('H', H);
%                 %       x = P.interiorPoint;
%                 %       if ~x.isStrict
%                 %         continue % Not full-dim
%                 %       end
%                 
%                 reg = Polyhedron('H',H);                                
%                 if ~reg.isFullDim
%                     continue
%                 end
%             end
%             
%             R = -sol.obj;
%             % if the radius is lower than the smallest available tolerance
%             % for empty region, continue
%             if R < MPTOPTIONS.region_tol/2
%                 % Lower-dimensional pivot
%                 continue
%             end
%             
%             % don't do redundancy elimination here, because it is numerically
%             % sensitive and costs time
%             % Hth = Polyhedron('H', H).minHRep.H;
%             Hth = H;
        end
        
        % This is a valid pivot
        Iout = I;
        leave = fP(l);
        if(gammaVar == -1),
            if MPTOPTIONS.modules.solvers.plcp.debug>=1
                disp('gammaVar should not be empty here...')
            end
            continue
        end
        comp_leave = comp(Iout(leave),lc.n);
        Iout(leave)   = enter;
        

        [Adj, adj_flag] = lcp_pivot(Iout, comp_leave, lc, gammaVar, Hth, Adj, piv+1, f, HASH, regions, UNEX);
    end;
else
    % Pivot is unique and not a function of theta
    fPZ = find(PZ);
    if length(fPZ) > 1
        %iB = inv(lc.A(:,I));
        if MPTOPTIONS.modules.solvers.plcp.QRfactor
            iB = QRsolve(lc.A, eye(lc.n), I, f);
        else
            % beta = inv(A(:,I))
            iB   = lc.A(:,I) \ eye(lc.n);
        end
        %%%% CONSIDER THE -ALPHA TERM...
        %     fprintf('nnz(PZ) = %i\n', nnz(PZ));
        %     if nnz(PZ) > 4
        %       fprintf('Found one...\n');
        %     end
        leave = fPZ(lexmin([CR(PZ,end) iB(PZ,:)]./repmat(v(PZ),1,lc.n+1)));
    else
        leave = fPZ;
    end
    if gammaVar == -1, gammaVar = leave; end;
    %iLeave = I(leave);
    comp_leave = comp(I(leave),lc.n);
    I(leave)   = enter;
    
    %   % If this is an LP, then the next pivot will be a dual/primal pivot
    %   % Classify the leaving variable
    %   if iLeave < lc.numVars,           lType = 'Dual Slack';
    %   elseif iLeave < lc.n,             lType = 'Primal Slack';
    %   elseif iLeave < lc.n+lc.numVars, lType = 'Primal Var';
    %   elseif iLeave < 2*lc.n,           lType = 'Dual Var';
    %   else                               lType = 'gamma';
    %   end
    %
    %   if enter < lc.numVars,            eType = 'Dual Slack';
    %   elseif enter < lc.n,              eType = 'Primal Slack';
    %   elseif enter < lc.n+lc.numVars,  eType = 'Primal Var';
    %   elseif enter < 2*lc.n,            eType = 'Dual Var';
    %   else                               eType = 'gamma';
    %   end
    %
    %   fprintf('Pivot : enter %2i leave %2i | enter / leave : %s / %s \n', enter, iLeave, eType, lType);
    
    [Adj, adj_flag] = lcp_pivot(I, comp_leave, lc, gammaVar, Hth, Adj, piv+1, f, HASH, regions, UNEX);
end
end

function tf = isMatrixLexPos(A)

global MPTOPTIONS


P = cumsum(A >  MPTOPTIONS.lex_tol,2) > 0;
N = cumsum(A < -MPTOPTIONS.lex_tol,2) > 0;

tf = ~any(P(:)-N(:) < 0);

% I = ones(size(A,1),1);
%
% tf = true;
% for i=1:size(A,2)
%   if any(N(I,i)), tf = false; return; end
%   I(P) = false;
% end

% tf = true;
% for i=[1:size(A,1)]
%   if(lexpos(A(i,:)) == -1)
%     tf = false;
%     break
%   end
% end
%
% if tf ~= tfx
%   error('New lex method doesn''t work')
% end

end

function m = lexmin(t)
%
% Compute the row of t that is the lexicographic minimum
%

global MPTOPTIONS

if size(t,1) == 1
  m = 1;
  return
end

I = true(size(t,1),1);
J = find(sum(abs(t),1) > MPTOPTIONS.zero_tol); % non-zero columns
for j=1:length(J)
  i = J(j);
  % Remove from consideration those that are not minimal
  I(t(:,i) > min(t(I,i)) + MPTOPTIONS.lex_tol) = false;
  if nnz(I) == 1
    m = find(I);
    break
  end
end
if nnz(I) > 1
    % retry with absolutely 0 tolerance
    J = find(sum(abs(t),1) > 0); % non-zero columns
    for j=1:length(J)
        i = J(j);
        % Remove from consideration those that are not minimal
        I(t(:,i) > min(t(I,i))) = false;
        if nnz(I) == 1
            m = find(I);
            break
        end
    end

    if nnz(I)>1
        error('Lexmin is not unique!')
    end
end


% m = 1;
% for i=[2:size(t,1)]
%   if(lexgreater(t(m,:),t(i,:)) == 1)
%     m = i;
%   end;
% end;
%
% if mx ~= m
%   error('New lexmin doesn''t work')
% end

end

function tf = lexgreater(a,b)
%
% returns:
%    1: a > b
%   -1: a < b
%    0: a = b
%

tf = lexpos(a-b);
end

function tf = lexpos(a)
%
% returns:
%   1: a > 0
%  -1: a < 0
%   0: a = 0

global MPTOPTIONS

% if(all(abs(a) < zerotol))
%   tf = 0;
%   return;
% end;

% if(m == 1 || n == 1)
tf = 0;
for i=1:length(a)
  if (a(i) >  MPTOPTIONS.lex_tol)
      tf =  1;
      break;
  end;
  if (a(i) < -MPTOPTIONS.lex_tol)
      tf = -1;
      break;
  end;
end;
% else
%   tf = 1;
%   for i=[1:m]
%     if(lexpos(a(i,:)) == -1)
%       tf = -1;
%       break;
%     end;
%   end;
% end;
end

function c = comp(t,m)
% returns an index of a complementary variable to "t" working with
% dimension "m"
%
% if t is a vector, then this function returns complementary variable to
% each index in the vector

c = zeros(1,length(t));
for i=1:length(t)
    if (t(i) > m)
        c(i) = t(i)-m;
    else
        c(i) = t(i)+m;
    end
end
end


function [X,fnew] = QRsolve( A, B, J, f)
% solution of equation (Q*R)*X=B with the help of recursive QR
% factorization based on the current index set J and previous index set
% Jold
% inputs: A, B from  A*X=B matrix equation
%         J -new index set
%         structure f containting Q, R from the previous factorization
%         A(:,Jold) = Q*R given by Jold

Q = f.Q;
R = f.R;
Jold = f.I;

% which column has changed
v = (J~=Jold);

if any(v)
    % v must be column
    v=v(:);
    % pick the appropriate column of A-Aold
    deltaA = A(:,J)-A(:,Jold);
    u = deltaA(:,v);
    % update Q,R, factors
    [Q,R] = qrupdate(Q,R,u,double(v));
end
    
% solve upper triangular system
X = linsolve(R,Q'*B,struct('UT',true));

% return updated Q,R factors
fnew.Q = Q;
fnew.R = R;
fnew.I = J;

end

function tf = islexposbasis(iB,b)
%
% internal function used LCPM 
%

n = length(b);
D = [iB*b iB];
%D = iB*[b A];
for i=1:n
  if(lexpos(D(i,:)) < 1)
    tf = 0;
    return;
  end;
end;
tf = 1;
end


function [z,w,I,result,piv] = lcpm(M,q,Io,art,term)
% matlab implementation of standard lexicographic Lemke's algorithm
%
% [z,w,I,result,piv] = lcp(M,q,Io,art,term)
%
% Finds a solution (w,z) s.t.
%  w - Mz = q, w>=0, z>=0, w'z=0
%
% Io is an ACB (or empty) and art is the column for the artificial variable
%
% This is a MATLAB alternative to LCP function, always use LCP if possible.
%

global MPTOPTIONS

% some preliminary checks
validate_realmatrix(M);
validate_realvector(q);

if nargin<2
    error('At least 2 arguments "M, q" must be provided');
end

piv = 0;

[m,n] = size(M);

% make q column vector
q = q(:);

if(nargin < 3)
    A = [eye(m) -M -ones(m,1)];
    % Initial basis
    I = 1:m;
    term = 0;
elseif nargin < 4
    A = [eye(m) -M -ones(m,1)];
    if ~isempty(Io)
        I = Io;
    else
        I = 1:m;
    end
    term = 0;
else
    A = [eye(m) -M art];
    I = Io;
end;

% Bring in the artificial variable
iB = A(:,I) \ eye(m);
qb = iB*q;

% If q is already positive then this is a feasible basis
if(min(qb) >= 0)
    w = q;
    z = zeros(m,1);
    result = 1;
    return;
end;

% Choose the leaving variable
canleave = find(abs(qb-min(qb)) < MPTOPTIONS.zero_tol);
D = iB*[q eye(n)];
D = D(canleave,:);
leave = lexmin(D);
leave = canleave(leave);

% Bring the artificial variable in
enter = size(A,2);
% and qmin out
didleave = I(leave);

I(leave) = enter;
nI = setdiff(1:2*m,I);

% If we're not lex-positive then modify the problem to be lex-positive
iB = A(:,I) \ eye(m);
if(islexposbasis(iB,q) == 0)
    %     if (MPTOPTIONS.verbose >= 2 ) || (MPTOPTIONS.modules.solvers.plcp.debug >= 1)
    %         disp('Changing input to make initial basis lex-positive.');
    %     end
    A = iB*A;
    q = iB*q;
end;

piv = 0;

while (1)
    piv = piv + 1;
    
    B  = A(:,I);
    N  = A(:,nI);
    iB = B \ eye(m);
    
    % The artificial variable has left - this is feasible
    if (didleave == 2*m+1)
        result = 1;
        break;
    end;
    
    % The entering variable must the complement of the previous leaving variable
    if(didleave > m)
        enter = didleave-m;
    else
        enter = didleave+m;
    end;
    enter = find(nI==enter);
    
    % Choose the leaving variable by the lex-min rule
    v = iB*N(:,enter);
    canleave = find(v > MPTOPTIONS.zero_tol);
    if(isempty(canleave))
        result = -2; % Unbounded (i think...)
        break;
    end;
    
    % Test for special termination criteria
    if(term == 1)
        % Test if all of qzero is in canleave
        qzero = sort(find(iB*q < MPTOPTIONS.zero_tol));
        done = 1;
        
        for i=1:length(qzero)
            if(isempty(find(canleave==qzero(i), 1)))
                done = 0;
                break;
            end;
        end;
        
        if(done == 1)
            I = [I nI(enter)];
            I = setdiff(I,2*m+1);
            result = 1;
            break;
        end;
    end;
    
    n = size(B,1);
    D = iB*[q eye(n)];
    D = D(canleave,:)./repmat(v(canleave),1,size(D,2));
    leave = lexmin(D);
    leave = canleave(leave);
    
    %  fprintf('I = '); fprintf('%2i ',I);
    %  fprintf('enter = %2i leave = %2i\n',nI(enter),I(leave));
    
    % Make the pivot
    t = I(leave);
    I(leave) = nI(enter);
    nI(enter) = t;
    
    didleave = t;
    
    %I  = I(randperm(length(I)));%sort(I);
    %nI = sort(nI);
end;

x = zeros(2*m,1);
x(I) = iB*q;
w = x(1:m);
z = x(m+1:end);

end

function adj_list = verify_graph(regions,adj_list)
%
% check the adjacency list and correct if neighbors don't fit
%

global MPTOPTIONS
if isempty(MPTOPTIONS)
   MPTOPTIONS = mptopt;
end

validate_polyhedron(regions);
narginchk(2, 2);

if ~isa(adj_list,'cell')
    error('Provided adjacency list must be in a CELL format.');
end

tic
for i=1:length(regions)
    
    if toc > MPTOPTIONS.report_period && MPTOPTIONS.verbose>=0
        fprintf('progress: %d/%d\n', i, length(regions));
        tic;
    end
    
    % clear regions outside of the range first
    for j=1:length(adj_list{i})
        if ~isempty(adj_list{i}{j})
            adj_list{i}{j}(adj_list{i}{j}>length(regions)) = [];
        end
    end
    
    % check each facet
    n = length(adj_list{i});
    for j=1:n % loops thru facets of region i
               
        % get all neighbors to this facet
        lt = adj_list{i}{j};
        
        % if there is no neighbor in the list, do a hard check via
        % sequential search
        if isempty(lt)
            % find facet interior point
            xc = regions(i).chebyCenter(j);
            if xc.exitflag==MPTOPTIONS.OK
                % small step beyond facet
                g = regions(i).H(j,1:end-1);
                ng = norm(g);
                if ng>MPTOPTIONS.abs_tol
					th = xc.x + 1.3/ng*g'*MPTOPTIONS.region_tol;
					[isin, inwhich]=regions.isInside(th);
					if isin
                        % the point could be located in the current region,
                        % we must remove this region from the list
						lt = setdiff(inwhich,i);
						adj_list{i}{j} = lt;
					end
                end
            end
        end
        
        % for each neighbor check if it lists region i
        for k=1:numel(lt)
            % get all possible neighbors of region k
            bt = [adj_list{lt(k)}{:}];            
            % k does not list i, correct
            if ~ismember(i,bt)
                P = regions(i);
                Q = regions(lt(k));    
                % must try if P and Q are neighbors first
                [ts,iP,iQ] = P.isNeighbor(Q,j);
                % if the neighborhood fails try distance with some
                % tolerance
                if ~ts
                    facetP = P.getFacet(j);
                    xP = facetP.chebyCenter.x;
                    if isempty(xP)
                        % if the solver fails, normalize and try again
                        facetP.normalize;
                        xP = facetP.chebyCenter.x;
                        if isempty(xP)
                            % if the solver fails second time, solve
                            % feasibility problem
                            xP = Opt(facetP).solve.xopt;
                        end
                        if isempty(xP)
                            xP = NaN(facetP.Dim,1);
                        end
                    end
                    dt = Q.project(xP);
                    ts = dt.dist<MPTOPTIONS.rel_tol;
                    if isempty(ts)
                        ts = 0;
                    end
                    if ts
                        % find the facet of Q that is the closest to P
                        % distance of the point xP to the polytope Q
                        dtQ = zeros(size(Q.H,1),1);
                        for jj=1:size(Q.H,1)
                            % compute the distance to center of each facet of Q
                            facetQ = Q.getFacet(jj);
                            xQ = facetQ.chebyCenter.x;
                            if isempty(xQ),
                                facetQ.normalize;
                                xQ = facetQ.chebyCenter.x;
                                if isempty(xQ)
                                    % if the solver fails second time, solve
                                    % feasibility problem
                                    xQ = Opt(facetQ).solve.xopt;
                                end
                                if isempty(xQ)
                                    % all approaches failed, leave
                                    xQ = NaN(facetQ.Dim,1);
                                end
                            end
                            dtQ(jj) = norm(xP-xQ);
                        end                        
                        [mindtQ,iQ] = min(dtQ);
                    end
                end
                if ts
                    % regions are neighbors, correct
                    adj_list{lt(k)}{iQ} = [adj_list{lt(k)}{iQ}, i];
                else
                    % regions are not neighbors, remove region lt(k) from i
                    adj_list{i}{j}(adj_list{i}{j}==lt(k))=[];
                end
            end
        end
        
    end

end


% plot(regions);
% for i=1:length(regions)
%     xc = chebyCenter(regions(i));
%     text(xc.x(1),xc.x(2),num2str(i));
% end


end

function F=feasible_set(regions, adj, problem)
%
% construct feasible set from the optimizer using the adjacency list
%

% global options
global MPTOPTIONS

H = [];
for i = 1:length(adj)
	% faces over which there is no neighbor => faces of the feasible set
	boundary_idx = find(cellfun('isempty', adj{i}));
	
	% simple code, but heavily relies on correctness of the adjacency list
	if ~isempty(boundary_idx)
		H = [H; regions(i).H(boundary_idx, :)];
	end

end

F = Polyhedron(H(:, 1:end-1), H(:, end));

if MPTOPTIONS.modules.solvers.plcp.adjcheck
    if isEmptySet(F)
        warning('mpt_plcp:adjacencyList','Numerical problems in verifying the adjacency list. The adjacency list may be wrong. Trying a failback method. This could take a while.');
        % The feasible set is wrong. Terminate this approach and try
        % more robust way of detecting the feasible set.
        
        % fallback scenario
        H = [];
        for i = 1:length(adj)
            % faces over which there is no neighbor => faces of the feasible set
            boundary_idx = find(cellfun('isempty', adj{i}));
            
            if ~isempty(boundary_idx)
                for j=1:numel(boundary_idx)
                    % boundary facet
                    bfacet = regions(i).getFacet(boundary_idx(j));
                    xF = bfacet.chebyCenter.x;
                    if ~isempty(xF)
                        % check that a point beyond the boundary gives infeasible PLCP
                        xplus = xF + bfacet.Ae'/norm(bfacet.Ae)*MPTOPTIONS.rel_tol;
                        [~,~,~,exitflag]=lcp(problem.M,problem.q+problem.Q*xplus);
                        if exitflag~=MPTOPTIONS.OK
                            % check that a point inside is contained in the
                            % same region (to make sure that this is actually
                            % the facet of this region)
                            xminus = xF - bfacet.Ae'/norm(bfacet.Ae)*MPTOPTIONS.rel_tol;
                            if regions(i).contains(xminus)
                                H = [H; regions(i).H(boundary_idx(j),:)];
                                P = Polyhedron('H',H);
                                if isEmptySet(P)
                                    % The feasible set is wrong. Show the
                                    % warning and quit
                                    warning('mpt_plcp:feasibleSet','Numerical problems in computing the feasible set. Failback approach failed.');
                                    return;
                                end
                            end
                        end
                    end
                end
            end
            
        end
        F = Polyhedron(H(:, 1:end-1), H(:, end));
        
    end
end

end


%-----------------------------------------------------------
function [MAP, rewrite] = map_put(MAP, key, new_value)
% insers into HASH "value" under "key"
% returns the updated map as the first input
% returns rewrite=true if the key existed before, false otherwise
%
% implicitly assumes "key" is a vector of integers! (because int2str is so
% much faster than num2str)

% Note: containers.Map does not offer the option to retrieve the index of a
% key. In fact, keys are sorted alphabetically in containers.Map. Since
% adjacency list depends on indices, we store the value as a structure:
%   stored_value.value = user-provided value
%   stored_value.index = automatically incremented index of the key

key = int2str(key);
rewrite = MAP.isKey(key);
if rewrite
	% just update the stored value, do not change its index!
	old = MAP(key);
	old.value = new_value;
	MAP(key) = old;
else
	% insert new value, update its index
	MAP(key) = struct('value', new_value, 'index', MAP.Count+1);
end

end

%-----------------------------------------------------------
function idx = map_getIndex(MAP, key)
% returns index of a key. idx=[] if the key does not exist in the map
%
% implicitly assumes "key" is a vector of integers! (because int2str is so
% much faster than num2str)

key = int2str(key);
if MAP.isKey(key)
	value = MAP(key);
	idx = double(value.index); % make sure the index is a double
else
	idx = [];
end

end

%-----------------------------------------------------------
function update_boundary(region, index)
% updates the stack of boundary half-spaces

global BNDH

use_caching = true;

% Before calling is_boundary_facet() we check whether the half-space
% is already in BNDH. If it is, we can exit quickly. Note that such a check
% requires the half-space to be normalized, e.g. to h'*x<=1
%
% An another reason is to reduce the number of (possibly redundant)
% boundary half-spaces.

% half-space to check
hs = region.H(index, :);

if use_caching
    % normalize the candidate half-space to h'*x<=+/-1
    % avoid dividing by zero
    if abs(hs(end))>1e-6
        % note that we must divide by abs(hs(end)) as not to reverse the
        % inequality if hs(end) is negative
        hs = hs./abs(hs(end));
    end
    
    % is the half-space already included? (helps to reduce the number of LCPs
    % to be solved)
    included = any(sum(abs(BNDH-repmat(hs, size(BNDH, 1), 1)), 2)<1e-8);
else
    included = false;
end

if ~included
    % check the boundary status by solving an LCP
    if is_boundary_facet(region, index)
        BNDH = [BNDH; hs];
    end
end

end


%-----------------------------------------------------------
function isboundary = is_boundary_facet(region, index)
% checks whether a given facet is indeed a boundary of the feasible set
%
% A facet is boundary if:
% * the LCP problem is infeasible when solved for a point accross the facet
% * the LCP problem is feasible when solved for an interior point of the
%   region, close to the investigated facet

global MPTOPTIONS INPUT_LCP

% if enabled, we also check whether the problem is feasible on an inward
% point of the investigated facet
check_inward = false;

% compute a point on the facet
facet = region.chebyCenter(index);
if facet.exitflag == MPTOPTIONS.OK
    % compute the point accross the j-th facet
    step_size = MPTOPTIONS.rel_tol*10;
    dx = region.A(index, :)'/norm(region.A(index, :)')*step_size;
    % solve the problem for the new point
    lcpsol = INPUT_LCP.solve(facet.x+dx);
    % infeasible = candidate boundary facet, feasible = not boundary
    isboundary = (lcpsol.exitflag ~= MPTOPTIONS.OK);
    if check_inward && isboundary
        % check the inwards direction, must be feasible, otherwise we have
        % a fake boundary
        lcpsol = INPUT_LCP.solve(facet.x-dx);
        
        % infeasible outwards, feasible inwards => boundary facet
        % infeasible outwards, infeasible inwards => not boundary
        isboundary = (lcpsol.exitflag==MPTOPTIONS.OK);
    end
else
    % no point on the facet found = not boundary
    isboundary = false;
end

end % function

