function [Pn,Fi,Gi,activeConstraints,Phard,details]=mpt_mpqp(Matrices,Options)
%MPT_MPQP Explicitly solves the given quadratic program (QP)
%
% [Pn,Fi,Gi,activeConstraints,Phard,details]=mpt_mpqp(Matrices,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Solves the following QP as a multiparametric program:
% min_U  0.5 U' H U + (x' F + Cf) U + x' Y x + Cx x + Cc
% subj. to  GU <= W + E x    (constraints)
%        bndA*x<= bndb       (bound exploration space)
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% Matrices - a struct with all the parameters which are needed.
%             See description above for explanation.
%    Matrices.G=G;
%    Matrices.E=E;             
%    Matrices.W=W;
%    Matrices.H=H;
%    Matrices.F=F;
%    Matrices.Y=Y;
%    Matrices.Cf=Cf;
%    Matrices.Cx=Cx;
%    Matrices.Cc=Cc;
%    Matrices.bndA=bndA;
%    Matrices.bndb=bndb;
%
% Options.verbose     - Optional: level of verbosity
% Options.lpsolver    - Optional: which LP solver to use (help mpt_solveLP)
% Options.qpsolver    - Optional: which QP solver to use (help mpt_solveQP)
% Options.step_size   - Optional: length of step over a facet; Making this
%                          value larger often mitigates numerical problems but
%                          may produce small gaps in the partition.
%                          Default is 1e-4;
%
% Options.debug_level  
%           Due to numerical problems tiny regions are sometimes difficult to
%           calculate, i.e. are not identified at all. This may create "gaps"
%           in the computed control law. For the exploration, these will be
%           jumped over and the exploration in the state space will continue.
%           "debug_level" can have three values:
%       
%           0: No debug done
%           1: A tolerance is given to find gap in the region partition,
%              small empty regions inside the region partition will be discarded.
%              Note that this is generally not a problem, since the feedback law 
%              is continuous and can therefore be interpolated easily.
%              Correction to the calculation of the outer hull.
%           2: Zero tolerance to find gap in the region partition, empty regions
%              if they exist, will be detected, i.e. the user will be notified.
%              Correction to the calculation of the outer hull.
%
% Options.zero_tol    - everything below this value is considered zero
%                       (default is 1e-13)
%
% Note: If Options is missing or some of the fiels are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% Pn,Fi,Gi           - for region Pn(i) the optimal input is U=Fi{i}*x+Gi{i} 
% activeConstraints  - Cell Array which stores the active constraints 
%                      of the optimizer in each region.
% Phard              - The set of feasible states (i.e. union of
%                      all regions) as Phard.H*x<=Phard.K
% details            - structure with fields Ai, Bi, Ci in which the value
%                      function associated to each region is stored
%                      V = x' Ai{i} x + Bi{i} x + Ci{i}
%
% see also MPT_OPTCONTROL, MPT_CONSTRUCTMATRICES, MPT_MPLP

% Copyright is with the following author(s):
%
% (C) 2007 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%     kvasnica@control.ee.ethz.ch
% (C) 2005 Johan Loefberg, Automatic Control Laboratory, ETH Zurich,
%     loefberg@control.ee.ethz.ch
% (C) 2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%     kvasnica@control.ee.ethz.ch
% (C) 2004 Mato Baotic, Automatic Control Laboratory, ETH Zurich,
%     baotic@control.ee.ethz.ch
% (C) 2003 Pascal Grieder, Automatic Control Laboratory, ETH Zurich,
%     grieder@control.ee.ethz.ch
% (C) 2003 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%     kvasnica@control.ee.ethz.ch
% (C) 2002 Kari Unneland, Automatic Control Laboratory, ETH Zurich,
% (C) 2002 Mato Baotic, Automatic Control Laboratory, ETH Zurich,
%     kvasnica@control.ee.ethz.ch

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

narginchk(1, 2);

nu=size(Matrices.H,1);
nx=size(Matrices.E,2);

% set defaults for Matrices.F, Matrices.Cf, Matrices.Y, Matrices.Cx, Matrices.Cc
if(~isfield(Matrices,'F') | isempty(Matrices.F))
    Matrices.F=zeros(nx,nu);
end
if(~isfield(Matrices,'Cf') | isempty(Matrices.Cf))
    Matrices.Cf=zeros(1,nu);
end
if(~isfield(Matrices,'Y') | isempty(Matrices.Y))
    Matrices.Y=zeros(nx,nx);
end
if(~isfield(Matrices,'Cx') | isempty(Matrices.Cx))
    Matrices.Cx=zeros(1,nx);
end
if(~isfield(Matrices,'Cc') | isempty(Matrices.Cc))
    Matrices.Cc=0;
end
if ~isfield(Matrices, 'bndA'),
    Matrices.bndA = [];
    Matrices.bndb = [];
end

% make sure all matrices are in full format. this is important to do, because we
% call mpt_solvelpi which doesn't convert sparse matrices to full for solvers
% which don't support them
Matrices.F  = full(Matrices.F);
Matrices.Y  = full(Matrices.Y);
Matrices.H  = full(Matrices.H);
Matrices.G  = full(Matrices.G);
Matrices.E  = full(Matrices.E);
Matrices.W  = full(Matrices.W);
Matrices.Cf  = full(Matrices.Cf);
Matrices.Cx  = full(Matrices.Cx);
Matrices.Cc  = full(Matrices.Cc);
Matrices.bndA = full(Matrices.bndA);
Matrices.bndb = full(Matrices.bndb);

if ~isfield(Matrices, 'PbndIncluded'),
    % only include bndA and bnbd if not added in mpt_constructMatrices
    Matrices.E=[Matrices.E; -Matrices.bndA];
    Matrices.W=[Matrices.W; Matrices.bndb];
    Matrices.G=[Matrices.G; zeros(size(Matrices.bndA,1),size(Matrices.G,2))];
end

if size(Matrices.G, 2)==0,
    % we don't have anything to optimize over, return one region and associated
    % cost
    Pn = polytope(-Matrices.E, Matrices.W);
    Phard = Pn;
    details.Ai{1} = Matrices.Y;
    details.Bi{1} = Matrices.Cx;
    details.Ci{1} = Matrices.Cc;
    Fi{1} = zeros(0, size(Matrices.E, 2));
    Gi{1} = zeros(0, 0);
    activeConstraints = [];
    return
end    

if isfield(Matrices, 'constraints_reduced'),
    constraints_reduced = Matrices.constraints_reduced;
else
    constraints_reduced = 0;
end

if ~constraints_reduced,
    % we must remove redundant constraints otherwise it can happen that when
    % identifying active constraints, we choose several of them which are identical,
    % which in turn leads to problems when constructing critical regions.
    P = polytope([Matrices.G Matrices.E], Matrices.W);
    if ~isfulldim(P),
        % problem is infeasible
        Pn=polytope;
        Phard=polytope;
        Fi=[]; Gi=[]; activeConstraints=[];
        details.Ai = {}; details.Bi = {}; details.Ci = {};
        disp('mpt_mpqp:  Infeasible optimization problem from the begining');
        return
    end
    [H, K] = double(P);
    Matrices.G = H(:, 1:size(Matrices.G, 2));
    Matrices.E = H(:, size(Matrices.G, 2)+1:end);
    Matrices.W = K;
end

e = eig(Matrices.H);
if min(e) == 0
    disp('WARNING: Lack of strict convexity may lead to troubles in mpQP solver.')
elseif min(e) < -1e-8
    disp('WARNING: Problem is not positive semidefinite! Your mpQP solution may be completely wrong.')
elseif min(e) < 1e-5
    disp('WARNING: QP is close to being (negative) semidefinite, may lead to troubles in mpQP solver.')
end

[G,S,W,H,F,Hinv,GHinv,bndA,bndb,nx,nu,Matrices]=sub3_extract(Matrices);

global mptOptions;

if ~isstruct(mptOptions),
    mpt_error;
end

if(nargin<2)
    Options=[];
end
if ~isfield(Options,'verbose'),
    Options.verbose=mptOptions.verbose;
end
if(~isfield(Options,'useSymmetry'))
    Options.useSymmetry=0;  
end
if(~isfield(Options,'lpsolver'))
    Options.lpsolver=mptOptions.lpsolver;          %0: NAG ; 1: LINPROG ; 2: CPLEX; 3:CDD
end
if(~isfield(Options,'qpsolver'))
    Options.qpsolver=mptOptions.qpsolver;          %0: NAG ; 1: LINPROG ; 2: CPLEX
end
if(~isfield(Options,'debug_level'))
    Options.debug_level=mptOptions.debug_level;
end
if(~isfield(Options,'step_size'))
    Options.step_size = mptOptions.step_size; % How far to go from a border when finding a new region        
end
if(~isfield(Options,'abs_tol'))
    Options.abs_tol = mptOptions.abs_tol;
end
if(~isfield(Options,'rel_tol'))
    Options.rel_tol = mptOptions.rel_tol;
end
if(~isfield(Options,'center'))
    center = zeros(nx,1);
else
    center = Options.center;
end
if ~isfield(Options, 'statusbar'),
    Options.statusbar = 0;
end
if ~isfield(Options, 'zero_tol'),
    Options.zero_tol = 1e-13;     % everything below this value is considered zero
end

zero_tol = Options.zero_tol;

useSymmetry=Options.useSymmetry;  
lpsolver=Options.lpsolver;       
DEBUG=Options.debug_level;
stepSize=Options.step_size;


if(DEBUG>2)
    DEBUG=2;
end
if(DEBUG==2)
    if Options.verbose>0,
        disp('mpt_mpqp:')
        disp('!!!!!!!!!!!!!!!!!!!!!!ATTENTION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        disp('!!!!  DEBUG LEVEL 2 WILL INCREASE RUNTIME SIGNIFICANTLY   !!!!')
        disp('!!!!!!!!!!!!!!!!!!!!!!ATTENTION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    end
end
EMPTY_ROW_TOL = Options.abs_tol;  % Tolerance for declaring that a row is empty
CONSTR_TOL = Options.rel_tol;     % Tolerance for declaring that constr. is redundant

%----------------------------------------------------------------------
if(DEBUG==2)
    TOLERANCE=0;
else
    TOLERANCE=stepSize*10;
end

%-------------------------------------------------------------
% The limits on the states are considered as hard constraints
%--------------------------------------------------------------
hardA=bndA; % hard constraints
hardb=bndb; % hard constraints


activeConstraints={};   % list of active constraints
PEmpty= polytope; % normalized polyhedron description
Pn = PEmpty;
Fi={};      % control law
Gi={};      % control law
constraintStorage=[];   %storage structure for active constraints
xB=[]; %structure for storing all the points xBeyond
xBRegion=[]; %structure for storing region assosiated to point

options=[];


bboxOpt.noPolyOutput = 1;
BBoxes = struct('bmin', [], 'bmax', []);

%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                               
%%      FIND THE STARTING POINT                   
%% (IF ORIGINAL PROBLEM IS FEASIBLE)        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%THE FIRST REGION WILL BE THE UNCONSTRAINED ONE        
ii=[]; %no active constraints
%--------------------------------------------------------
% Calculating control law and region definition
%--------------------------------------------------------
nR = 1; % Number of regions 
[Pn, Fi, Gi, kr, fulldim] = sub1_computelaw(ii,nu,nx,Matrices,Options);
%%keptrows{nR} = getkeptrows(Pn);
keptrows{nR} = kr{1};
isemptypoly = ~fulldim;


if(isemptypoly)
    %-----------------------------------------------------
    %origin is infeasible, solve LP to get feasible point
    %-----------------------------------------------------
    if(~isempty(bndA))
        crA=[G -S; zeros(size(bndA,1),nu) bndA];
        crb=[W; bndb];
    else
        crA=[G -S];
        crb=W;
    end
    Pcr = polytope(crA, crb);
    [aux,R]=chebyball(Pcr,Options);
    if R<=0,
        Pn=polytope;
        Phard=polytope;
        Fi=[]; Gi=[]; activeConstraints=[];
        details.Ai = {}; details.Bi = {}; details.Ci = {};
        disp('mpt_mpqp:  Infeasible optimization problem from the begining');
        return
    end
    xFeasible=aux(nu+1:nu+nx);
    %================================================%
    % Solve optimization problem for point xFeasible %
    %================================================%
    options=[];

    % removing redundant constraints helps to improve robustness for
    % certain QP solvers (e.g. for quadprog)
    [A, b] = double(polytope(G, W+S*xFeasible));
    [zopt,lambda,how,exitflag,Vz] = mpt_solveQP(H, zeros(nu,1), ...
        A, b, [], [], [], Options.qpsolver, options);
    ii=find(abs(G*zopt-W-S*xFeasible)<=zero_tol);   % active constraints
    %ii=find(lambda>0);
    if exitflag<=0,
        error('mpt_mpqp: No controllable/feasible state found!!!');
    end
    %--------------------------------------------------------
    % Calculating control law and region definition
    %--------------------------------------------------------
    [Pn,Fi,Gi,kr,fulldim]=sub1_computelaw(ii,nu,nx,Matrices,Options);
    if ~fulldim,
        Pn=polytope;
        Phard=polytope;
        Fi=[]; Gi=[]; activeConstraints=[];
        details.Ai = {}; details.Bi = {}; details.Ci = {};
        disp('mpt_mpqp:  Infeasible optimization problem from the begining');
        return
    end
    for reg_ctr=1:length(Fi)
        keptrows{nR} = kr{reg_ctr};
        isemptypoly = ~isfulldim(Pn);
        if(isemptypoly)
            error('mpt_mpqp: Initial region is empty despite feasible constraints !!!');
        end
        if(length(Fi)>reg_ctr)
            nR=nR+1;
        end
    end
end

[d, bmin, bmax] = pelemfun(@bounding_box, Pn, bboxOpt);
BBoxes.bmin = [BBoxes.bmin [bmin{:}]];
BBoxes.bmax = [BBoxes.bmax [bmax{:}]];

activeConstraints{nR}=ii;
nii=length(ii);
constraintStorage=mpt_addToStorage(ii,constraintStorage);

%%%%%%%%%%%%%%%%%%%%%%%%
%%                    %%
%% DO THE EXPLORATION %%
%%                    %%
%%%%%%%%%%%%%%%%%%%%%%%%
region=1;

Tree = {};
treeCtr = 0;

if Options.statusbar,
    statbar.handle = mpt_statusbar('Computing...');
    statbar.progress = 0;
    Options.verbose = -1;
end

while region<=nR,
    %----------------------------
    % Writing out current region 
    %----------------------------
    if Options.statusbar & mod(region, 10)==0,
        statbar.progress = statbar.progress + region / 10000;
        if isempty(mpt_statusbar(statbar.handle, statbar.progress, 0, 0.9)),
            mpt_statusbar;
            error('Break...');
        end
    end
    if (mod(region,20)==0) & Options.verbose >= 1,
		fprintf('regions: %4i, unexplored: %i \n', ...
			region, nR-region);
    end 
    %-----------------------------------------------------
    % For each border in the region we are, check if 
    % there exist a new region outside
    %-----------------------------------------------------
    start=length(find(keptrows{region}<length(bndb)+1)); %ignore facets on outter bound
    for border=(start+1):nconstr(Pn(region))
        nRstart = nR;
        %-----------------------------------------------
        % Create a point close to the border
        %-----------------------------------------------
        [xBorder,RBorder]=facetcircle(Pn(region),border,Options);
        if(RBorder<0),
            if Options.verbose>1,
                disp('mpt_mpqp: Polyhedron does not exist => EMPTY REGION WAS STORED !!')
            end
            continue
            % Checkpoint !Should be taken into account!
        end %(RBorder<0)
        [Hp,Kp] = double(Pn(region));
        xBeyond=xBorder+stepSize*Hp(border,:)';
        
        if(~isempty(bndA) & ~all(bndA*xBeyond-bndb<=0))
            % Nothing happens, outside bounds
        else
            %--------------------------------------------------------------------
            % Check wheter point already exist inside a polyhedra
            %---------------------------------------------------------------------
            [nMatch,RegionStore]=sub6_RedundantPolyhedron(nR, Pn, ...
                region, xBeyond, BBoxes);
            %-----------------------------------------------------------------
            % Polyhedras should not overlap
            %-----------------------------------------------------------------
            if (nMatch>=1)
                if nMatch>1
                    %--------------------------------------------------------------
                    % When to overlapping regions are found:
                    % Check how big the overlap is and decide if it is significant 
                    %--------------------------------------------------------------
                    sub5_CheckOverlap(RegionStore,Pn,CONSTR_TOL,Options);
                end %nMatch>1,
            else
                %===========================================================================%
                % Solve optimization problem for point xBeyond and find active constraints
                %============================================================================%
                [zopt,lambda,how,exitflag,Vz]=mpt_solveQP(H,zeros(nu,1),G,W+S*xBeyond,[],[],[],Options.qpsolver,options);
                if exitflag==1
                    ii=find(abs(G*zopt-W-S*xBeyond)<=zero_tol);   % active constraints
                end
                %ii=find(lambda>0);
                %------------------------------------------------------------------
                % If the active constraints are equal to the region we came from, 
                % NOT POSSIBLE: increase stepsize and try again
                %--------------------------------------------------------------------
                kk=1;
                outsidebounds=0;
                
                while(mpt_iscombequal(ii,activeConstraints{region}) & exitflag>0 & ~outsidebounds)
                    xBeyond=xBorder+stepSize*2^kk*Hp(border,:)';
                    if(kk<5 & DEBUG==2)
                        xB(size(xB,1)+1,:)=xBeyond';
                        xBRegion(size(xB,1),:)=region;
                    end
                    %----------------------------------------------------------------
                    % Check if one is outside bounds: if end while
                    %---------------------------------------------------------------
                    if(~isempty(bndA) & ~all(bndA*xBeyond-bndb<=0))
                        outsidebounds=1;
                    else
                        [zopt,lambda,how,exitflag,Vz]=mpt_solveQP(H,zeros(nu,1),G,W+S*xBeyond,[],[],[],Options.qpsolver,options);
                        ii=find(abs(G*zopt-W-S*xBeyond)<=zero_tol);
                        %ii=find(lambda>0);
                        kk=kk+1;
                    end
                    if(kk>10000)
                        error(['mpt_mpqp: Maximum stepsize reached / Region: ',num2str(region)])
                    end 
                end% (ifa_iscombequal(ii,activeConstraints{region}) & kk<=30 & exitflag>0 &...
                % (nHard==0|all(hardA*xBeyond <= hardb))&all(bndA*xBeyond-bndb<=0))
                %-----------------------------------------------------------------------------
                % Check wheter active constraints already are covered
                %-----------------------------------------------------------------------------
                if(outsidebounds==0 & exitflag>0)
                    [constraintStorage,noNewFacetFound]=mpt_checkStorage(constraintStorage,ii);
                end
                %-----------------------------------------------------------------------------
                % If you have not found a new region, the exitflag and outsidebounds decide
                % whether you should add it as an outer bound
                %----------------------------------------------------------------------------        
                if (exitflag<=0 | outsidebounds==1 | noNewFacetFound),
                    if(exitflag<=0 & outsidebounds==0)
                        hardA(length(hardb)+1,:)=Hp(border,:);
                        hardb(length(hardb)+1,:)=Kp(border,:);
                    end
                    %----------------------------------------------------
                    % Add new region
                    %----------------------------------------------------
				elseif isequal(how, 'ok')
                    % Add control law
                    [Pcr,Fii,Gii,kr,fulldim]=sub1_computelaw(ii,nu,nx,Matrices,Options);
                    if ~fulldim
                        if Options.verbose > 1,
                            disp('WARNING: Empty controller region despite active constraints !')
                        end
                        break
					end
                    for reg_ctr=1:length(Fii)
                        nR = nR + 1;
                        activeConstraints{nR}=ii; %store active constraints
                        Fi{nR}=Fii{reg_ctr};
                        Gi{nR}=Gii{reg_ctr};
						Pn = [Pn Pcr(reg_ctr)];
                        [d, bmin, bmax] = bounding_box(Pn(nR), bboxOpt);
                        BBoxes.bmin = [BBoxes.bmin bmin];
                        BBoxes.bmax = [BBoxes.bmax bmax];

                        keptrows{nR} = kr{reg_ctr};
                    end
                    
                    %---------------------------------
                    % Remove redundant constraints
                    %-------------------------------
                    kk=1;
                    outsidebounds=0;
                    isemptypoly = ~isfulldim(Pcr);
                    while(isemptypoly & kk<50 & ~outsidebounds) %this is theoretically not possible !!
                        if(kk<5 & DEBUG==2)
                            xB(size(xB,1)+1,:)=xBeyond';
                            xBRegion(size(xB,1),:)=nR;
                        end
                        % Increase step size 
                        [Hp,Kp] = double(Pn(region));
                        xBeyond=xBeyond+stepSize*2^kk*Hp(border,:)';
                        % Check if one are outside bounds
                        if (~isempty(bndA) & ~all(bndA*xBeyond-bndb<=0))
                            outsidebounds=1;
                        else
                            % Check wheter point already exist inside a polyhedra
                            %%% nMatch=0;
                            [nMatch,RegionStore]=sub6_RedundantPolyhedron((nR-1), ...
                                Pn,region,xBeyond,BBoxes);
                            if (nMatch>=1)
                                outsidebounds=1;
                                if nMatch>1
                                    %--------------------------------------------------------------
                                    % When to overlapping regions are found:
                                    % Check how big the overlap is and decide if it is significant 
                                    %--------------------------------------------------------------
                                    sub5_CheckOverlap(RegionStore,Pn,CONSTR_TOL,lpsolver,DEBUG); 
                                end %nMatch>1,
                            else
                                %===========================================================================%
                                % Solve optimization problem for point xBeyond and find active constraints
                                %============================================================================%
                                [zopt,lambda,how,exitflag,Vz]=mpt_solveQP(H,zeros(nu,1),G,W+S*xBeyond,[],[],[],Options.qpsolver,options);
                                
                                if(exitflag<=0)
                                    hardA(length(hardb)+1,:)=Hp(border,:);
                                    hardb(length(hardb)+1,:)=Kp(border,:);
                                    outsidebounds=1;
                                elseif(mpt_iscombequal(find(lambda>0),ii))
                                    %Nothing happens, increase stepsize and try again
                                else
                                    ii=find(abs(G*zopt-W-S*xBeyond)<=zero_tol); 
                                    %ii=find(lambda>0);
                                    %Check if new region is found
                                    [constraintStorage,noNewFacetFound]=mpt_checkStorage(constraintStorage,ii);
                                    if(noNewFacetFound)
                                        outsidebounds=1;
                                    else
                                        activeConstraints{nR}=ii; % Store active constraints
                                        % Adding the control law
                                        [Pcr,Fi{nR},Gi{nR},kr,fulldim]=sub1_computelaw(ii,nu,nx,Matrices,Options);
										if fulldim
											for reg_ctr=1:length(Fii)
												Fi{nR}=Fii{reg_ctr};
												Gi{nR}=Gii{reg_ctr};
												Pn(nR) = Pcr(reg_ctr);
												[d, bmin, bmax] = bounding_box(Pn(nR), bboxOpt);
												BBoxes.bmin = [BBoxes.bmin bmin];
												BBoxes.bmax = [BBoxes.bmax bmax];
												
												keptrows{nR} = kr{reg_ctr};
												if(length(Fii)>reg_ctr)
													nR=nR+1;
												end
											end
										end
                                        
                                        isemptypoly = ~isfulldim(Pcr);
                                        if(isemptypoly==0)
                                            nMatch=0;
                                            [nMatch,RegionStore]= ...
                                                sub6_RedundantPolyhedron((nR-1), ...
                                                Pn,nR,xBeyond,BBoxes);
                                            if(nMatch>0)
                                                isemptypoly=1;
                                            end
                                        end
                                    end
                                end
                                
                            end
                        end
                        kk=kk+1;
                        if(kk==50)
                            error(['mpt_mpqp Maximum stepsize reached / Region: ',num2str(region)])
                        end
                    end 
                    if(~isemptypoly)
                        constraintStorage=mpt_addToStorage(activeConstraints{nR},constraintStorage);
                        %---------------------------------------------------------------
                        % If the constraints are symmetric the mirrored region is added
                        %---------------------------------------------------------------
                        if(useSymmetry)
                            nR=nR+1;
                            [activeConstraints{nR}]=sub7_SymmetricRegion(ii);
                            keptrows{nR}=keptrows{nR-1};
                            Pn(nR) = -Pn(nR-1);
                            [d, bmin, bmax] = bounding_box(Pn(nR), bboxOpt);
                            BBoxes.bmin = [BBoxes.bmin bmin];
                            BBoxes.bmax = [BBoxes.bmax bmax];
                            
                            Fi{nR} = Fi{nR-1};
                            Gi{nR} = -Gi{nR-1};
                            constraintStorage=mpt_addToStorage(activeConstraints{nR},constraintStorage);
                        end
                    else
                        %---------------------------------------------------------------
                        % Remove empty polyhedron
                        %--------------------------------------------------------------
                        Pn(nR) = [];
                        activeConstraints{nR}=[];
                        Fi{nR}=[];
                        Gi{nR}=[];
                        BBoxes.bmin(:, nR) = [];
                        BBoxes.bmax(:, nR) = [];
                        nR=nR-1;
                    end
                end
            end % nMatch==0
        end
        if(DEBUG==2)
            xB(size(xB,1)+1,:)=xBeyond'; % For each border, store explored point
            xBRegion(size(xB,1),:)=nR;
        end
    end % border
    region=region+1;
end % END EXPLORATION
%-----------------------------------------------------------------
% Remove outer constraints which do not belong to the outer hull
% ifa_adjustCellSize removes "empty" entries from the cell arrays
%-----------------------------------------------------------------

% MK: with the change made on Aug 31, 2012 it should no longer be necessary
% to call mpt_adjustCellSize, since we make sure that empty polytopes are
% not added
%[Pn,Fi,Gi,activeConstraints]=mpt_adjustCellSize(Pn,Fi,Gi,activeConstraints);

%!!!!!OBS OBS OBS!!!!!!
%The two next functions has to occur in the following way
%-----------------------------------------------------------------
% -Remove false outer bounds
% -Check if the outer hull contains empty regions
%-----------------------------------------------------------------

if Options.statusbar,
    mpt_statusbar(statbar.handle, 1);
end

if(DEBUG==1|DEBUG==2)
    [hardA,hardb,NoSolutionToPoint]=sub4_fixouterhull(hardA,hardb,xB,Pn,lpsolver,TOLERANCE,stepSize,xBRegion,Options);
end
if isempty(hardA),
    Phard = union(Pn);
    if length(Phard) > 1 | ~isfulldim(Phard)
        % union is non-convex or empty
        if Options.verbose >= 0,
            fprintf('mpt_mpqp: Warning: Phard is returned as an empty polytope. Please report this case to mpt@control.ee.ethz.ch\n');
        end
        Phard = polytope;
    end
else
    Phard = polytope(hardA, hardb);
end
%------------------------------------------------
% Remove redundant constraints in the outer hull
%-----------------------------------------------

%-------------------------------------------------------------------
% Check if one can find points with no solution in the outer hull
%------------------------------------------------------------------
if Options.verbose>0,
    disp(sprintf('mpt_mpqp: %d regions', nR));
end
Phard = set(Phard, 'minrep', 0);
Phard = reduce(Phard);
%%keptrows{nR} = getkeptrows(Phard);

for ii=1:length(Gi),
    Gi{ii} = Gi{ii} - Hinv*Matrices.Cf';
end

% Compute value function associated to each region
% (if output 'details' is expected)
%-------------------------------------------------
if nargout>=6,
    Ai = {};
    Bi = {};
    Ci = {};
    %The resulting value function is V(x)=x'*Ai{i}*x+Bi{i}*x+Ci{i} if x \in P(i)
    Matrices.H=Matrices.H/2; %to extract proper cost without 1/2
    for i=1:length(Fi),
        Ai{i}= Fi{i}'*Matrices.H*Fi{i}+ F*Fi{i} + Matrices.Y;
        Bi{i}= (2*Gi{i}'*Matrices.H+Matrices.Cf)*Fi{i}+Gi{i}'*Matrices.F'+Matrices.Cx;
        Ci{i}= Gi{i}'*Matrices.H*Gi{i}+Matrices.Cf*Gi{i}+Matrices.Cc;
    end
    details.Ai = Ai;
    details.Bi = Bi;
    details.Ci = Ci;
end    

if Options.statusbar,
    mpt_statusbar;
end

return


%------------------------------------------------------------------------------------------
% SUBFUNCTION 1
%------------------------------------------------------------------------------------------
function [Pcr,Fcontrol,Gcontrol,Krows,fulldim] = sub1_computelaw(ii,nu,nx,Matrices,Options)
%----------------------------------------------------------------
% Function which compute control law and region definition
%---------------------------------------------------------------
nii=length(ii);
if nii==0,
    %---------------------------------
    % Control law
    %----------------------------------
    Fcontrol{1}=zeros(nu,nx); 
    Gcontrol{1}=zeros(nu,1);  
    %----------------------------------- 
    % Region definiton
    %-----------------------------------
    crA=[Matrices.bndA; -Matrices.S];      
    crb=[Matrices.bndb;  Matrices.W];      
else
    Gt=Matrices.G(ii,:);
    if nii>nu,
        if Options.verbose>-1,
            disp('mpt_mpqp: Degeneracy')
        end
        %---------------------------------
        %Deal with primal degeneracy
        %---------------------------------
        [Gt,keptrows]=mpt_getFullRankSubset(Gt,1);
        if(iscell(Gt))
            Pcr=polytope;
            Fc = {};
            Gc = {};
            Kc = {};
            ctr=0;
            for k=1:length(keptrows)
                [Ptemp,Ft,Gt,krows] = sub1_computelaw(ii(keptrows{k}),nu,nx,Matrices,Options);
                if(isfulldim(Ptemp))
                    if(Pcr>=Ptemp)
                        %region is already covered
                    else
                        ctr=ctr+1;
                        Pcr=[Pcr Ptemp];
                        Fc{ctr}=Ft{1};
                        Gc{ctr}=Gt{1};
                        if iscell(krows),
                            Kc{ctr} = krows{1};
                        else
                            Kc{ctr} = krows;
                        end
                    end
                end
            end
            [Pret,how] = union(Pcr); %merge pieces
            if(how)
                %union is convex
                Pcr = Pret;
                Krows{1} = length(Matrices.bndb)+(1:nconstr(Pcr));
                Fcontrol{1}=Fc{1};
                Gcontrol{1}=Gc{1};
            else
                Fcontrol=Fc;
                Gcontrol=Gc;
                Krows=Kc;
			end
			fulldim=isfulldim(Pret);
            return
        else
            ii=ii(keptrows);
        end
    end

    % we can have serious problems if Gt contains any multiple identical rows,
    % therefore we remove them:
    [aa,bb] = unique(Gt, 'rows');
    % also restore order of constraints because unique changes it!
    Gt = Gt(bb, :);  
    ii = ii(bb);
        
    Wt=Matrices.W(ii,:);
    St=Matrices.S(ii,:);
    GHGinv=inv(Gt*Matrices.Hinv*Gt');
    tmat=Matrices.Hinv*Gt'*GHGinv;
    %-----------------------------------
    % Control law
    %-----------------------------------
    Fcontrol{1}=tmat*St;      
    Gcontrol{1}=tmat*Wt;      
    %-----------------------------------
    % Region definition
    %----------------------------------
    crA=[Matrices.bndA;  Matrices.G*Fcontrol{1}-Matrices.S;  GHGinv*St];   
    crb=[Matrices.bndb; -Matrices.G*Gcontrol{1}+Matrices.W; -GHGinv*Wt];   
    crA(length(Matrices.bndb)+ii,:)=[];       %remove active constraints from Gz<=W+Sx, since the inequality holds by definition
    crb(length(Matrices.bndb)+ii,:)=[];       %remove active constraints from Gz<=W+Sx, since the inequality holds by definition
end

% Switch to the uopt description: U=Z-Hinv*F'*x; %
Fcontrol{1}=Fcontrol{1}-Matrices.Hinv*Matrices.F';
Pcr = polytope(crA, crb, 0, 2);   % do not reduce the polytope yet
if ~isfulldim(Pcr),
	fulldim = false;
    kprows=[];
else
    [Pcr, kprows] = reduce(Pcr);       % perform reduction in order to access kept rows
	fulldim = true;
end
Krows{1} = kprows;
return

%-----------------------------------------------
% SUBFUNCTION 3
%-----------------------------------------------
function [G,S,W,H,F,Hinv,GHinv,bndA,bndb,nx,nu,Matrices]=sub3_extract(Matrices)
%------------------------------------------------------------------------------
% Function which extract the variables out of the structure Matrices
%------------------------------------------------------------------------------
G=Matrices.G; 
H=Matrices.H;
W=Matrices.W;
F=Matrices.F;
% if(isfield(Matrices,'S'))
%     S=Matrices.S;
% else
%     S=Matrices.E + G*inv(H)*F';
% end

if(isfield(Matrices,'bndA'))
    bndA=Matrices.bndA;
    bndb=Matrices.bndb;
else
    bndA=[];
    bndb=[];
    Matrices.bndA=[];
    Matrices.bndb=[];
end

H = (H + H')/2;
Hinv=inv(H);
Hinv = (Hinv + Hinv')/2;
GHinv=G*Hinv;
S=Matrices.E + GHinv*F';
if isfield(Matrices,'Cf'),
    W=W+GHinv*Matrices.Cf';
end
nx=size(S,2);
nu=size(H,1);

Matrices.S=S;
Matrices.W=W;
Matrices.Hinv=Hinv;
Matrices.GHinv=GHinv;
Matrices.nx=nx;
Matrices.nu=nu;

%--------------------------------------------------
% SUBFUNCTION 4
%--------------------------------------------------
function [hardA,hardb,NoSolutionToPoint2]=sub4_fixouterhull(hardA,hardb,xB,Pn,lpsolver,TOLERANCE,stepSize,xBRegion,Options);
%------------------------------------------------------------------------------------
% Function which:
% -checks if all the regions lies inside the computed outer hull,
%  if not, the bound which one are outside will be removed
%- "TOLERANCE" is decided from the debug features 
% (DEBUG = 1, tolerance is equal to stepsize
%  DEBUG = 2, tolerance is strictly zero     )
%------------------------------------------------------------------------------------
if isempty(hardA),
    % will be handled outside
    NoSolutionToPoint2 = [];
    return
end

len_Pn = length(Pn);
x = cell(1,len_Pn);
R = cell(1,len_Pn);
for region=1:len_Pn
    [x{region},R{region}]=chebyball(Pn(region),Options);
    remove=find(hardA*x{region}-hardb>0);
    if(~isempty(remove))
        %-----------------------------------------------------------
        % Remove the false outer bounds from the boundary structure
        %------------------------------------------------------------
        hardA(remove,:)=[];
        hardb(remove,:)=[];
    end
end

%------------------------------------------------------------------
% Check if there exist points which are not assosiated with a region
% inside the outer hull
%------------------------------------------------------------------
RemoveRegion=find(xBRegion>length(Pn));
xBRegion(RemoveRegion)=[]; %remove points not associated to region anymore
xB(RemoveRegion,:)=[];
NoSolutionToPoint=[];
isinOpt.abs_tol = 2*TOLERANCE;
isinOpt.fastbreak = 1;
if(~isempty(xB))
    for i=1:size(xB,1)
        xBeyond=xB(i,:)';
        index=xBRegion(i);
        FoundInsideRegion=0;
        if isinside(Pn(index),xBeyond,isinOpt),
            FoundInsideRegion=1;
            if(TOLERANCE==0 & ~all(hardA*xBeyond-hardb<=TOLERANCE))
                remove=find(hardA*xBeyond-hardb>0);
                if(~isempty(remove))
                    %-----------------------------------------------------------
                    % Remove the false outer bounds from the boundary structure
                    %------------------------------------------------------------
                    hardA(remove,:)=[];
                    hardb(remove,:)=[];
                end
            end
        else
            for regions=1:length(Pn)
                if isinside(Pn(regions), xBeyond, isinOpt),
                    FoundInsideRegion=1;
                    if(TOLERANCE==0 & ~all(hardA*xBeyond-hardb<=TOLERANCE))
                        remove=find(hardA*xBeyond-hardb>0);
                        if(~isempty(remove))
                            %-----------------------------------------------------------
                            % Remove the false outer bounds from the boundary structure
                            %------------------------------------------------------------
                            hardA(remove,:)=[];
                            hardb(remove,:)=[];
                        end
                    end
                end
            end
        end
        if(FoundInsideRegion==0)
            if(all(hardA*xBeyond-hardb<=TOLERANCE))
                NoSolutionToPoint(size(NoSolutionToPoint,1)+1,:)=xBeyond';
            end
        end
    end
end
%---------------------------------------------------------------
% If a point without solution is found, check if the gap is big,
% otherwise it is ok
%-------------------------------------------------------------
NoSolutionToPoint2=[];
isinOpt.abs_tol = 0;
isinOpt.fastbreak = 1;

if(TOLERANCE==0 & ~isempty(NoSolutionToPoint))
    nx=size(xBeyond,1);
    for i=1:size(NoSolutionToPoint,1)
        nMatch=0;
        for k=1:nx
            Border=zeros(1,nx);
            Border(k)=1;
            xBeyond1=NoSolutionToPoint(i,:)'+stepSize*100*Border';
            xBeyond2=NoSolutionToPoint(i,:)'-stepSize*100*Border';
            %---------------------------------------------------
            % Check if point is inside another polyhedra
            %---------------------------------------------------
            [dummy1, inwhich] = isinside(Pn, xBeyond1, isinOpt);
            nMatch1 = length(inwhich);
            [dummy1, inwhich] = isinside(Pn, xBeyond2, isinOpt);
            nMatch2 = length(inwhich);
            if(nMatch1==0 | nMatch2==0)
                nMatch=0;
            end
        end
        if(nMatch==0)
            NoSolutionToPoint2=[NoSolutionToPoint(i,:);NoSolutionToPoint2];
        end
    end
end


if(~isempty(NoSolutionToPoint2))        
    if Options.verbose>0
        disp('****************************************************************************************')
        disp('mpt_mpqp: Points where QP have no feasible solution inside the hull:')
        disp(num2str(NoSolutionToPoint2));
        disp('****************************************************************************************')
    end
end

%-------------------------------------------------------------------
% SUBFUNCTION 5
%------------------------------------------------------------------
function [overlap]=sub5_CheckOverlap(RegionStore,Pn,CONSTR_TOL,Options)
%------------------------------------------------------------------
% When overlapping regions are found:
% Checks how big the overlap of two regions are, and gives a message
% if it is significantd
%--------------------------------------------------------------------
overlap=[];
overlaplimit=1e-4; 
Options.reduce_intersection=0;
for t=1:(length(RegionStore)-1)
    for s=t+1:length(RegionStore)
        Pt = Pn(RegionStore(t));
        Ps = Pn(RegionStore(s));
        %compatibility%crP = intersect(Pt,Ps,Options);
		crP = intersect(Pt, Ps);
        [x,R]=chebyball(crP,Options);
        overlap(t)=R;
        if(R>overlaplimit)
            if Options.verbose>-1,
                disp('mpt_mpqp: Polyhedra should not overlap');
                disp(['          Region ',num2str(RegionStore(t)),' and ',num2str(RegionStore(s)),'  Chebychew radius: ',num2str(R)]);
            end
        else
            %verbose nothing
        end
    end
end

%------------------------------------------------------------------
% SUBFUNCTION 6
%-------------------------------------------------------------------
% Checks if a point lies inside one or more polyhedrons
%-------------------------------------------------------------------
function [nMatch,RegionStore]=sub6_RedundantPolyhedron(nR,Pn,region,xBeyond,BBoxes)
%-------------------------------------------------------------------
% Checks if a point lies inside one or more polyhedrons
%-------------------------------------------------------------------
nMatch=0;
RegionStore=[];
isinOpt.abs_tol = 0;
isinOpt.fastbreak = 1;

% use bounding boxes for pruning. "ind" will be a logical array of 1/0
% indicies of possible candidate regions
%
% we introduce a tolerance of 1e-4 since bounding boxes are not always
% precise
LOWER = BBoxes.bmin;
UPPER = BBoxes.bmax;
%repX = repmat(xBeyond, 1, size(LOWER,2));
repX = xBeyond(:, ones(1, size(LOWER, 2))); % same as repmat, but faster
myfind = find(all((repX >= LOWER-1e-4) & (repX <= UPPER+1e-4), 1));
[ii,iwhere] = isinside(Pn(myfind), xBeyond, isinOpt);
iwhere = myfind(iwhere);

for ii=1:length(iwhere),
    if iwhere(ii)~=region & iwhere(ii)<=nR,
        nMatch=nMatch+1;
        RegionStore(length(RegionStore)+1)=iwhere(ii);
    end
end

%---------------------------------------------------------------------
% SUBFUNCTION 7
%-----------------------------------------------------------------------
% Calculating constraints which should be added when the system is symmetric
%-----------------------------------------------------------------------
function [ii2]=sub7_SymmetricRegion(ii)
ii2=[];
for j=1:length(ii)
    if(mod(ii(j),2)==1)
        ii2(j) = ii(j)+1;
    else
        ii2(j) = ii(j)-1;
    end
end
ii2=ii2';

%---------------------------------------------------------------------
% SUBFUNCTION 8
%-----------------------------------------------------------------------
function [constraintStorage,noNewFacetFound]=mpt_checkStorage(constraintStorage,activeConstraint)   
% [constraintStorage,noNewFacetFound]=mpt_checkStorage(constraintStorage,activeConstraint)   
%
%-----------------------------------
%   DESCRIPTION:
%-----------------------------------
% Check if region (ie active constraints) already exists
% Internal routine
%
%-----------------------------------
%    INPUT:
%-----------------------------------   
%  constraintStorage:           storage structure
%  activeConstraint:            constraint to be added to structure
%
%-----------------------------------
%   OUTPUT:
%-----------------------------------
%  constraintStorage:           storage structure
%  noNewFacetFound:             combination/activeConstraints already stored
%
%(C) 2003 Pascal Grieder, Automatic Control Laboratory, ETH Zurich, grieder@control.ee.ethz.ch
%
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


if(isempty(activeConstraint))
    activeConstraint = 0; 
end

search=1;
noNewFacetFound = 0;

try 
    if(isempty(constraintStorage{length(activeConstraint)}{max(min(activeConstraint),1)}))
        search=0;
    end
catch
    search=0;
end

if(search)
    storageSegment =  constraintStorage{length(activeConstraint)}{max(min(activeConstraint),1)};
    stLength       =  length(storageSegment);
    i=1;
    while(i<=stLength)
        if(length(activeConstraint)==length(storageSegment{i}))
            if(ismember(activeConstraint, storageSegment{i}))
                %combination already stored
                noNewFacetFound = 1;
                return
            end
        end
        i=i+1;
    end
end%if search



%---------------------------------------------------------------------
% SUBFUNCTION 9
%-----------------------------------------------------------------------
function [Pn,Fi,Gi,activeConstraints] = mpt_adjustCellSize(Pn,Fi,Gi,activeConstraints)
% [Pn,Fi,Gi,activeConstraints] = mpt_adjustCellSize(Pn,Fi,Gi,activeConstraints)
%
%-----------------------------------
%   DESCRIPTION:
%-----------------------------------
% this function removes "empty" entries from the cell arrays
% Internal routine.
%
%-----------------------------------
%    INPUT:
%-----------------------------------   
% Pn                 - polytope array defining polyhedral partition
% Fi,Gi              - cell arrays defining the PWA control law
% activeConstraints  - cell array containing active constraints for each region
%
%-----------------------------------
%   OUTPUT:
%-----------------------------------
% Pn                 - polytope array defining polyhedral partition
% Fi,Gi              - cell arrays defining the PWA control law
% activeConstraints  - cell array containing active constraints for each region
%
%(C) 2003 Pascal Grieder, Automatic Control Laboratory, ETH Zurich, grieder@control.ee.ethz.ch
%
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

ctr=1;
PPn = polytope;
FFn = {};
GGn = {};
constr = {};
for i=1:length(Pn)
    if isfulldim(Pn(i)),
        PPn = [PPn Pn(i)];
        GGn{ctr}=Gi{i};
        FFn{ctr}=Fi{i};
        constr{ctr}=activeConstraints{i};
        ctr=ctr+1;
    end
end
clear Pn Fi Gi activeConstraints
Pn=PPn;
Fi=FFn;
Gi=GGn;
activeConstraints=constr;



%---------------------------------------------------------------------
% SUBFUNCTION 10
%-----------------------------------------------------------------------
function [constraintStorage]=mpt_addToStorage(activeConstraint,constraintStorage)   
% [constraintStorage]=mpt_addToStorage(activeConstraint,constraintStorage,horizon)
%
%-----------------------------------
%   DESCRIPTION:
%-----------------------------------
% adds the set of active constraints to the storage structure
% Internal routine
%
%-----------------------------------
%    INPUT:
%-----------------------------------   
% activeConstraints
% constraintStorage
%
%-----------------------------------
%   OUTPUT:
%-----------------------------------
% constraintStorage
%
%(C) 2003 Pascal Grieder, Automatic Control Laboratory, ETH Zurich, grieder@control.ee.ethz.ch
%
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
if(isempty(activeConstraint))
    activeConstraint = 0; 
end

try
    if(isempty(constraintStorage{length(activeConstraint)}{max(min(activeConstraint),1)}))
        segmentLength=0;
    else
        segmentLength =  length(constraintStorage{length(activeConstraint )}{max(min(activeConstraint),1)});
    end
catch
    segmentLength=0;
end


if(length(activeConstraint )>0)
    constraintStorage{length(activeConstraint)}{max(min(activeConstraint),1)}{segmentLength+1}= activeConstraint';
end
