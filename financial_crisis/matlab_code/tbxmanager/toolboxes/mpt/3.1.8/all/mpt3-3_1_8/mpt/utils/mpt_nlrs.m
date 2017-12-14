function [V,R,H,K] = mpt_nlrs(method,VP,Options)
%MPT_LRS Matlab implementation of the LRS algorithm
%
% [V,R] = mpt_nlrs(method, VP, Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% enumerates extreme points or extreme facets of a polytope, works only
% with full-dimensional polyhedra
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% method   - string variable
%            'extreme' - enumeration of extreme points
%            'hull'    - computation of convex hull
% VP       - input argument - a Polyhedron object
%
% ---------------------------------------------------------------------------
% OUTPUT
% ---------------------------------------------------------------------------
% V        - output argument (either set of extreme points or a Polyhedron)
% R        - rays (only for extreme point computation)
%
%
% see also POLYHEDRON/EXTREME, HULL

% ---------------------------------------------------------------------------
% LITERATURE
% ---------------------------------------------------------------------------
% based on code published in:
%
% CONVEX HULL CALCULATIONS: A Matlab IMPLEMENTATION AND CORRECTNESS PROOFS
% FOR THE LRS-ALGORITHM
%
% http://www.mat.uc.pt/preprints/ps/p0326.pdf
%
% Authors: ALEXANDER KOVACEC AND BERNARDETE RIBEIRO
% Corresponding Author:
% Alexander Kovacec, Dep. Matematica, Univ. Coimbra, 3001-454 Coimbra,
% Portugal. kovacec@mat.uc.pt.
%
% Pre-Publicacoes do Departamento de Matematica, Universidade de Coimbra,
% Preprint Number 03/26
%

% Copyright is with the following author(s):
%
% 2012, Revised by M. Herceg, Automatic Control Laboratory, ETH Zurich,
%         herceg@control.ee.ethz.ch
%(C) 2004 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%         kvasnica@control.ee.ethz.ch
%(C) ALEXANDER KOVACEC AND BERNARDETE RIBEIRO

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


if nargin<3
    Options = [];
end
if ~isfield(Options,'lrsClearGlobals'),
    Options.lrsClearGlobals = 1;
end

if Options.lrsClearGlobals
    clear global lrsDictionary lrsBasis NlrsBasis
end

if ~ischar(method),
    error('mpt_nlrs: The first argument must be a string!');
end

if ~isa(VP,'Polyhedron')
    error('mpt_nlrs: The second argument must be a polyhedron object.');
end

if numel(VP)==0
    V = [];
    R = [];
    H = [];
    K = [];
    return
end

if numel(VP)>1
    error('mpt_nlrs: The polyhedron must be of length one.');    
end
if ~isFullDim(VP)
    error('mpt_nlrs: The polyhedron must be full-dimensional.');
end


if strcmpi(method,'extreme'),
    V = [];
    R = [];
    H = VP.A;
    K = VP.b;
    initv = initvert(H,K);
    if ~isempty(initv),
        [V,R] = htovr(H,K,initv);
        V = V';
        R = R';
    end
    
elseif strcmpi(method,'hull'),
    [H,K] = vrtoh(VP.V',VP.R');
    V = VP.V;
    R = VP.R;
    
else
    error(['mpt_nlrs: uknown method ''' method '''!']);
end

if Options.lrsClearGlobals
    clear global lrsDictionary lrsBasis NlrsBasis
end

end



function [V,R]=htovr(H,b,initv)
% [V,R]=htovr(H,b,initv)
% in: matrices H,b, and vertex-information as produced by
% initv=initvert(H,b).
% out: a collection of vertices and rays with information about
% their origin.
global lrsDictionary lrsBasis NlrsBasis
[lrsDictionary, lrsBasis, NlrsBasis]=tratodic(H,b,initv);
%[V,R]=lrs(lrsDictionary, lrsBasis, NlrsBasis);
[V,R] = lrs;

end

function initv=initvert(H,b)
% initv=initvert(H,b) in: an mxn-matrix H and an m-column b.
% out: if exists, a pair initv=[x0,I] of vertex x0 of the
% polyhedron Hx<=b; and a n-tuple I of positive integers
% so that H(I,:) has rank n and H(I,:)*x0=b; otherwise
% x0=[], I=[], and appropriate messages.
initv = [];
[m,n]=size(H);
x0=[];
I=[];
if rank(H)<n
    %disp('polyhedron Hx<=b has no vertex');
    return;
end
cnt=1; cl=10;
while 1
    if cnt<=cl
        p=randperm(m);
        p=p(1:n);
        z=zeros(1,m);
        z(p)=ones(1,n);
        cnt=cnt+1;
    end
    if cnt==cl+1
        z=[ones(1,n), zeros(1,m-n)];
        cnt=cl+2;
        %disp('Systematic vertex search begins. This may take time.');
    end
    T=H(find(z),:);
    if rank(T)==n
        x=T\(b(find(z)));
        u=b-H*x;
        f=(u>=-10*eps);
        if all(f)
            x0=x;
            I=z;
            %disp('vertex found');
            initv=[x0, (find(I))'];
            return;
        end
    end
    if cnt==cl+2
        if z(m-n+1:m)==1
            %disp('polyhedron Hx<=b is empty');
            return;
        end
        z=nextset(z);
    end
end
end

function boolean=lexmin(v)
% boolean=lexmin: in: v=0 or v in NlrsBasis. out: boolean=0 or boolean=1.
% if v=0: boolean=1 iff lrsBasis is lexmin for a basic feasible solution.
% if v in NlrsBasis: boolean=1 iff lrsBasis is lexmin for a geometric ray
% represented by v.
global lrsDictionary lrsBasis NlrsBasis
[sD1,sD2]=size(lrsDictionary);
b=lrsDictionary(:,sD2);
boolean=1;
if v==0
    for s=NlrsBasis
        for t=find(lrsBasis>s)
            if b(t)==0 & lrsDictionary(t,s)~=0
                boolean=0;
                return;
            end,
        end,
    end,
end
if ~(v==0)
    for s=NlrsBasis
        for t=find(lrsBasis>s)
            if b(t)==0 & lrsDictionary(t,s)~=0 & lrsDictionary(t,v)==0
                boolean=0;
                return;
            end,
        end,
    end,
end
end

function I=lmv(M)
% I=lmv(M) finds, given a real matrix M indices of rows that are
% lexicographically smallest.
j=1;
I=1:size(M,1);
n=size(M,2);
while (length(I)>=2)&(j<=n)
    col=M(I,j);
    m=min(col);
    sI=find(col==m);
    I=I(sI);
    j=j+1;
end
end

%function [vertices, rays]=lrs(lrsDictionary,lrsBasis,NlrsBasis)
function [vertices, rays]=lrs
% [vertices, rays]=lrs(lrsDictionary,lrsBasis,NlrsBasis): given a valid
% triple (D,B,N) as input, as obtained from tratodic.m, returns
% all vertices and extreme rays of associated polytope in
% 'vertices' and 'rays' in form origin ray origin ray ... '.
global lrsDictionary lrsBasis NlrsBasis
vertices = [];
rays = [];
[sD1,sD2]=size(lrsDictionary);
j=1;
lNb=sD2-sD1-1;
while 1
    while j<=lNb
        v=NlrsBasis(j);
        [brev, u]=reverse(v);
        if ~brev
            if u==0&lexmin(v)
                rays=[rays, lrsDictionary(2:lNb+1,sD2) -lrsDictionary(2:lNb+1, v)];
            end
            j=j+1;
        end
        if brev
            pivot(u,v);
            if lexmin(0),
                vertices=[vertices, lrsDictionary(2:lNb+1,sD2)];
            end
            j=1;
        end
    end
    [r,j]=slctpivo;
    if isempty(j),
        vertices=[vertices,lrsDictionary(2:lNb+1,sD2)];
        return;
    end
    pivot(r,NlrsBasis(j));
    j=j+1;
end
end

function r=lxminrat(s)
% r=lxminrat(s): given s in NlrsBasis and lexpositive lrsBasis calculates
% integer r=lexminratio(lrsBasis,s) in sense of [A, p185c6...11].
% in particular r==0 iff s represents a ray, an r in lrsBasis otherwise.
global lrsDictionary lrsBasis NlrsBasis
[sD1, sD2]=size(lrsDictionary);
lNb=sD2-sD1-1;
a=lrsDictionary(:,s);
I=find(a(lNb+2:sD1)>0)+(lNb+1);
if isempty(I)
    r=0;
    return;
end
D=lrsDictionary(:,[sD2, 1:sD1]);
Dtilde=diag(1./a(I))*D(I,:);
i=lmv(Dtilde);
t=I(i);
r=lrsBasis(t);

end

function nextx=nextset(x)
%input: a {0,1}-n-tuple x i.e. characteristic vector of
% subset of 1...n
%output: lexicographically next n-tuple with same number
% of ones as x; x itself if x is lexicographically
% last element.
lengthx=length(x);
last0=max(find(x==0));
xinit=x(1:last0-1);
xfin1s=x(last0+1:lengthx);
j=max(find(xinit==1));
if isempty(j)
    nextx=x;
    return;
end
nextx=x;
nextx([j,j+1])=[0 1];
nextx=[nextx(1:j+1),xfin1s];
nextx=[nextx, zeros(1,lengthx-length(nextx))];

end


function pivot(r,s)
% pivot(r,s) given r in lrsBasis and s in NlrsBasis, pivots lrsDictionary
% for lrsBasis to that for new lrsBasis lrsBasis-r+s according to
% [A, p183c-1]; updates lrsBasis and NlrsBasis
global lrsDictionary lrsBasis NlrsBasis
sD1=size(lrsDictionary,1);
a=lrsDictionary(:,s);
t=find(lrsBasis==r);
tn=find(NlrsBasis==s);
if ~isempty(t),
    rowt=lrsDictionary(t,:)/a(t);
    lrsDictionary(t,:)=rowt;
    for i=[1:t-1,t+1:sD1]
        lrsDictionary(i,:)=lrsDictionary(i,:)-a(i)*rowt;
    end
    lrsBasis(t)=s; NlrsBasis(tn)=r;
    lrsDictionary=lrsDictionary.*(abs(lrsDictionary)>100*eps);
    lrsDictionary(:,s)=zeros(sD1,1);
    lrsDictionary(t,s)=1;
end
end


function [brev,u]=reverse(v)
% u=reverse(v)
% in: v in NlrsBasis. out: u=lxminrat(v) and brev\in {0,1}.
% brev==1: if v represents an edge and the lexicographic pivot
% rule applied to B-u+v generates a pivot back to B.
% brev==0: in all other cases.
global lrsDictionary lrsBasis NlrsBasis
u=lxminrat(v);
if u==0
    brev=0;
    return,
end
i=find(lrsBasis==u);
a=lrsDictionary(:,v);
a0=a(1);
ai=a(i);
w0=lrsDictionary(1,:);
wi=lrsDictionary(i,:);
wbar=w0-(a0*wi/ai);
J=NlrsBasis(NlrsBasis<u);
if isempty(J)
    wbarg0=1;
else,
    wbarg0=all(wbar(J)>=-10*eps);
end
if (w0(v)>0 & wbarg0)
    brev=1;
else,
    brev=0;
end
end

function [xk,maxim]=simplex(c,H,b,initv)
% [xk,maxim]=simplex(c,H,b,initv)
% in: d-row c, mxd-matrix H, m-column b, dx2-column
% initv=initvert(H,b), so that P={x: Hx<=b} is a polytope
% having initv(:,1) as a vertex.
% out: vertex xk of P and real maxim=c*xk=max{cx: Hx<=b}
% or a message informing about the unboundedness of
% the problem.
global lrsDictionary lrsBasis NlrsBasis
[lrsDictionary, lrsBasis, NlrsBasis]=tratodic(H,b,initv);
[sD1,sD2]=size(lrsDictionary);
lNb=sD2-sD1-1;
u=c*lrsDictionary(2:lNb+1,sD1+1:sD2);
lrsDictionary(1,:)=[1 zeros(1,sD1-1), u(1:lNb), u(sD2-sD1)];
[r,j]=slctpivo;
while ~isempty(j)
    if r==0
        %disp('The problem is unbounded');
        return;
    end
    pivot(r,NlrsBasis(j));
    [r,j]=slctpivo;
end
xk=lrsDictionary(2:lNb+1,sD2);
maxim=lrsDictionary(1, sD2);
end

function [r,j]=slctpivo
% [r,j]=slctpivo: if measuring a non-optimal lrsDictionary, positive
% indices r in lrsBasis and j are returned so that r,s=NlrsBasis(j)
% reflect the lex pivot selection. A lrsDictionary is optimal
% iff j=[] is returned.
global lrsDictionary lrsBasis NlrsBasis
s=min(find(lrsDictionary(1,:)<0));
if ~isempty(s),
    j=find(NlrsBasis==s);
else
    j = [];
end
if isempty(j)
    r=0;
    return
end
r=lxminrat(s);
end

function [lrsDictionary, lrsBasis,NlrsBasis]=tratodic(H,b,initv)
%[lrsDictionary, lrsBasis,NlrsBasis]=tratodic(H,b,initv).
%in: real mxd-matrix H, m-column b, and dx2-matrix initv=[v,I],
% with I a column, of d distinct positive integers so that
% Hx<=b defines a nonempty polyhedron having v as a vertex,
% H(I,:) has rank d, and H(I,:)*v=b(I).
% If existing, initv can be produced by initvert.m
% out: a lrsDictionary digestible by lrs-function in lrs.m
[m,d]=size(H);
I=zeros(m,1);
I(initv(:,2))=ones(d,1);
H=[H(find(1-I),:); H(find(I),:)];
b=[b(find(1-I),:); b(find(I))];
Hupp=H(1:m-d,:);
bupp=b(1:m-d);
Hlow=H(m-d+1:m,:);
blow=b(m-d+1:m);
invHlow=inv(Hlow);
lrsDictionary=[eye(d) zeros(d,m-d) invHlow invHlow*blow; zeros(m-d,d) eye(m-d) -Hupp*invHlow bupp-Hupp*invHlow*blow];
lrsDictionary=[ 1 zeros(1,m) ones(1,d) 0; zeros(m,1) lrsDictionary ];
lrsDictionary=lrsDictionary.*(abs(lrsDictionary)>100*eps);
lrsBasis=1:(m+1); NlrsBasis=m+2:m+d+1;
end

function [H,b]=vrtoh(V,R)
% [H,b]=vrtoh(V,R) in: vertices as columns of V, ray(directions)
% as columns in R
% out: IF vertices define a full dimensional body then a
% H-description Hx<=b of the polyhedron. Otherwise
% an error maessage.
global lrsDictionary lrsBasis NlrsBasis
E=V-V(:,1)*ones(1,size(V,2));
E(:,1)=[];
E=[E,R];
if rank(E)<size(E,1)
    error('mpt_nlrs: VR-representation is not full dimensional');
end
V=-V';
m=size(V,1);
V=[-ones(m,1),V];
R=-R';
m=size(R,1);
R=[zeros(m,1),R];
A=[V;R];
b=zeros(size(A,1),1);
initv=initvert(A,b);
[lrsDictionary, lrsBasis, NlrsBasis]=tratodic(A,b,initv);
%[V,R]=lrs(lrsDictionary, lrsBasis, NlrsBasis);
[V,R] = lrs;
Rays=R(:,2:2:size(R,2));
Rays=Rays';
H=-Rays(:,2:size(Rays,2));
b=Rays(:,1);
end