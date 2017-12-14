%===============================================================================
%
% Title:        findFAS
%                                                             
% Project:      Transformation of HYSDEL model into PWA model
%
% Purpose:      find feedback arc set (FAS)
%
% Input:        A: adjacency matrix
%
% Output:       FAS: feedback arc set: FAS(i,:) contains the i-th feedback arc 
%                    from vertex FAS(i,1) to vertex FAS(i,2)
%               s_fas: vertex sequence corresponding to A_fas (the i-th entry of 
%                      s_fas contains the number of the i-th vertex and 
%                      corresponds to the i-th row and column of A, the 
%                      vertexTable maps the numbers of the vertices to their 
%                      variable names)
%
% Comments:     The algorithm is based on 
%               'A fast and effective heuristic for the feedback arc set problem', 
%               P. Eades, X. Lin, W.F. Smyth, 
%               Information Processing Letters 47 (1993) 319-323
%
% Authors:      Tobias Geyer <geyer@control.ee.ethz.ch>
                                                                 
% History:      date        subject                                       
%               2003.01.??  initial version 
%
% Contact:      Tobias Geyer
%               Automatic Control Laboratory
%               ETH Zentrum, CH-8092 Zurich, Switzerland
%
%               geyer@aut.ee.ethz.ch
%
%               Comments and bug reports are highly appreciated
%
%===============================================================================


function [FAS, s_fas] = findFAS(A);

% for debugging:
%A = [0 1 1 0; 0 0 0 1; 0 1 0 1; 1 0 0 0];

% vertex sequence
s = 1:length(A);

% store the original digraph
A_org = A;
s_org = s;


% find feedback arc set (FAS)
s1 = [];
s2 = [];
while length(A) > 0
    
    % remove all sinks
    sinks = findSinks(A, s);
    while ~isempty(sinks)
        u = sinks(1);                           % choose the first one
        s2 = [u s2];                            % add to s2
        [A, s] = remVertex(u, A, s);            % and remove from the digraph
        sinks = findSinks(A, s);
    end;
    
    % remove all sources
    sources = findSources(A, s);
    while ~isempty(sources)
        u = sources(1);                         % choose the first one
        s1 = [s1 u];                            % add to s1
        [A, s] = remVertex(u, A, s);            % and remove from the digraph
        sources = findSources(A, s);
    end;
    
    % remove vertex u for which delta(u) is a maximum
    u = findMaxVertex(A, s);
    s1 = [s1 u];                                % add to s1
    [A, s] = remVertices(u, A, s);              % and remove from the digraph
    
end;

s_fas = [s1 s2];


% now, we have the following:
% the original digraph A_org with the original vertex sequence s_org
% and the new vertex sequence from the FAS algorithm s_fas
% the new digraph A_fas is missing

% build A_fas by permuting A_org and s_org according to s_org and s_fas
A_fas = A_org;
for i=1:length(A_org)
    if s_org(i) ~= s_fas(i)
        % i.e. the entry of s_org at position i is wrong
        j = find( s_org == s_fas(i) );
        [A_fas, s_org]  = permuteA( A_fas, s_org, i, j );
    end;
end;

% check whether the permuted sequence and the one resulting from the FAS 
% algorithm are the same - if not: error
if ~all(s_fas == s_org), error('Error in sequence permutation.'); end;


% derive the feedback arc set (FAS) 
% and the ajacency matrix witout feedback arcs (A_wofa)
% remark: the FAS corresponds to ones in the lower triangle of A_fas
% structure of the FAS:
% every row holds a feedback arc, where the first entry is the number of 
% the source vertex and the second entry the number of the sink vertex
% (the following always holds: number of sink < number of source vertex)
FAS = []; 
A_wofa = A_fas;
for r=1:length(A_fas)
    for c=1:r
        if A_fas(r,c) == 1
            % we have found a feedback arc
            if r == c
                % we have found a self-loop: error
                error('found self-loop in adjacency matrix')
            end;
            FAS(end+1,:) = [s_fas(r), s_fas(c)];
            A_wofa(r, c) = 0;
        end;
    end;
end;




function v = findSinks(A, s)
% find the vertex numbers of the sinks in the digraph
% a sink has outdegree 0, i.e. the v-th row has only zeros
a_entries = find( sum(A,2) == 0);
v = s(a_entries);
return

function v = findSources(A, s)
% find the vertex numbers of the sources in the digraph
% a source has indegree 0, i.e. the v-th column has only zeros
a_entries = find( sum(A,1) == 0);
v = s(a_entries);
return

function v = findMaxVertex(A, s)
% find the vertex number for which delta=outdegree-indegree
% is a maximum
delta_entries = sum(A,2)' - sum(A,1);
[d_max, entry] = max(delta_entries);
v = s(entry);
return

function [A, s] = remVertex(v, A, s)
% remove the vertex with the number v from the digraph
% by deleting the corresponding row and column of A and by removing
% the corresponding entry in s
v_entry = find( s==v );     
A(v_entry,:) = [];
A(:,v_entry) = [];
s(v_entry) = [];
return

function [A, s] = permuteA(A, s, i, j)
% swap the i-th and j-th entries in s and
% modify A accordingly (swapping the i-th and j-th rows and columns yields
% the desired result)
buf = s(i);   s(i)   = s(j);   s(j)   = buf;
buf = A(:,i); A(:,i) = A(:,j); A(:,j) = buf;
buf = A(i,:); A(i,:) = A(j,:); A(j,:) = buf;
return




% old and not used anymore

function [A, A_num] = remVertices(v_num, A, A_num);
% remove the vertices with the numbers v_num from the digraph
% by deleting the corresponding rows and columns of A and by removing
% the corresponding entries in A_num
for i=1:length(v_num)
    v_entry = find( A_num==v_num(i) );
    A(v_entry,:) = [];
    A(:,v_entry) = [];
    A_num(v_entry) = [];
end;
return

% remove all sinks
u = findSinks(A, A_num);                    % find all sinks
s2 = [u s2];                                % add to s2
[A, A_num] = remVertices(u, A, A_num);      % and remove from the digraph

% remove all sources
u = findSources(A, A_num);                  % find all sources
s1 = [s1 u];                                % add to s1
[A, A_num] = remVertices(u, A, A_num);      % and remove from the digraph