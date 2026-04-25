function [c4n,n4e,nodes_line,n4sDb] ...
    = refineRGB_with_line(c4n,n4e,nodes_line,sides_line,n4sDb,ms)
%% refineRGB
% [c4n,n4e,n4sDb,n4sNb,Pe4e,Pn4n]
%            = refineRGB(c4n,n4e,n4sDb,n4sNb,ms)
%
% Refine marked sides of a grid ("red/green/blue").
% Marked sides are given with their node numbers in n4ms
% (a [#Marked-sides 2] vector with n4s-structure).
% A closure algorithm is applied to n4ms first to ensure that
% the reference sides of all marked elements are marked.
%
% ms can be a simple list ([? x 1]!!) of number or marked
% (numbering as in n4s), or a list of node pairs ([? x 2]),
% i.e., a (row-)subset of n4s.

nr_nodes_old = size(c4n,1);
n4s = computeN4s(n4e);
s4n = computeS4n(n4e,n4s);
s4e = computeS4e(n4e);

if size(ms,2)==2  % ms is really n4ms ==> transform it
    ms = SidesFromN4s(ms,s4n);  % Unclosed marked sides list
end
ms = closure(ms,s4e);  % (Closed) list of marked sides
n4ms = n4s(ms,:);
mid4ms = computeMid4s(c4n,n4ms);
% s4ms(k)==j : k-th marked side is side j (of all sides)
s4ms = SidesFromN4s(n4ms,s4n);
% nNew4s(k)==j>0 : Old side k is marked and, its center will
%                  be the new node j
% nNew4s(k)==0 : Old side k is not marked and remains
% New nodes are simply attached to the old c4n in the order
% given by ms (==> n4ms ==> mid4ms)
nNew4s = zeros(size(n4s,1),1);
nNew4s(s4ms) = (1:length(s4ms))+nr_nodes_old;
c4n = [c4n;mid4ms];  % ** Update c4n with new nodes
% nNew4e(k,m)==j>0 : Side m (1,2 or 3) of (old) element k is
%                    marked and will generate the new node j
% nNew4e(k,m)==0 : Side m of (old) element k is not marked
nNew4e = reshape(nNew4s(s4e),[],size(s4e,2));
% ** Core refinement algorithm
% Get lists of numbers of (old) elements:
% e0    : Not to be refined (all sides remain untouched)
% er    : Red (all sides are marked, thus to be refined)
% eg    : Green (only reference side is marked)
% eb/eB : Blue (reference and 2nd/3rd side are marked)
e0 = find(~any(nNew4e,2));
er = find(all(nNew4e,2));
eg = find(and(nNew4e(:,1),~any(nNew4e(:,[2 3]),2)));
eb = find(and(all(nNew4e(:,[1 2]),2),~nNew4e(:,3)));
eB = find(and(all(nNew4e(:,[1 3]),2),~nNew4e(:,2)));
%        n3        Element [n1 n2 n3] (as of a row in n4e)
%       /  \       has the sides/new nodes [s1 s2 s3] (as
%     s3    s2     given by a row in nNew4e), in the order
%    /        \    depicted here.  Each sj can be a new
%   n1 --s1-- n2   node, depending of the refinement type.
n4e = [               n4e(e0,:)                   % Untouched
    n4e(er,1) nNew4e(er,1) nNew4e(er,3)    % Red
    nNew4e(er,1)    n4e(er,2) nNew4e(er,2)
    nNew4e(er,3) nNew4e(er,2)    n4e(er,3)
    nNew4e(er,[2 3 1])
    n4e(eg,3)    n4e(eg,1) nNew4e(eg,1)    % Green
    n4e(eg,2)    n4e(eg,3) nNew4e(eg,1)
    n4e(eb,3)    n4e(eb,1) nNew4e(eb,1)    % Blue b
    nNew4e(eb,1)    n4e(eb,2) nNew4e(eb,2)
    n4e(eb,3) nNew4e(eb,1) nNew4e(eb,2)
    n4e(eB,2)    n4e(eB,3) nNew4e(eB,1)    % Blue B
    nNew4e(eB,1)    n4e(eB,3) nNew4e(eB,3)
    n4e(eB,1) nNew4e(eB,1) nNew4e(eB,3) ];
clear('e0','er','eg','eb','eB');
% Refine line
nodes_line = refine_line(ms,nr_nodes_old,sides_line,nodes_line);
% ** Refine boundary
n4sDb = refineBoundary(n4sDb,s4n,nNew4s);
end

function ms = closure(ms,s4e)
% Append more sides to list of marked sides (ms) until the
% reference sides of all elements that have at least one
% marked side are marked.
m4s = zeros(max(s4e(:)),1);
m4s(ms) = 1;  % m4s(k)==1 if side k marked, ==0 otherwise
while true
    todo4e = reshape(m4s(s4e),[],size(s4e,2));
    todo4e = and(any(todo4e(:,[2 3]),2),~todo4e(:,1));
    % Now: todo4e(k)==1 : Element k has a marked side, but its
    %      reference side is not marked, so we need to mark it
    %      also, and then check again.
    if ~any(todo4e)  % If all elements are consistent => done!
        break;
    end
    sTodo = s4e(find(todo4e),1);  % List of sides to be marked
    sTodo = unique(sTodo);
    m4s(sTodo) = 1;  % Mark those sides
end
ms = find(m4s);
end

function n4sB = refineBoundary(n4sB,s4n,nNew4s)
% Refine the boundary given by it nodes in n4sB
% s4n must be the full s4n matrix of the old grid.
if isempty(n4sB)
    return
end
% Get new (middle) node for each boundary side
nNew4sB = SidesFromN4s(n4sB,s4n);
nNew4sB = nNew4s(nNew4sB)';
%nNew4sB = reshape(nNew4sB,1,[]);  % TODO Is this needed
% nNew4sB(j)=k>0 : Boundary side j will become new node k
%           =k=0 : Boundary side j won't be modified
n4sB = n4sB';
n4sB = [n4sB(1,:) ; nNew4sB ; nNew4sB ; n4sB(2,:)];
n4sB = reshape(n4sB,[],1);
n4sB = n4sB(find(n4sB));
n4sB = reshape(n4sB,2,[]);
n4sB = n4sB';
end

function nodes_line = refine_line(ms,nr_nodes_old,sides_line,nodes_line)

% get the local indices (with respect to ms) of each side that is on the
% line and marked. This is also the index of each new node on the line
% (indexing with respect to all the new nodes). Therefore, their row index
% in c4n is local_index + size(c4n_old,1)
new_nodes_line_local = find(ismember(ms,sides_line));
new_nodes_line = new_nodes_line_local + nr_nodes_old;

% update nodes_line vector by appending the new nodes
nodes_line = [nodes_line;new_nodes_line];
end

function s = SidesFromN4s(n4s,s4n)
% Get list of sides from n4s structure
% s4n must refer to full grid!
% (This is a copy from refineUniformRed.m)
s = s4n(n4s(:,1),n4s(:,2));
s = diag(s);
s = full(s);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2009-2015
% Numerical Analysis Group
% Prof. Dr. Carsten Carstensen
% Humboldt-University
% Departement of Mathematics
% 10099 Berlin
% Germany
%
% This file is part of AFEM.
%
% AFEM is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% AFEM is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
