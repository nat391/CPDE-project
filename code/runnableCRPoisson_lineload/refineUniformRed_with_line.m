function [c4nNew,n4eNew,n4sDbNew,n4s_lineNew] = refineUniformRed_with_line(c4n,n4e,n4s,n4sDb,n4s_line)
%% refineUniformRed - Refine every element the "red" way.
%   refineUniformRed(c4n, n4e, n4sDb, n4sNb) refines a given mesh uniformly
%       using the red refinement.
% input:  c4n               nr_nodes x 3 matrix from loadGeometry
%         n4e               nr_elems x 3 matrix from loadGeometry
%         n4s               nr_sides x 2 matrix from computeN4s
%         n4sDb             nr_dirichletSides x 2 matrix from loadGeometry
%         n4s_line          nr_sides_on_line x 2 matrix giving the two
%                           nodes for each side on the line
% output: c4nNew            updated c4n matrix with new nodes
%         n4eNew            updated n4e matrix with new elements
%         n4sDbNew          updated n4sDb matrix with new dirichlet nodes
%         n4s_lineNew       updated n4s_line matrix with new nodes on line

    %% Preliminary work.
    nrNodes = size(c4n,1);
    nrElems = size(n4e,1);
    nrSides = size(n4s,1);
    nr_sidesSlit = size(n4s_line,1);
    newNodes4s = sparse(n4s(:,1),n4s(:,2),(1:nrSides)'+ nrNodes, ...
                        nrNodes,nrNodes);
    newNodes4s = newNodes4s + newNodes4s';

    %% Compute coordinates for new nodes.
    % As every element is refined "red", there is a new node on each side.
    mid4s = computeMid4s(c4n,n4s);
    c4nNew = [c4n;mid4s];

    %% red refinement
    n4eNew = zeros(4*nrElems,3);
    for curElem = 1 : nrElems
        curNodes = n4e(curElem,:);
        curNewNodes = [newNodes4s(curNodes(1),curNodes(2));
                       newNodes4s(curNodes(2),curNodes(3));
                       newNodes4s(curNodes(3),curNodes(1));
                      ];
        % Generate new elements.
        n4eNew(4*(curElem-1)+1:4*curElem,:) = ...
                  [ curNodes(1)    curNewNodes(1) curNewNodes(3);
                    curNewNodes(1) curNodes(2)    curNewNodes(2);
                    curNewNodes(2) curNewNodes(3) curNewNodes(1);
                    curNewNodes(3) curNewNodes(2) curNodes(3);
                  ];  
    end
   
    %% refinement of  Dirichlet boundary
    n4sDbNew = zeros(2*size(n4sDb,1),2);
    for curSide = 1 : size(n4sDb,1)
        curNodes = n4sDb(curSide,:);
        curNewNodes = newNodes4s(curNodes(1),curNodes(2));
        % Generate new Dirichlet boundary sides.
        n4sDbNew(2*(curSide-1)+1:2*curSide,:) = ...
                  [ curNodes(1)    curNewNodes;
                    curNewNodes    curNodes(2);
                  ];  
    end

    %% refinement of the line
    
    n4s_lineNew = zeros(2*nr_sidesSlit,2);
    for curSide = 1 : size(n4s_line,1)
        curNodes = n4s_line(curSide,:);
        curNewNodes = newNodes4s(curNodes(1),curNodes(2));
        % Generate new line sides.
        n4s_lineNew(2*(curSide-1)+1:2*curSide,:) = ...
                  [ curNodes(1)    curNewNodes;
                    curNewNodes    curNodes(2);
                  ];  
    end

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
