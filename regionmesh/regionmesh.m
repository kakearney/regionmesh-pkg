function A = regionmesh(method, src, varargin)
%REGIONMESH Create mesh for use in cost-distance function
%
% A = regionmesh('grid', src, 'xypoly', xy)
% A = regionmesh('mesh2d', src, 'xypoly', xy)
% A = regionmesh('manual', src, 'vmesh', v, 'fmesh', f);
%
% This function creates a triangular mesh filling a polygonal region.
% Additional points can also be added to the mesh.  This function is
% intended to be used with the costdistance function.
%
% Input variables:
%
%   method: Method to use to generate the background mesh
%           'grid':     Fill polygon with evenly-spaced grid nodes
%           'mesh2d':   Use mesh2d function (can be very time-consuming)
%           'manual':   User provides background mesh
%
%   src:    ns x 2 array, x and y coords extra points to add to mesh in
%           addition to the background
%
% Optional input variables (parameter/value pairs):
%
%   xypoly: np x 2 array, x and y (or lon and lat) coords of bounding
%           polygon, following Mapping Toolbox NaN-delimited conventions
%           (clockwise-oriented external contours and counterclockwise
%           intertal contours, i.e. holes) ('grid','mesh2d')
%
%   ngrdx:  Number of grid nodes in x direction ('grid') [100]
%
%   ngrdy:  Number of grid nodes in y direction ('grid') [100]
%
%   offset: Random noise can be added to the grid node positions.  This can
%           be useful so that the resulting Delaunay triangulation
%           diagonals don't always go in the same direction. Noise is
%           applied as a fraction of the grid spacing, set by this
%           parameter. ('grid') [0.1]  
%
%           For the mesh2d option, this is the distance between each added
%           point and its partner point (required to create closed polygon
%           out of the point)
%
%   vmesh:  nv x 2 array, x and y coordinates of background mesh ('manual')
%
%   fmesh:  nf x 3 array, triangulation matrix of background mesh
%           ('manual')
%
% Output variables:
%
%   A:      1 x 1 structure with the following fields
%
%           vmesh:  nvert x 2 array, x and y coordinates of final mesh
%                   nodes
%
%           fmesh:  nface x 2 array, triangulation matrix for final mesh
%
%           loc:    indices of vmesh nodes corresponding to the extra point
%                   nodes (i.e. src) 

% Copyright 2013 Kelly Kearney

%-------------------------
% Parse input
%-------------------------

p = inputParser;
p.addParameter('xypoly', []);
p.addParameter('ngrdx', 100);
p.addParameter('ngrdy', 100);
p.addParameter('offset', 0.1);
p.addParameter('plotflag', false);
p.addParameter('fmesh', []);
p.addParameter('vmesh', []);
p.parse(varargin{:});
Opt = p.Results;

% Check for necessary inputs for each method

switch method
    case {'grid', 'mesh2d'}
        if isempty(Opt.xypoly)
            error('Must supply xypoly input');
        end
    case 'manual'
        if isempty(Opt.fmesh) || isempty(Opt.vmesh)
            error('Must supply fmesh and vmesh inputs');
        end
end
        

%-------------------------
% Build mesh
%-------------------------

switch method
      
    case 'grid'
       
        % Start with coastline nodes, and build constraints

        [nx, ny] = polysplit(Opt.xypoly(:,1), Opt.xypoly(:,2));
        nnode = cellfun(@length, nx);
        coastxy = [cat(1, nx{:}) cat(1, ny{:})];
        connect = arrayfun(@(x) [(1:x-1)' (2:x)'; x 1], nnode, 'uni', 0);
        nodesum = [0; cumsum(nnode)];
        connect = cellfun(@(a,b) a+b, connect, num2cell(nodesum(1:end-1)), 'uni', 0);
        connect = cat(1, connect{:});

        % Add grid nodes

        nx = linspace(min(Opt.xypoly(:,1)), max(Opt.xypoly(:,1)), Opt.ngrdx);
        ny = linspace(min(Opt.xypoly(:,2)), max(Opt.xypoly(:,2)), Opt.ngrdy);

        [nx, ny] = meshgrid(nx, ny);

        % Slight offsets help the issue that Delaunay triangulation for a
        % grid is not unique, and will favor one diagonal over the other 

        xjig = (rand(size(nx))-0.5) * Opt.offset * (nx(1,2)-nx(1,1));
        yjig = (rand(size(ny))-0.5) * Opt.offset * (ny(2,1)-ny(1,1));

        nx = nx + xjig;
        ny = ny + yjig;

        isin = inpolygons(nx, ny, Opt.xypoly(:,1), Opt.xypoly(:,2));
        gridxy = [nx(isin) ny(isin)];

        % Add source and sink nodes, if applicable

        if ~isempty(src)
            tf = ismember(src, [coastxy; gridxy], 'rows');
            srcxy = src(~tf,:);
        else
            srcxy = zeros(0,2);
        end

        nodes = [coastxy; gridxy; srcxy];

        % Triangulate

        S = warning('off', 'MATLAB:DelaunayTri:DupPtsConsUpdatedWarnId');
        Tri = DelaunayTri(nodes, connect);
        warning(S);
        isin = inOutStatus(Tri);

        vmesh = Tri.X;
        fmesh = Tri.Triangulation(isin,:);
        
    case 'grid2' % Never mind, still a work in progress
        
        
        
        % Start with coastline nodes, and build constraints

        [nx, ny] = polysplit(Opt.xypoly(:,1), Opt.xypoly(:,2));
        nnode = cellfun(@length, nx);
        coastxy = [cat(1, nx{:}) cat(1, ny{:})];
        connect = arrayfun(@(x) [(1:x-1)' (2:x)'; x 1], nnode, 'uni', 0);
        nodesum = [0; cumsum(nnode)];
        connect = cellfun(@(a,b) a+b, connect, num2cell(nodesum(1:end-1)), 'uni', 0);
        connect = cat(1, connect{:});

        % Add grid nodes

        nx = linspace(min(Opt.xypoly(:,1)), max(Opt.xypoly(:,1)), Opt.ngrdx);
        ny = linspace(min(Opt.xypoly(:,2)), max(Opt.xypoly(:,2)), Opt.ngrdy);

        [nx, ny] = meshgrid(nx, ny);

        % Slight offsets help the issue that Delaunay triangulation for a
        % grid is not unique, and will favor one diagonal over the other 

%         xjig = (rand(size(nx))-0.5) * Opt.offset * (nx(1,2)-nx(1,1));
%         yjig = (rand(size(ny))-0.5) * Opt.offset * (ny(2,1)-ny(1,1));
% 
%         nx = nx + xjig;
%         ny = ny + yjig;

        % 8-connected neighbors
        
        [r,c] = size(nx);                         % Get the matrix size
        diagVec1 = repmat([ones(c-1,1); 0],r,1);  % Make the first diagonal vector
                                                  %   (for horizontal connections)
        diagVec1 = diagVec1(1:end-1);             % Remove the last value
        diagVec2 = [0; diagVec1(1:(c*(r-1)))];    % Make the second diagonal vector
                                                  %   (for anti-diagonal connections)
        diagVec3 = ones(c*(r-1),1);               % Make the third diagonal vector
                                                  %   (for vertical connections)
        diagVec4 = diagVec2(2:end-1);             % Make the fourth diagonal vector
                                                  %   (for diagonal connections)
        adj = diag(diagVec1,1)+...                % Add the diagonals to a zero matrix
              diag(diagVec2,c-1)+...
              diag(diagVec3,c)+...
              diag(diagVec4,c+1);
        adj = adj+adj.';                          % Add the matrix to a transposed
                                                  %   copy of itself to make it
                                                  %   symmetric
              
        isin = inpolygons(nx, ny, Opt.xypoly(:,1), Opt.xypoly(:,2));
        adj(~isin,:) = 0;
        adj(:,~isin) = 0;
        
        gridxy = [nx(isin) ny(isin)];
        adj = adj(isin,isin);

        % Add source and sink nodes, if applicable

        if ~isempty(src)
            tf = ismember(src, [coastxy; gridxy], 'rows');
            srcxy = src(~tf,:);
        else
            srcxy = zeros(0,2);
        end

        nodes = [coastxy; gridxy; srcxy];
        
        S = warning('off', 'MATLAB:DelaunayTri:DupPtsConsUpdatedWarnId');
        Tri = DelaunayTri(nodes, connect);
        warning(S);
        isin = inOutStatus(Tri);

        vmesh = Tri.X;
        fmesh = Tri.Triangulation(isin,:);
        
        
        blah
        
    case 'mesh2d'
        

        % Make "polygons" out of additional points
        
        check = inpolygons(src(:,1), src(:,2), Opt.xypoly(:,1), Opt.xypoly(:,2));
        if ~all(check)
            error('At least one point not in the bounding polygon');
        end
        
        nsrc = size(src,1);
        isin = false(nsrc,1);
        [lttmp, lntmp] = deal(zeros(nsrc,1));
        while ~all(isin)
            az = rand(sum(~isin),1)*360;
            len = Opt.offset;
            
            [lttmp(~isin), lntmp(~isin)] = reckon(src(~isin,2), src(~isin,1), Opt.offset, rand(sum(~isin),1)*360);
            
            isin = inpolygons(lntmp, lttmp, Opt.xypoly(:,1), Opt.xypoly(:,2));
        end

        xnew = [nan(nsrc,1) src(:,1) lntmp src(:,1)]';
        ynew = [nan(nsrc,1) src(:,2) lttmp src(:,2)]';
        
        
        % Convert to triangulation format
        
        xpoly = [Opt.xypoly(:,1); xnew(:)];
        ypoly = [Opt.xypoly(:,2); ynew(:)];
        
        [nx, ny] = polysplit(xpoly, ypoly);

        nnode = cellfun(@length, nx);
        node = [cat(1, nx{:}) cat(1, ny{:})];
        
        connect = arrayfun(@(x) [(1:x-1)' (2:x)'; x 1], nnode, 'uni', 0);
        nodesum = [0; cumsum(nnode)];
        connect = cellfun(@(a,b) a+b, connect, num2cell(nodesum(1:end-1)), 'uni', 0);
        connect = cat(1, connect{:});
        
        % Create mesh
        
        [vmesh, fmesh] = mesh2d(node, connect);
        
%         
%         % Add additional points
%         
%         ntri = size(fmesh,1);
% 
%         vlat = vmesh(:,2);
%         vlon = vmesh(:,1);
%         
%         plon = [vlon(fmesh) nan(ntri,1)]';
%         plat = [vlat(fmesh) nan(ntri,1)]';
%         [plon, plat] = poly2cw(plon(:), plat(:));
%         
%         ns = size(src,1);
% 
%         [isin, loc] = inpolygons(src(:,1), src(:,2), plon, plat);
%         [unqloc, idx] = aggregate(cell2mat(loc), (1:ns)');
%         ismissing = unqloc == 0;
%         if any(ismissing)
%             idxsrcmiss = idx{ismissing};
%             unqloc = unqloc(~ismissing);
%             idx = idx(~ismissing);
%         end
%         nu = length(unqloc);
%         [vnew, fnew] = deal(cell(nu,1));
%         for ii = 1:nu
%             lt = [vlat(fmesh(unqloc(ii),:)); src(idx{ii},2)];
%             ln = [vlon(fmesh(unqloc(ii),:)); src(idx{ii},1)];
%             dt = DelaunayTri(ln, lt);
% 
%             vnew{ii} = dt.X;
%             fnew{ii} = dt.Triangulation;
% 
%         end
% 
%         fmesh(unqloc,:) = [];
% 
%         for ii = 1:length(vnew)
%             [tf,loc] = ismember(vnew{ii}, vmesh, 'rows');
%             nnew = sum(~tf);
%             newidx = (1:nnew) + size(vmesh,1);
%             loc(~tf) = newidx;
% 
%             vmesh = [vmesh; vnew{ii}(~tf,:)];
%             fmesh = [fmesh; loc(fnew{ii})];
%         end
%         
%         if any(ismissing)
%             vmesh = [vmesh; src(idxsrcmiss,:)];
%         end

        
    case 'manual'
        
        vmesh = Opt.vmesh;
        fmesh = Opt.fmesh;
        
        % Add additional points
        
        ntri = size(fmesh,1);

        vlat = vmesh(:,2);
        vlon = vmesh(:,1);
        
        
        plon = [vlon(fmesh) nan(ntri,1)]';
        plat = [vlat(fmesh) nan(ntri,1)]';
        [plon, plat] = poly2cw(plon(:), plat(:));
        
        ns = size(src,1);

        [isin, loc] = inpolygons(src(:,1), src(:,2), plon, plat);
        [unqloc, idx] = aggregate(cell2mat(loc), (1:ns)');
        ismissing = unqloc == 0;
        if any(ismissing)
            idxsrcmiss = idx{ismissing};
            unqloc = unqloc(~ismissing);
            idx = idx(~ismissing);
        end
        nu = length(unqloc);
        [vnew, fnew] = deal(cell(nu,1));
        for ii = 1:nu
            lt = [vlat(fmesh(unqloc(ii),:)); src(idx{ii},2)];
            ln = [vlon(fmesh(unqloc(ii),:)); src(idx{ii},1)];
            dt = DelaunayTri(ln, lt);

            vnew{ii} = dt.X;
            fnew{ii} = dt.Triangulation;

        end

        fmesh(unqloc,:) = [];

        for ii = 1:length(vnew)
            [tf,loc] = ismember(vnew{ii}, vmesh, 'rows');
            nnew = sum(~tf);
            newidx = (1:nnew) + size(vmesh,1);
            loc(~tf) = newidx;

            vmesh = [vmesh; vnew{ii}(~tf,:)];
            fmesh = [fmesh; loc(fnew{ii})];
        end
        
        if any(ismissing)
            vmesh = [vmesh; src(idxsrcmiss,:)];
        end
         
end

%-------------------------
% Output
%-------------------------

A.fmesh = fmesh;
A.vmesh = vmesh;
if ~isempty(src)
    [tf, A.loc] = ismember(src, A.vmesh, 'rows');
else
    A.loc = [];
end


%-------------------------
% Plot, for checks
%-------------------------

if Opt.plotflag
    figure;
    triplot(fmesh, vmesh(:,1), vmesh(:,2));
    hold on;
    if ~isempty(Opt.xypoly)
        plot(Opt.xypoly(:,1), Opt.xypoly(:,2), 'k');
    end
    if ~isempty(src)
        plot(src(:,1), src(:,2), 'r.');
    end
    axis equal;
end

