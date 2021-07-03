% finds SO(3) minimum perturbation coordinates for an n-dim system
% based heavily on Ross Hatton's reference point optimization

% inputs:
    % grid_points: cell array containing vector dimension(s) of points
    % A_orig: 
        % local connection tensoral object, of corresponding
        % dimensionality
% outputs:
    % beta: n-dim coordinate transform
    % A_opt:
        % local connection tensoral object, in optimal coordinates

function [beta, A_opt] = optimize_so3(grid_points, A_orig)
    % dimensionality
    n_dim = length(grid_points);
    
    %% hypercube mesh setup
    % generate hypercube mesh to work over
    grid = cell(1, n_dim);
    [grid{:}] = ndgrid(grid_points{:});
    [nodes, cubes] = hypercube_mesh(grid);
    % generate shape functions
    [basis_fns, d_basis_fns] = hypercube_element_shape_functions(n_dim);
    % generate points, weights for quadrature integration
    [quad_points, quad_weights] = hypercube_quadrature_init(n_dim);
    % variables for clarity
    n_nodes = size(nodes, 1);
    n_points = size(quad_points{1},1);
    n_vertices = size(cubes, 2);
    
    %% compute values of basis fns/derivatives at quad vertices
    % matrix of basis values at quadpoints (n_points rows of 2^n_dim columns)
    % for ea. quad point, value of surrounding points in a hypercube
    basis_cell = cellfun(@(fn) fn(quad_points),basis_fns,...
                         'UniformOutput', false);
    basis_at_quad = cat(2, basis_cell{:});
    % 2d cell of basis derivatives
    d_basis_cell = cellfun(@(fn) fn(quad_points), d_basis_fns,...
                           'UniformOutput', false);
    % convert to cell of derivative matrices (n_dim matrices)
    % value of derivative of basis, for each quad point's surroundings,
    % along each dimension
    d_basis_at_quad = cell(1, n_dim);
    for i = 1:n_dim
        % concat all (2^n_dim) derivative points into matrix
        % stored in cell corresp. to dimension
        d_basis_at_quad{i} = cat(2, d_basis_cell{:,i});
    end
    
    %% scale quadrature points according to problem coordinates
    % big assumption: regular grid spacing
    % preallocate space
    quad_points_act = cell(n_dim, 1); % points in each dimension
    quad_deriv_points_act = cell(n_dim, n_dim); % points w.r.t. ea. dim.
    for i = 1:n_dim
        % SCALES ea. quadpoint by size of cube in ea. dim
        % assumes a regular grid, and can be done because linear basis
        quad_points_act{i} = basis_at_quad * nodes(cubes(1,:), i);
        
        % do derivatives for this dim
        for j =1:n_dim
            % derivatives along all other (jth) dimensions for ith dim
            quad_deriv_points_act{i, j} = d_basis_at_quad{j} * nodes(cubes(1,:), i);
        end
    end
    
    %% scale basis function derivatives using the Jacobian
    % get Jacobian from first quadpoint (invariant across regular grid)
    J_all = celltensorconvert(quad_deriv_points_act); %nxn cell -> cell of nxn
    J = J_all{1};
    % scale derivatives w.r.t. problem coordinates
    d_basis_mat = cell2mat(d_basis_cell);
    d_basis_mat_act = (J\d_basis_mat')'; %effectively "divides" by J
    % split into (n_bases x n_dim) 2d cell of n_quad_points vectors
    d_basis_cell_act = mat2cell(d_basis_mat_act,...
                                n_points*ones(n_vertices,1),...
                                ones(1,n_dim));
    % finally, construct rescaled cell of derivative matrices
    for i = 1:n_dim %loop along ea. dimension
        d_basis_at_quad{i} = cat(2, d_basis_cell_act{:,i});
    end
    
    %% collect ingredients for calculus over the field
    % _dim concatenates ea. dimension into a single matrix
    
    % quadrature weights (for integration across dimensions)
    quad_weights_dim = repmat(quad_weights', 1, n_dim); %1x(n_points)n_dim
    
    % nabla (del, gradient)
    % derivative of each basis fn, in each direction, at each vertex
    del_dim = cat(1, d_basis_at_quad{:}); %(n_points)n_dim x n_vertices
    
    % rho
    rho = basis_at_quad;
    
    % gradient of basis functions
    grad_rho_dim = del_dim; % equal; all fns over field are linear w.r.t. basis fns
    
    % grad_rho dot nabla
    grad_rho_dot_del_dim = cell(1, n_vertices); %gradient at each vertex
    for i = 1:n_vertices
        grad_rho_dot_del_dim{i} = repmat(grad_rho_dim(:,i), 1, n_vertices) .* del_dim;
    end
    
    %% construct linear system from objective functions
    % LHS:
        % 3 equations, evaluated at n_nodes points
        % assigning 3 weights (x, y, z) for each node\
    LHS = zeros(3*n_nodes);
    % RHS: vector result of an area integral
    RHS = zeros(3*n_nodes, 1);
    
    % organization:
        % rows are blocked out by equation
        % columns are blocked out by weights (x, y, z)
    
    % add contributions from each node (rows of LHS, RHS)
    for n = 1:n_nodes
        cubes_containing_node = cubes(any(cubes==n, 2),:);
        
        % contributions are integral values at this node, summed over all
        % adjacent hypercubes
        for c = 1:size(cubes_containing_node, 1)
            % TODO: codify objective function
        end
    end
end
    