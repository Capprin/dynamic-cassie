% finds SO(3) minimum perturbation coordinates for an n-dim system

% objective function based on Ross Hatton's SO(3) coordinate optimization
% quadrature rules based on Becker Carey Oden 5.3.2

% inputs:
    % grid_points: cell array containing vector dimension(s) of points
    % A_orig: 
        % local connection tensoral object
        % ought to be a (n_dim)x(n_dim)x...x(n_dim) cell of 3x... matrices
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
    bases_at_quad = cat(2, basis_cell{:});
    % 2d cell of basis derivatives
    d_basis_cell = cellfun(@(fn) fn(quad_points), d_basis_fns,...
                           'UniformOutput', false);
    % convert to cell of derivative matrices (n_dim matrices)
    % value of derivative of basis, for each quad point's surroundings,
    % along each dimension
    d_bases_at_quad = cell(1, n_dim);
    for i = 1:n_dim
        % concat all (2^n_dim) derivative points into matrix
        % stored in cell corresp. to dimension
        d_bases_at_quad{i} = cat(2, d_basis_cell{:,i});
    end
    
    %% scale quadrature points according to problem coordinates
    % big assumption: regular grid spacing
    % preallocate space
    bases_at_quad_act = cell(n_dim, 1); % points in each dimension
    quad_deriv_points_act = cell(n_dim, n_dim); % points w.r.t. ea. dim.
    for i = 1:n_dim
        % SCALES ea. value by size of cube in ea. dim
        % assumes a regular grid, and can be done because linear basis
        bases_at_quad_act{i} = bases_at_quad * nodes(cubes(1,:), i);
        
        % do derivatives for this dim
        for j =1:n_dim
            % derivatives along all other (jth) dimensions for ith dim
            quad_deriv_points_act{i, j} = d_bases_at_quad{j} * nodes(cubes(1,:), i);
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
        d_bases_at_quad{i} = cat(2, d_basis_cell_act{:,i});
    end
    
    %% collect ingredients for calculus over the field
    % _dim concatenates ea. dimension into a single matrix
    
    % quadrature weights (for integration across dimensions)
    quad_weights_dim = repmat(quad_weights', 1, n_dim); %1x(n_points)n_dim
    
    % nabla (del, gradient)
    % derivative of each basis fn, in each direction, at each vertex
    del_dim = cat(1, d_bases_at_quad{:}); %(n_points)(n_dim) x n_vertices
    
    % gradient of basis functions
    grad_rho_dim = del_dim; % equal; all fns over field are linear w.r.t. basis fns
    
    %% integrate ingredients, simplifying math later
    
    % gradient dot product (correct)
    grad_rho_dot_del = cell(1, n_vertices); %values at each vertex
    for i = 1:n_vertices
        gradient = repmat(grad_rho_dim(:,i), 1, n_vertices) .* del_dim;
        grad_rho_dot_del{i} = quad_weights_dim * gradient;
    end
    
    % gradient scaled by rho (probably correct)
    grad_rho_dot_rho = cell(n_dim, n_vertices);
    for i = 1:n_dim
        for j = 1:n_vertices
            gradient_scaled = repmat(d_bases_at_quad{i}(:,j), 1, n_vertices) .* bases_at_quad;
            grad_rho_dot_rho{i,j} = quad_weights' * gradient_scaled;              
        end
    end
    
    % rho scaled by rho (maybe correct)
    rho_dot_rho_all = zeros(n_vertices, n_vertices, n_points);
    for i = 1:n_points
        % multiply all bases at this point, generating (2^N)x(2^N) matrix
        bases_scaled = bases_at_quad(i,:)' * bases_at_quad(i,:);
        % integrate
        rho_dot_rho_all(:,:,i) = quad_weights(i) * bases_scaled;
    end
    rho_dot_rho_mat = sum(rho_dot_rho_all, 3);
    rho_dot_rho = mat2cell(rho_dot_rho_mat, ones(1,n_vertices), n_vertices)';
    
    %% construct linear system from objective functions
    % organization:
        % rows are blocked out by equation
        % columns are blocked out by weights (x, y, z)
    % LHS:
        % 3 equations, evaluated at n_nodes points
        % assigning 3 weights (x, y, z) for each node\
    LHS = zeros(3*n_nodes);
    % RHS: vector result of an area integral
    RHS = zeros(3*n_nodes, 1);
    
    % reorganize LC (by rotational component, and by node, for later)
    % each (n_nodes)x(n_dim) organized in x,y,z pages
    A_nodes = cat(3, grid_to_columns(A_orig(1,:)),...
                     grid_to_columns(A_orig(2,:)),...
                     grid_to_columns(A_orig(3,:)));
    
    % add contributions from each node (rows of LHS, RHS)
    for n = 1:n_nodes
        cubes_containing_node = cubes(any(cubes==n, 2),:);
        
        % contributions are integral values at this node, summed over all
        % adjacent hypercubes
        for c = 1:size(cubes_containing_node, 1)
            % identify corner with node of interest
            corner_idx = cubes_containing_node(c, :) == n;
            
            % get LC at cube vertices (breaking into components)
            A_x = A_nodes(cubes_containing_node(c,:),:,1);
            A_y = A_nodes(cubes_containing_node(c,:),:,2);
            A_z = A_nodes(cubes_containing_node(c,:),:,3);
            
            % generally, need:
                % rho (and its gradient) at this node
                % rho (and its gradient) everywhere in the hypercube
            
            % build each term (everywhere in cube)
            
            % get value for corner of interest, and slot into matrix
            % should end up with a column for every node (confirm pls)
                % LHS(1:n_nodes) = [x_xcols, x_ycols, x_zcols]
                % LHS(n_nodes+1,2*n_nodes) = [y_xcols, y_ycols, y_zcols]
                % LHS(2*n_nodes+1:end) = [z_xcols, z_ycols, z_zcols]
        end
    end
end
    