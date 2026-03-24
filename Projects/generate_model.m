function model = generate_model(geom,mater,plotopt)
% Example :
% model = generate_model(geom,mater,1)



%% Generating the geometry and the mesh (spatial discretization)
[nodes, conn] = gen_mesh(geom, plotopt);


%% Generating cell's matrices using Finite Element Method (integration)
[K, M] = gen_mat(nodes, conn, mater);

% Checking for potential errors/issues
disp('Size of global stiffness matrix:'); disp(size(K));
disp('Any NaN or Inf in Kg?'); disp(any(isnan(K(:))) || any(isinf(K(:))));
disp('Trace of Mg (≈ total mass):'); disp(trace(M));

% Numerical infos
nod_edge = sum(nodes(:,1) < 1e-10); % Number of nodes on the interface
dofs_edge = 2*nod_edge; % Number of DOFs on the interface
dofs_tot=length(M); % Total number of DOFs

uL=1:dofs_edge; % Left side DOF indexes
uR=dofs_tot+1-(dofs_edge:-1:1); % Right side DOF indexes
ui=(uL(end)+1):(uR(1)-1); % Inner DOF indexes

dofs=struct('uL',uL,'ui',ui,'uR',uR);
mesh=struct('nodes',nodes,'conn',conn);

%% Final data structure, containing all unit-cell information
model=struct('geom',geom,'mater',mater,'mesh',mesh,'dofs',dofs,'K',K,'M',M)

end

% How to call elements from this data structure 'model' : 
% K_LR = model.K([model.dofs.uL, model.dofs.uR]); % Submatrix K_LR
% d = model.geom.lx; % Periodicity of the structure



function [K_global, M_global] = gen_mat(nodes, conn, mater)
    % ASSEMBLE_GLOBAL_MATRICES   Assemble global stiffness and mass matrices
    %   from element-level contributions for 2D linear triangular elements.
    %
    % Inputs:
    %   nodes     N × 2 matrix: global node coordinates [x y]
    %   conn      M × 3 matrix: element connectivity (node indices, 1-based)
    %   mater     Structure containing material properties
    % Outputs:
    %   K_global  (2N × 2N) global stiffness matrix
    %   M_global  (2N × 2N) global consistent mass matrix
    %
    % Assumptions:
    %   - Uses the function genelem() that returns [Ke, Me] for each element
    %   - Plane stress, isotropic material (as defined in genelem)
    %   - DOFs are ordered: u1, v1, u2, v2, ..., uN, vN
    
    % Total number of nodes and elements
    n_nodes = size(nodes, 1);
    n_elem  = size(conn, 1);
    
    % Total DOFs = 2 per node
    ndof = 2 * n_nodes;
    
    % Initialize global matrices (sparse is much more efficient for large models)
    K_global = sparse(ndof, ndof);
    M_global = sparse(ndof, ndof);
    
    % Loop over all elements
    for e = 1:n_elem
        
        % Get the three node indices of current element (1-based)
        node_ids = conn(e, :);
        
        % Extract coordinates of the three nodes (3×2)
        nodes_elem = nodes(node_ids, :);
        
        % Get element stiffness and mass matrices (both 6×6)
        [Ke, Me] = genelem(nodes_elem,mater);
        
        % Local-to-global DOF mapping
        % For node i: ux = 2*i-1, uy = 2*i
        dofs = zeros(6,1);
        for i = 1:3
            n = node_ids(i);
            dofs(2*i-1) = 2*n - 1;   % u_x
            dofs(2*i)   = 2*n;       % u_y
        end
        
        % Assemble into global matrices
        for i = 1:6
            for j = 1:6
                gi = dofs(i);
                gj = dofs(j);
                
                K_global(gi, gj) = K_global(gi, gj) + Ke(i,j);
                M_global(gi, gj) = M_global(gi, gj) + Me(i,j);
            end
        end
    end

    K_global=K_global*(1+mater.eta*1i);
    
    % Optional: convert to full matrices if preferred (usually not needed)
    % K_global = full(K_global);
    % M_global = full(M_global);
    
    fprintf('Global assembly completed:\n');
    fprintf('  → %d nodes, %d elements, %d DOFs\n', n_nodes, n_elem, ndof);
    fprintf('  → K_global: %d × %d  (nnz = %d)\n', ndof, ndof, nnz(K_global));
    fprintf('  → M_global: %d × %d  (nnz = %d)\n', ndof, ndof, nnz(M_global));
end



function [nodes, conn] = gen_mesh(geom, plotopt)
    % Create a 2D triangular mesh for a rectangle with a centered circular hole (void).
    % Ensures left and right side meshes are compatible (same y-coordinates).
    % Removes any elements whose centroid lies inside the hole.
    % Nodes are returned sorted: first by x, then by y.
    %
    % Inputs: geom : structure with fields
    %   lx: length of rectangle along x (horizontal)
    %   ly: length of rectangle along y (vertical)
    %   r:  radius of the circular hole
    %
    % Outputs:
    %   nodes: N x 2 matrix of node coordinates [x y] — **sorted by x then y**
    %   conn:  M x 3 matrix of triangle connectivity (node indices, updated)

    lx=geom.lx; 
    ly=geom.ly;
    r=geom.r;

    % Center of the circle
    cx = lx / 2;
    cy = ly / 2;

    % Safety check
    if r >= min(lx, ly)/2 - 1e-6
        error('Radius too large; circle must fit strictly inside the rectangle.');
    end

    % Approximate mesh size h
    min_dim = min([lx, ly, 2*r]);
    h = min_dim / 5;

    % Number of points on circle
    N_circle = max(20, round(2 * pi * r / h));
    h = 2 * pi * r / N_circle;  % refined h for circle

    % Number of divisions on sides
    Nx = max(4, round(lx / h));
    Ny = max(4, round(ly / h));

    % Define consistent spacing
    x_horz = linspace(0, lx, Nx + 1)';
    y_vert = linspace(0, ly, Ny + 1)';

    % Outer boundary points (counter-clockwise)
    bottom_P = [x_horz,          zeros(Nx+1,1)];
    right_P  = [lx*ones(Ny-1,1), y_vert(2:Ny)];
    top_P    = [x_horz(end:-1:1),ly*ones(Nx+1,1)];
    left_P   = [zeros(Ny-1,1),   y_vert(Ny:-1:2)];

    outer_ordered_P = [bottom_P; right_P; top_P; left_P];

    % Inner circle points (clockwise)
    theta = linspace(0, 2*pi, N_circle+1);
    theta(end) = [];
    circle_P = [cx + r*cos(theta)', cy + r*sin(theta)'];

    % Internal points on regular grid (outside the circle)
    xvec = h/2 : h : lx - h/2;
    yvec = h/2 : h : ly - h/2;
    [X, Y] = meshgrid(xvec, yvec);
    dist = sqrt((X - cx).^2 + (Y - cy).^2);
    keep = dist >= r + 1e-9;   % strict outside
    internal_P = [X(keep), Y(keep)];

    % All points together
    P = [outer_ordered_P; circle_P; internal_P];

    % Connectivity constraints (rigid edges)
    N_outer   = size(outer_ordered_P,1);
    N_circ_s  = N_outer + 1;
    N_circ_e  = N_outer + N_circle;

    C_outer = [(1:N_outer-1)', (2:N_outer)'; N_outer, 1];
    C_inner = [(N_circ_s:N_circ_e-1)', (N_circ_s+1:N_circ_e)'; N_circ_e, N_circ_s];
    C = [C_outer; C_inner];

    % Constrained Delaunay triangulation
    DT = delaunayTriangulation(P, C);

    % --- Remove elements inside the hole ---
    tri = DT.ConnectivityList;
    centroids = (P(tri(:,1),:) + P(tri(:,2),:) + P(tri(:,3),:)) / 3;
    dist_cent = sqrt((centroids(:,1) - cx).^2 + (centroids(:,2) - cy).^2);
    keep_tri = dist_cent >= r - 1e-9;
    conn_clean = tri(keep_tri, :);

    % Get unsorted nodes and connectivity
    nodes_unsorted = DT.Points;
    conn_unsorted  = conn_clean;

    % === SORT NODES by x then by y ===
    [~, sortIdx] = sortrows(nodes_unsorted, [1 2]);  % column 1 = x, column 2 = y

    % Apply sorting to nodes
    nodes = nodes_unsorted(sortIdx, :);

    % Remap connectivity to the new node ordering
    % Create a map: old node index → new node index
    old2new = zeros(size(nodes_unsorted,1), 1);
    old2new(sortIdx) = 1:length(sortIdx);

    % Update connectivity
    conn = old2new(conn_unsorted);


if plotopt==1
    % === Visualization ===
    figure('Color', 'w');
    triplot(conn, nodes(:,1), nodes(:,2), ...
        'Color', [0.2 0.2 0.8], ...   % nice blue lines
        'LineWidth', 1.0);
    hold on;
    theta_plot = linspace(0, 2*pi, 200);
    plot(cx + r*cos(theta_plot), cy + r*sin(theta_plot), ...
         'r--', 'LineWidth', 1.8, 'DisplayName', 'Hole boundary');
    axis equal tight;
    box on;
    grid on;
    xlabel('x'); ylabel('y');
    title(sprintf('Mesh – Rectangle [%.1f × %.1f] with circular hole r = %.2f', lx, ly, r));
    set(gca, 'FontSize', 11);
    legend('show', 'Location', 'bestoutside');
    hold off;

    % Optional: quick check that nodes are sorted
    % is_sorted = all(diff(nodes(:,1)) >= -1e-12) && ...  % non-decreasing x
    %             all(nodes(1:end-1,1) <= nodes(2:end,1) + 1e-12);
    % fprintf('Nodes sorted by x then y: %s\n', mat2str(is_sorted));
end
end





function [K, M] = genelem(nodes,mater)
    % Generate 2D triangular finite element mass and stiffness matrices
    % for an isotropic material.
    % Input:
    %   nodes: 3x2 matrix of node coordinates [x1 y1; x2 y2; x3 y3]
    %   mater     Structure containing material properties
    % Outputs:
    %   K: 6x6 stiffness matrix
    %   M: 6x6 consistent mass matrix
    % Assumptions:
    %   Plane stress condition
    
    E =  mater.E;
    nu = mater.nu;
    rho = mater.rho;
    t = mater.t;
    
    x1 = nodes(1,1); y1 = nodes(1,2);
    x2 = nodes(2,1); y2 = nodes(2,2);
    x3 = nodes(3,1); y3 = nodes(3,2);
    
    % Element area
    A = 0.5 * abs(x1*(y2 - y3) + x2*(y3 - y1) + x3*(y1 - y2));
    
    % Check for degenerate triangle
    if A < 1e-10
        error('Degenerate triangle with zero area.');
    end
    
    % Beta and gamma terms
    beta = [y2-y3, y3-y1, y1-y2];
    gamma = [x3-x2, x1-x3, x2-x1];
    
    % Strain-displacement matrix B (3x6)
    B = (1/(2*A)) * [beta(1) 0 beta(2) 0 beta(3) 0; ...
                     0 gamma(1) 0 gamma(2) 0 gamma(3); ...
                     gamma(1) beta(1) gamma(2) beta(2) gamma(3) beta(3)];
    
    % Material matrix D for plane stress (3x3)
    D = E / (1 - nu^2) * [1 nu 0; nu 1 0; 0 0 (1-nu)/2];
    
    % Stiffness matrix K = t * A * B' * D * B (6x6)
    K = t * A * B' * D * B;
    
    % Consistent mass submatrix (3x3)
    mm = (A / 12) * [2 1 1; 1 2 1; 1 1 2];
    
    % Mass matrix M (block diagonal, 6x6)
    M = rho * t * blkdiag(mm, mm);
end



