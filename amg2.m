clc; clear; close all;

% --- Setup Problem ---
n = 30; % Grid resolution (System size N = n^2 = 900)
A1 = gallery('poisson', n);
b = ones(size(A1, 1), 1);
tol = 1e-6;
max_cycles = 100;

% --- Setup Multigrid Hierarchy & Solve ---
AMG = amg_setup(A1);
[x, res] = amg_solve(AMG, b, max_cycles, tol);

% --- Visualization ---
X = reshape(x, n, n);
figure('Name', 'AMG Solution Contour');
contourf(X, 20);
colorbar;
title('AMG Solution Contours');

figure('Name', 'Convergence History');
semilogy(res, '-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
title('AMG Residual History');
xlabel('Cycle');
ylabel('||r|| / ||b||');
grid on;

%% ==========================================
%% --- AMG FUNCTIONS ---
%% ==========================================

function AMG = amg_setup(A)
    theta = 0.25;
    AMG.A{1} = A;
    level = 1;
    
    % Build levels until small coarse size or max level reached
    while size(AMG.A{level}, 1) > 20 && level < 8
        Alevel = AMG.A{level};
        
        % Classical C/F Coarsening
        [S, ST] = connections(Alevel, theta);
        [C, F]   = coarsening(S, ST, size(Alevel, 1));
        
        if isempty(C) || isempty(F)
            break; 
        end
        
        % Construct Prolongation P and Galerkin Coarse Matrix Ac
        P  = interpolate(Alevel, S, C, F);
        Ac = P' * Alevel * P;
        
        fprintf('Level %d: Size = %d -> %d coarse nodes\n', level, size(Alevel, 1), length(C));
        
        AMG.P{level}   = P;
        AMG.A{level+1} = Ac;
        level = level + 1;
    end
    AMG.levels = level;
end

function [S, ST] = connections(A, theta)
    n = size(A, 1);
    S = cell(n, 1);
    ST = cell(n, 1);
    
    for i = 1:n
        S{i} = []; ST{i} = [];
    end
    
    % Identify strong off-diagonal connections
    for i = 1:n
        [~, cols, vals] = find(A(i, :));
        mask = (cols ~= i) & (vals < 0);
        cols = cols(mask);
        vals = vals(mask);
        
        if isempty(vals), continue; end
        
        strength = -vals;
        alpha = max(strength);
        strong_idx = cols(strength >= theta * alpha);
        
        S{i} = strong_idx;
        for j = strong_idx
            ST{j} = [ST{j}, i];
        end
    end
end

function [C, F] = coarsening(S, ST, n)
    state = zeros(n, 1); % 0: undecided, 1: C, -1: F
    lambda = zeros(n, 1);
    
    for i = 1:n
        lambda(i) = length(ST{i});
    end
    
    % Pass 1: Primary C/F Assignment
    while any(state == 0)
        score = lambda;
        score(state ~= 0) = -Inf;
        [~, i] = max(score);
        state(i) = 1; % C-point
        
        for j = ST{i}
            if state(j) == 0
                state(j) = -1; % F-point
                for k = ST{j}
                    if state(k) == 0
                        lambda(k) = lambda(k) + 1;
                    end
                end
            end
        end
    end
    
    % Pass 2: Ensure F-F connections share C-points
    F_nodes = find(state == -1)';
    for i = F_nodes
        for j = intersect(S{i}, F_nodes)
            if isempty(intersect(S{i}, intersect(S{j}, find(state == 1))))
                state(i) = 1;
                break;
            end
        end
    end
    
    C = find(state == 1);
    F = find(state == -1);
end

function P = interpolate(A, S, C, F)
    n = size(A, 1);
    nc = length(C);
    
    map = zeros(n, 1);
    map(C) = 1:nc;
    
    % Preallocate Sparse Triplet Lists
    i_idx = [];
    j_idx = [];
    v_val = [];
    
    % Identity for C-points
    for k = 1:nc
        i_idx(end+1) = C(k); 
        j_idx(end+1) = k;   
        v_val(end+1) = 1.0;  
    end
    
    % Standard Interpolation for F-points
    for i = F(:)'
        Ci = intersect(S{i}, C);
        Fi = intersect(S{i}, F);
        
        if isempty(Ci), continue; end
        
        off_diag = find(A(i, :));
        off_diag(off_diag == i) = [];
        weak_nodes = setdiff(off_diag, S{i});
        
        down = A(i, i) + sum(A(i, weak_nodes));
        if abs(down) < 1e-12, down = A(i, i); end
        
        for j = Ci(:)'
            up = A(i, j);
            for m = Fi(:)'
                k_C = intersect(S{m}, C);
                denom_2 = sum(A(m, k_C));
                if abs(denom_2) > 1e-12
                    up = up + A(i, m) * (A(m, j) / denom_2);
                end
            end
            
            i_idx(end+1) = i;       
            j_idx(end+1) = map(j);   
            v_val(end+1) = -up / down;
        end
    end
    
    P = sparse(i_idx, j_idx, v_val, n, nc);
end

function x = gauss_seidel(A, x, b, iter)
    DL = tril(A); 
    for i = 1:iter
        x = x + DL \ (b - A * x);
    end
end

function x = V_cycle(AMG, level, x, b)
    A = AMG.A{level};
    
    % Exact Direct Solve at Coarsest Level
    if level == AMG.levels
        x = A \ b;
        return;
    end
    
    % Pre-smoothing
    x = gauss_seidel(A, x, b, 2);
    
    % Residual Computation & Restriction
    r = b - A * x;
    r_c = AMG.P{level}' * r;
    
    % Coarse Grid Correction
    e_c = zeros(size(r_c));
    e_c = V_cycle(AMG, level + 1, e_c, r_c);
    
    % Prolongation & Error Correction
    x = x + AMG.P{level} * e_c;
    
    % Post-smoothing
    x = gauss_seidel(A, x, b, 2);
end

function [x, res] = amg_solve(AMG, b, cycles, tol)
    n = size(AMG.A{1}, 1);
    x = zeros(n, 1);
    bnorm = norm(b);
    res = zeros(cycles, 1);
    
    for k = 1:cycles
        x = V_cycle(AMG, 1, x, b);
        r = b - AMG.A{1} * x;
        rnorm = norm(r);
        res(k) = rnorm / bnorm;
        
        fprintf('Cycle %3d | Relative Residual: %.3e\n', k, res(k));
        
        if res(k) < tol
            fprintf('Converged at cycle %d!\n', k);
            res = res(1:k);
            break;
        end
    end
end