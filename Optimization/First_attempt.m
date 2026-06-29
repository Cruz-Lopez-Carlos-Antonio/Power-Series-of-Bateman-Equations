% ================================================================
% File: Superposition_Bateman_complete.txt
%
% Purpose:
%   Complete MATLAB implementation of a stepwise superposition solver
%   for the last nuclide concentration X_n(t) in a linear decay chain.
%
% Main features:
%   1) Local power-series propagation on short time windows.
%   2) Restart of initial conditions at each time block.
%   3) Superposition of contributions from each nonzero component.
%   4) Cached time-independent combinatorial structure.
%   5) Garrappa-style generalized Mittag-Leffler evaluation:
%         ml(z, alpha, beta, gamma)
%
% Notes:
%   - Save this file as: Superposition_Bateman.m
%   - Make sure Garrappa's ml(...) is on the MATLAB path.
%   - The current example targets the chain
%       [lambda2, lambda2, lambda3, lambda3, lambda3, ...
%        lambda4, lambda4, lambda4, lambda4]
%     and computes X_9(t).
% ================================================================

clear; clc;

% ------------------------------------------------
% USER INPUT
% ------------------------------------------------
% Replace these half-lives with the benchmark values you want to use.
half_lives = [ ...
    2, 2, ...          % lambda_2, lambda_2
    3, 3, 3, ...       % lambda_3, lambda_3, lambda_3
    4, 4, 4, 4 ...     % lambda_4, lambda_4, lambda_4, lambda_4
];

lambda = log(2) ./ half_lives;

% Global time vector where X_n(t) is requested
Time_vector = [1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100];

% Initial condition: only the first element of the full chain is populated
x10 = 6.023e23;

% Local propagation parameters
dt_block  = 1;                  % window length in seconds
Mmax      = 5;                 % local truncation level
tol       = 1e-10;              % optional stopping tolerance
alpha_eps = 0.999999999999;     % Garrappa workaround

% ------------------------------------------------
% SOLVE X_n(t) BY STEPWISE SUPERPOSITION
% ------------------------------------------------
[Xn_series, history] = bateman_superposition_solver( ...
    lambda, x10, Time_vector, dt_block, Mmax, tol, alpha_eps);

disp('Computed X_n(t):')
disp(table(Time_vector(:), Xn_series(:), ...
    'VariableNames', {'t','Xn_series'}))


% ================================================================
% MAIN DRIVER
% ================================================================
function [Xn_out, history] = bateman_superposition_solver( ...
    lambda, x10, Time_vector, dt_block, Mmax, tol, alpha_eps)
% BATEMAN_SUPERPOSITION_SOLVER
%
% Computes the last nuclide concentration X_n(t) by dividing time into
% windows of length dt_block and restarting the solution at each block.
%
% At each restart, the current state vector
%     [X1, X2, ..., Xn]
% is decomposed by superposition into elementary states:
%
%   state 1: [X1, 0,  0, ..., 0]
%   state 2: [0,  X2, 0, ..., 0]
%   ...
%   state n: [0,  0,  0, ..., Xn]
%
% and each elementary state is propagated only through its corresponding
% subchain:
%
%   j -> j+1 -> ... -> n
%
% The contributions are then added to reconstruct the full state
% at the end of the block.
%
% OUTPUT:
%   Xn_out  : approximation of X_n(t) at the requested times
%   history : structure storing block states and diagnostics

    lambda = lambda(:).';
    Time_vector = Time_vector(:).';

    n = numel(lambda);

    if n < 1
        error('The lambda vector must contain at least one element.');
    end

    if isempty(Time_vector)
        Xn_out = [];
        history = struct();
        return;
    end

    tmax = max(Time_vector);

    % Block edges: [0, dt_block, 2*dt_block, ...]
    block_edges = 0:dt_block:tmax;
    if block_edges(end) < tmax
        block_edges = [block_edges, tmax];
    end

    % Full state at the beginning of each block
    state0 = zeros(1, n);
    state0(1) = x10;

    % Output
    Xn_out = zeros(size(Time_vector));

    % Diagnostics
    history = struct();
    history.block_edges = block_edges;
    history.block_states = cell(numel(block_edges)-1, 2); % {start_state, end_state}

    % ------------------------------------------------------------
    % March block by block
    % ------------------------------------------------------------
    for b = 1:(numel(block_edges)-1)

        t_left  = block_edges(b);
        t_right = block_edges(b+1);

        % Requested times that fall inside the current block
        mask = (Time_vector > t_left) & (Time_vector <= t_right);
        t_req = Time_vector(mask);

        % Local times in the current block
        tau_req = t_req - t_left;

        % Always include tau = block length to update the block state
        tau_end = t_right - t_left;

        history.block_states{b,1} = state0;

        % Propagate the whole state by superposition
        [state_req, state_end] = propagate_block_by_superposition( ...
            lambda, state0, tau_req, tau_end, Mmax, tol, alpha_eps);

        % Save X_n(t) on requested points
        if ~isempty(t_req)
            Xn_out(mask) = state_req(:, n).';
        end

        % Store end state and move to next block
        history.block_states{b,2} = state_end;
        state0 = state_end;
    end
end


% ================================================================
% BLOCK PROPAGATION BY SUPERPOSITION
% ================================================================
function [state_req, state_end] = propagate_block_by_superposition( ...
    lambda, state0, tau_req, tau_end, Mmax, tol, alpha_eps)
% PROPAGATE_BLOCK_BY_SUPERPOSITION
%
% Given an initial state vector at the beginning of a block,
% decomposes it into elementary populated components and propagates
% each one through its corresponding subchain. The results are summed.
%
% INPUT:
%   lambda   : full chain decay constants [lambda1, ..., lambdan]
%   state0   : full initial state at current block
%   tau_req  : local times inside the block where the state is requested
%   tau_end  : final local time of the block
%
% OUTPUT:
%   state_req : matrix (#requested_times) x n
%   state_end : row vector 1 x n

    n = numel(lambda);
    tau_req = tau_req(:).';  % row vector

    % Allocate
    state_req = zeros(numel(tau_req), n);
    state_end = zeros(1, n);

    % ------------------------------------------------------------
    % Superposition over each nonzero initial component
    % ------------------------------------------------------------
    for j = 1:n

        xj0 = state0(j);

        if xj0 == 0
            continue;
        end

        % Subchain starts at node j and ends at node n
        lambda_sub = lambda(j:end);

        % --------------------------------------------------------
        % Requested local times inside the block
        % --------------------------------------------------------
        if ~isempty(tau_req)
            sub_req = local_subchain_solution( ...
                lambda_sub, xj0, tau_req, Mmax, tol, alpha_eps);
            % sub_req is expected to be (#times) x (n-j+1)

            state_req(:, j:end) = state_req(:, j:end) + sub_req;
        end

        % --------------------------------------------------------
        % End-of-block state
        % --------------------------------------------------------
        sub_end = local_subchain_solution( ...
            lambda_sub, xj0, tau_end, Mmax, tol, alpha_eps);
        % sub_end is expected to be 1 x (n-j+1)

        state_end(j:end) = state_end(j:end) + sub_end;
    end
end


% ================================================================
% LOCAL SUBCHAIN SOLUTION
% ================================================================
function sub_state = local_subchain_solution( ...
    lambda_sub, xj0, tau, Mmax, tol, alpha_eps)
% LOCAL_SUBCHAIN_SOLUTION
%
% Given a subchain
%     X_j -> X_{j+1} -> ... -> X_n
% with only X_j(0)=xj0 nonzero at the local restart time,
% returns the concentrations of all members of that subchain
% at the requested local time(s) tau.
%
% OUTPUT SHAPE:
%   - if tau is a vector of length Nt:
%         sub_state is Nt x length(lambda_sub)
%   - if tau is scalar:
%         sub_state is 1 x length(lambda_sub)
%
% STRATEGY:
%   For each target node r = 1,...,length(lambda_sub),
%   compute the local concentration of the last node of the prefix
%   chain lambda_sub(1:r). This gives the concentration of node r
%   of the subchain.

    lambda_sub = lambda_sub(:).';
    tau = tau(:).';   % row vector

    nr = numel(lambda_sub);
    Nt = numel(tau);

    sub_state = zeros(Nt, nr);

    for r = 1:nr
        lambda_prefix = lambda_sub(1:r);

        % Concentration of the r-th node of the subchain equals
        % the last-node concentration of the prefix chain.
        [xr, ~] = local_series_last_node( ...
            lambda_prefix, xj0, tau, Mmax, tol, alpha_eps);

        sub_state(:, r) = xr(:);
    end
end


% ================================================================
% HELPER 1: polynomial coefficients a0,...,an
% ================================================================
function a = bateman_poly_coeffs(lambda)
% BATEMAN_POLY_COEFFS
%
% Computes coefficients of
%   prod_i (s + lambda_i) = a0 s^n + a1 s^(n-1) + ... + an
%
% OUTPUT ORDERING:
%   a = [a0, a1, ..., an]

    lambda = lambda(:).';
    n = numel(lambda);

    a = zeros(1, n + 1);
    a(1) = 1;   % a0

    for j = 1:n
        for k = j:-1:1
            a(k + 1) = a(k + 1) + lambda(j) * a(k);
        end
    end
end


% ================================================================
% CACHED LOCAL POWER-SERIES PROPAGATOR
% ================================================================
function [xn, info] = local_series_last_node( ...
    lambda, x10, t, Mmax, tol, alpha_eps)
% LOCAL_SERIES_LAST_NODE
%
% Improved version of the local last-node propagator.
%
% Instead of regenerating all compositions and coefficients for each
% call and for each time value, this routine first builds a CACHE
% containing all time-independent data:
%
%   - compositions K for each m,
%   - beta(k),
%   - gamma(k),
%   - coefficient
%       (-1)^|k| * c^k / k!
%
% After that, the time loop only evaluates the Mittag-Leffler block
% and the power t^gamma.
%
% INPUTS
%   lambda    : [lambda1, ..., lambdan]
%   x10       : initial concentration of first node of local chain
%   t         : scalar or vector of local times
%   Mmax      : maximum truncation level
%   tol       : optional stopping tolerance
%   alpha_eps : value close to 1 used in Garrappa's ml()
%
% OUTPUTS
%   xn        : approximation to X_n(t)
%   info      : diagnostics

    if nargin < 5 || isempty(tol)
        tol = [];
    end

    lambda = lambda(:).';
    t = t(:).';   % row vector

    n = numel(lambda);

    if n < 1
        error('The lambda vector must contain at least one element.');
    end

    % ------------------------------------------------------------
    % Trivial one-node case
    % ------------------------------------------------------------
    if n == 1
        xn = x10 * exp(-lambda(1) * t);

        info = struct();
        info.a = [1, lambda(1)];
        info.mu = lambda(1);
        info.c = [];
        info.prefactor = x10;
        info.n = 1;
        info.alpha_eps = alpha_eps;
        info.cache = [];
        return;
    end

    % ------------------------------------------------------------
    % Build all time-independent data only once
    % ------------------------------------------------------------
    cache = build_series_cache(lambda, x10, Mmax);

    z = -cache.mu .* t;
    xn = zeros(size(t));

    level_contrib = zeros(numel(cache.levels), numel(t));

    % ------------------------------------------------------------
    % Sum level by level
    % ------------------------------------------------------------
    for ell = 1:numel(cache.levels)
        level = cache.levels{ell};

        Xm = zeros(size(t));

        for q = 1:numel(level.coeff)
            coeff_q = level.coeff(q);
            beta_q  = level.beta(q);
            gamma_q = level.gamma(q);
            m_q     = level.m;

            % Mittag-Leffler derivative via Garrappa
            Eder = factorial(m_q) .* ml(z, alpha_eps, beta_q + alpha_eps*m_q, m_q + 1);

            term = coeff_q .* (t .^ gamma_q) .* Eder;
            Xm = Xm + term;
        end

        level_contrib(ell, :) = cache.prefactor .* Xm;
        xn = xn + level_contrib(ell, :);

        % Optional early stopping
        if ~isempty(tol)
            denom = max(1, max(abs(xn)));
            if max(abs(level_contrib(ell, :))) / denom < tol
                level_contrib = level_contrib(1:ell, :);
                break;
            end
        end
    end

    info = struct();
    info.a = cache.a;
    info.mu = cache.mu;
    info.c = cache.c;
    info.prefactor = cache.prefactor;
    info.n = n;
    info.alpha_eps = alpha_eps;
    info.cache = cache;
    info.level_contrib = level_contrib;
end


% ================================================================
% NEW HELPER: build cache for time-independent series data
% ================================================================
function cache = build_series_cache(lambda, x10, Mmax)
% BUILD_SERIES_CACHE
%
% Precomputes all time-independent objects needed by the local series
% propagator for the last node of a chain.
%
% OUTPUT:
%   cache.a         = polynomial coefficients [a0, ..., an]
%   cache.mu        = a1/a0
%   cache.c         = [an/a0, a_{n-1}/a0, ..., a2/a0]
%   cache.prefactor = x10 * prod(lambda(1:n-1))/a0
%   cache.levels    = cell array, one entry per m = 0,...,Mmax
%
% Each cache.levels{m+1} contains:
%   .m      = truncation level
%   .K      = all compositions for this m
%   .beta   = vector beta(k)
%   .gamma  = vector gamma(k)
%   .coeff  = vector (-1)^|k| c^k / k!
%
% Notes:
%   The time-independent quantities are computed once here, then reused
%   for every evaluation time in the local block.

    lambda = lambda(:).';
    n = numel(lambda);

    a = bateman_poly_coeffs(lambda);
    a0 = a(1);
    a1 = a(2);

    mu = a1 / a0;
    c  = a(end:-1:3) / a0;
    prefactor = x10 * prod(lambda(1:n-1)) / a0;

    % Weights for beta(k), gamma(k)
    w_beta  = (n-1):-1:1;
    w_gamma = n:-1:2;

    levels = cell(Mmax + 1, 1);

    p = n - 1;   % number of Diophantine variables

    for m = 0:Mmax
        K = compositions_nonnegative(m, p);   % rows = solutions
        Nr = size(K, 1);

        beta_vec  = zeros(Nr, 1);
        gamma_vec = zeros(Nr, 1);
        coeff_vec = zeros(Nr, 1);

        for r = 1:Nr
            kvec = K(r, :);

            beta_vec(r)  = n + sum(w_beta  .* kvec);
            gamma_vec(r) = (n - 1) + sum(w_gamma .* kvec);

            % coefficient = (-1)^m * c^k / k!
            log_kfact = sum(gammaln(kvec + 1));

            coeff_prod = 1;
            zero_hit = false;

            for q = 1:numel(c)
                if c(q) == 0
                    if kvec(q) > 0
                        zero_hit = true;
                        break;
                    end
                else
                    coeff_prod = coeff_prod * c(q)^kvec(q);
                end
            end

            if zero_hit
                coeff_vec(r) = 0;
            else
                coeff_vec(r) = (-1)^m * coeff_prod * exp(-log_kfact);
            end
        end

        level = struct();
        level.m = m;
        level.K = K;
        level.beta = beta_vec;
        level.gamma = gamma_vec;
        level.coeff = coeff_vec;

        levels{m + 1} = level;
    end

    cache = struct();
    cache.a = a;
    cache.mu = mu;
    cache.c = c;
    cache.prefactor = prefactor;
    cache.levels = levels;
end


% ================================================================
% FASTER COMPOSITIONS GENERATOR
% ================================================================
function K = compositions_nonnegative(m, p)
% COMPOSITIONS_NONNEGATIVE
%
% Faster generator of all nonnegative integer p-tuples summing to m.
%
% Strategy:
%   A p-tuple [x1,...,xp] with x1+...+xp=m is determined by choosing
%   p-1 separator positions among the integers 1,...,m+p-1
%   ("stars and bars").
%
% Number of rows:
%   nchoosek(m+p-1, p-1)
%
% Output:
%   K : matrix with one composition per row

    if p < 1
        error('The number of parts p must be at least 1.');
    end

    if p == 1
        K = m;
        return;
    end

    if m == 0
        K = zeros(1, p);
        return;
    end

    bars = nchoosek(1:(m + p - 1), p - 1);
    Nb = size(bars, 1);

    K = zeros(Nb, p);

    for r = 1:Nb
        b = bars(r, :);
        y = [0, b, m + p];
        K(r, :) = diff(y) - 1;
    end
end
