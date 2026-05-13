% Requeriments
%   - Garrappa's ml(...) is on the MATLAB path
% ================================================================

clear; clc;
half_lives = [2,2,3,3,3,4,4,4,4];
lambda = log(2) ./ half_lives;
Time_vector = [1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100];
x10 = 6.023e23;
dt_block  = 1;
Mmax      = 5;
tol       = 1e-10;
alpha_eps = 0.999999999999;
verbose   = false;

% Cache for vector/scalar Mittag-Leffler evaluations
MLCache = containers.Map('KeyType','char','ValueType','any');

% Cache for series metadata associated with each lambda-prefix
SeriesCache = containers.Map('KeyType','char','ValueType','any');

% Precompute factorials once
FactCache = factorial(0:Mmax);

tic
[Xn_series, history, MLCache, SeriesCache] = bateman_superposition_solver( ...
    lambda, x10, Time_vector, dt_block, Mmax, tol, alpha_eps, ...
    MLCache, SeriesCache, FactCache, verbose);

disp('Computed X_n(t):')
disp(table(Time_vector(:), Xn_series(:), 'VariableNames', {'t','Xn_series'}))

filename = 'Bateman_superposition_results_optimized_autonomous.txt';
fid = fopen(filename,'w');
fprintf(fid,'t\tXn_series\n');
for i = 1:length(Time_vector)
    fprintf(fid,'%d\t%.16e\n', Time_vector(i), Xn_series(i));
end
fclose(fid);

fprintf('Results written to %s\n', filename);
fprintf('ML cache size     = %d entries\n', MLCache.Count);
fprintf('Series cache size = %d entries\n', SeriesCache.Count);
toc

function [Xn_out, history, MLCache, SeriesCache] = bateman_superposition_solver( ...
    lambda, x10, Time_vector, dt_block, Mmax, tol, alpha_eps, ...
    MLCache, SeriesCache, FactCache, verbose)

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
    block_edges = 0:dt_block:tmax;
    if block_edges(end) < tmax
        block_edges = [block_edges, tmax];
    end

    state0 = zeros(1, n);
    state0(1) = x10;
    Xn_out = zeros(size(Time_vector));

    history = struct();
    history.block_edges = block_edges;
    history.block_states = cell(numel(block_edges)-1, 2);

    for b = 1:(numel(block_edges)-1)
        if verbose
            fprintf('Processing time block %d / %d\n', b, numel(block_edges)-1);
        end

        t_left  = block_edges(b);
        t_right = block_edges(b+1);

        mask = (Time_vector > t_left) & (Time_vector <= t_right);
        t_req = Time_vector(mask);

        tau_req = t_req - t_left;
        tau_end = t_right - t_left;

        history.block_states{b,1} = state0;

        [state_req, state_end, MLCache, SeriesCache] = propagate_block_by_superposition( ...
            lambda, state0, tau_req, tau_end, Mmax, tol, alpha_eps, ...
            MLCache, SeriesCache, FactCache, verbose);

        if ~isempty(t_req)
            Xn_out(mask) = state_req(:, n).';
        end

        history.block_states{b,2} = state_end;
        state0 = state_end;
    end
end


function [state_req, state_end, MLCache, SeriesCache] = propagate_block_by_superposition( ...
    lambda, state0, tau_req, tau_end, Mmax, tol, alpha_eps, ...
    MLCache, SeriesCache, FactCache, verbose)

    n = numel(lambda);
    tau_req = tau_req(:).';

    state_req = zeros(numel(tau_req), n);
    state_end = zeros(1, n);

    for j = 1:n
        xj0 = state0(j);
        if xj0 == 0
            continue;
        end

        lambda_sub = lambda(j:end);

        if ~isempty(tau_req)
            [sub_req, MLCache, SeriesCache] = local_subchain_solution( ...
                lambda_sub, xj0, tau_req, Mmax, tol, alpha_eps, ...
                MLCache, SeriesCache, FactCache, verbose);
            state_req(:, j:end) = state_req(:, j:end) + sub_req;
        end

        [sub_end, MLCache, SeriesCache] = local_subchain_solution( ...
            lambda_sub, xj0, tau_end, Mmax, tol, alpha_eps, ...
            MLCache, SeriesCache, FactCache, verbose);
        state_end(j:end) = state_end(j:end) + sub_end;
    end
end


function [sub_state, MLCache, SeriesCache] = local_subchain_solution( ...
    lambda_sub, xj0, tau, Mmax, tol, alpha_eps, ...
    MLCache, SeriesCache, FactCache, verbose)

    lambda_sub = lambda_sub(:).';
    tau = tau(:).';

    nr = numel(lambda_sub);
    Nt = numel(tau);
    sub_state = zeros(Nt, nr);

    for r = 1:nr
        lambda_prefix = lambda_sub(1:r);
        [xr, ~, MLCache, SeriesCache] = local_series_last_node( ...
            lambda_prefix, xj0, tau, Mmax, tol, alpha_eps, ...
            MLCache, SeriesCache, FactCache, verbose);
        sub_state(:, r) = xr(:);
    end
end


function a = bateman_poly_coeffs(lambda)
    lambda = lambda(:).';
    n = numel(lambda);

    a = zeros(1, n + 1);
    a(1) = 1;

    for j = 1:n
        for k = j:-1:1
            a(k + 1) = a(k + 1) + lambda(j) * a(k);
        end
    end
end


function [xn, info, MLCache, SeriesCache] = local_series_last_node( ...
    lambda, x10, t, Mmax, tol, alpha_eps, ...
    MLCache, SeriesCache, FactCache, verbose)

    if nargin < 5 || isempty(tol)
        tol = [];
    end

    lambda = lambda(:).';
    t = t(:).';
    n = numel(lambda);

    if n < 1
        error('The lambda vector must contain at least one element.');
    end

    if n == 1
        xn = x10 * exp(-lambda(1) * t);
        info = struct();
        info.a = [1, lambda(1)];
        info.mu = lambda(1);
        info.c = [];
        info.base_prefactor = 1;
        info.prefactor = x10;
        info.n = 1;
        info.alpha_eps = alpha_eps;
        info.cache = [];
        return;
    end

    [cache, SeriesCache] = get_or_build_series_cache(lambda, Mmax, SeriesCache);

    z = -cache.mu .* t;
    xn = zeros(size(t));
    level_contrib = zeros(numel(cache.levels), numel(t));
    prefactor = x10 * cache.base_prefactor;

    for ell = 1:numel(cache.levels)
        if verbose
            fprintf('   series level %d / %d\n', ell, numel(cache.levels));
        end

        level = cache.levels{ell};
        if isempty(level.coeff)
            continue;
        end

        m_q = level.m;
        fact_m = FactCache(m_q + 1);
        beta_eff_vec = level.beta + alpha_eps * m_q;

        % --- Get Mittag-Leffler values only for UNIQUE beta-values ---
        [ubeta, ~, ibeta] = unique(beta_eff_vec);
        Ebeta = cell(numel(ubeta), 1);
        for ib = 1:numel(ubeta)
            [Ebeta{ib}, MLCache] = get_cached_ml_eval( ...
                z, alpha_eps, ubeta(ib), m_q + 1, fact_m, MLCache);
        end

        % --- Build powers only for UNIQUE gamma-values ---
        [ugamma, ~, igamma] = unique(level.gamma);
        Tgamma = cell(numel(ugamma), 1);
        for ig = 1:numel(ugamma)
            Tgamma{ig} = t .^ ugamma(ig);
        end

        Xm = zeros(size(t));
        for q = 1:numel(level.coeff)
            Eder = Ebeta{ibeta(q)};
            Tpow = Tgamma{igamma(q)};
            Xm = Xm + level.coeff(q) .* Tpow .* Eder;
        end

        level_contrib(ell, :) = prefactor .* Xm;
        xn = xn + level_contrib(ell, :);

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
    info.base_prefactor = cache.base_prefactor;
    info.prefactor = prefactor;
    info.n = n;
    info.alpha_eps = alpha_eps;
    info.cache = cache;
    info.level_contrib = level_contrib;
end


function [cache, SeriesCache] = get_or_build_series_cache(lambda, Mmax, SeriesCache)
    key = lambda_key(lambda);
    if isKey(SeriesCache, key)
        cache = SeriesCache(key);
    else
        cache = build_series_cache(lambda, Mmax);
        SeriesCache(key) = cache;
    end
end


function [val, MLCache] = get_cached_ml_eval(z, alpha, beta, gamma, fact_m, MLCache)
    key = ml_vector_cache_key(z, alpha, beta, gamma);

    if isKey(MLCache, key)
        val = MLCache(key);
    else
        val = fact_m .* ml(z, alpha, beta, gamma);
        MLCache(key) = val;
    end
end


function key = lambda_key(lambda)
    key = sprintf('%.16e,', lambda);
end


function key = ml_vector_cache_key(z, alpha, beta, gamma)
    if isscalar(z)
        key = sprintf('S|%.16e|%.16e|%.16e|%.16e', z, alpha, beta, gamma);
    else
        key = sprintf('V|%.16e|%.16e|%.16e|%d|', alpha, beta, gamma, numel(z));
        key = [key, sprintf('%.16e,', z)];
    end
end


function cache = build_series_cache(lambda, Mmax)
    lambda = lambda(:).';
    n = numel(lambda);

    a = bateman_poly_coeffs(lambda);
    a0 = a(1);
    a1 = a(2);
    mu = a1 / a0;
    c  = a(end:-1:3) / a0;

    % Independent of x10 and therefore safe to cache.
    base_prefactor = prod(lambda(1:n-1)) / a0;

    w_beta  = (n-1):-1:1;
    w_gamma = n:-1:2;
    levels = cell(Mmax + 1, 1);
    p = n - 1;

    for m = 0:Mmax
        % Generate K and log_kfact_vec on the fly
        K = generate_multi_indices(m, p);
        
        % Calculate log_kfact_vec using the properties of the Gamma function
        log_kfact_vec = sum(gammaln(K + 1), 2);
        Nr = size(K, 1);

        beta_vec  = zeros(Nr, 1);
        gamma_vec = zeros(Nr, 1);
        coeff_vec = zeros(Nr, 1);
        sgn = (-1)^m;

        for r = 1:Nr
            kvec = K(r, :);
            beta_vec(r)  = n + sum(w_beta  .* kvec);
            gamma_vec(r) = (n - 1) + sum(w_gamma .* kvec);

            coeff_prod = 1;
            zero_hit = false;
            for q = 1:numel(c)
                if c(q) == 0
                    if kvec(q) > 0
                        zero_hit = true;
                        break;
                    end
                elseif kvec(q) ~= 0
                    coeff_prod = coeff_prod * c(q)^kvec(q);
                end
            end

            if zero_hit
                coeff_vec(r) = 0;
            else
                coeff_vec(r) = sgn * coeff_prod * exp(-log_kfact_vec(r));
            end
        end

        % Remove exact zeros before aggregation.
        nz = (coeff_vec ~= 0);
        beta_vec  = beta_vec(nz);
        gamma_vec = gamma_vec(nz);
        coeff_vec = coeff_vec(nz);

        if isempty(coeff_vec)
            level = struct();
            level.m = m;
            level.beta = zeros(0,1);
            level.gamma = zeros(0,1);
            level.coeff = zeros(0,1);
            levels{m + 1} = level;
            continue;
        end

        % Aggregate repeated (beta,gamma) pairs.
        BG = [beta_vec, gamma_vec];
        [BGuniq, ~, ic] = unique(BG, 'rows');
        coeff_aggr = accumarray(ic, coeff_vec, [], @sum);

        % Remove tiny coefficients generated by cancellation.
        keep = abs(coeff_aggr) > 0;
        BGuniq = BGuniq(keep, :);
        coeff_aggr = coeff_aggr(keep);

        level = struct();
        level.m = m;
        level.beta = BGuniq(:, 1);
        level.gamma = BGuniq(:, 2);
        level.coeff = coeff_aggr;
        levels{m + 1} = level;
    end

    cache = struct();
    cache.a = a;
    cache.mu = mu;
    cache.c = c;
    cache.base_prefactor = base_prefactor;
    cache.levels = levels;
end

function K_matrix = generate_multi_indices(m, p)
    % Generates all weak compositions of 'm' into 'p' parts.
    % Each row of K_matrix is a vector [k_0, k_1, ..., k_{p-1}] such that sum(k) = m.
    
    if p == 1
        K_matrix = m;
        return;
    end
    
    if m == 0
        K_matrix = zeros(1, p);
        return;
    end
    
    % Use the "stars and bars" combinatorial method
    combinations = nchoosek(1:(m + p - 1), p - 1);
    num_rows = size(combinations, 1);
    K_matrix = zeros(num_rows, p);
    
    % Distribute values based on the positions of the "bars"
    K_matrix(:, 1) = combinations(:, 1) - 1;
    for i = 2:(p - 1)
        K_matrix(:, i) = combinations(:, i) - combinations(:, i - 1) - 1;
    end
    
    % FIX: The combinations matrix only has p-1 columns.
    K_matrix(:, p) = (m + p - 1) - combinations(:, p - 1);
end
