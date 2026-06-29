% ================================================================
% File: Superposition_Bateman_with_DioCache_and_MLCache.txt
%
% Complete MATLAB implementation of a stepwise superposition solver
% for the last nuclide concentration X_n(t) in a linear decay chain,
% using:
%   1) a precomputed Diophantine dictionary (DioCache),
%   2) an in-memory cache for Mittag-Leffler evaluations (MLCache).
%
% Save this file as:
%   Superposition_Bateman_with_DioCache_and_MLCache.m
%
% Make sure:
%   - Garrappa's ml(...) is on the MATLAB path
%   - dio_cache_p14_m5.mat already exists
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

load('dio_cache_p14_m5.mat','DioCache');

MLCache = containers.Map('KeyType','char','ValueType','double');

tic
[Xn_series, history, MLCache] = bateman_superposition_solver( ...
    lambda, x10, Time_vector, dt_block, Mmax, tol, alpha_eps, DioCache, MLCache);

disp('Computed X_n(t):')
disp(table(Time_vector(:), Xn_series(:), 'VariableNames', {'t','Xn_series'}))

filename = 'Bateman_superposition_results.txt';
fid = fopen(filename,'w');
fprintf(fid,'t\tXn_series\n');
for i = 1:length(Time_vector)
    fprintf(fid,'%d\t%.16e\n', Time_vector(i), Xn_series(i));
end
fclose(fid);

fprintf('Results written to %s\n', filename);
fprintf('ML cache size = %d entries\n', MLCache.Count);
toc


function [Xn_out, history, MLCache] = bateman_superposition_solver( ...
    lambda, x10, Time_vector, dt_block, Mmax, tol, alpha_eps, DioCache, MLCache)

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
        fprintf('Processing time block %d / %d\n', b, numel(block_edges)-1);

        t_left  = block_edges(b);
        t_right = block_edges(b+1);

        mask = (Time_vector > t_left) & (Time_vector <= t_right);
        t_req = Time_vector(mask);

        tau_req = t_req - t_left;
        tau_end = t_right - t_left;

        history.block_states{b,1} = state0;

        [state_req, state_end, MLCache] = propagate_block_by_superposition( ...
            lambda, state0, tau_req, tau_end, Mmax, tol, alpha_eps, DioCache, MLCache);

        if ~isempty(t_req)
            Xn_out(mask) = state_req(:, n).';
        end

        history.block_states{b,2} = state_end;
        state0 = state_end;
    end
end


function [state_req, state_end, MLCache] = propagate_block_by_superposition( ...
    lambda, state0, tau_req, tau_end, Mmax, tol, alpha_eps, DioCache, MLCache)

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
            [sub_req, MLCache] = local_subchain_solution( ...
                lambda_sub, xj0, tau_req, Mmax, tol, alpha_eps, DioCache, MLCache);
            state_req(:, j:end) = state_req(:, j:end) + sub_req;
        end

        [sub_end, MLCache] = local_subchain_solution( ...
            lambda_sub, xj0, tau_end, Mmax, tol, alpha_eps, DioCache, MLCache);
        state_end(j:end) = state_end(j:end) + sub_end;
    end
end


function [sub_state, MLCache] = local_subchain_solution( ...
    lambda_sub, xj0, tau, Mmax, tol, alpha_eps, DioCache, MLCache)

    lambda_sub = lambda_sub(:).';
    tau = tau(:).';

    nr = numel(lambda_sub);
    Nt = numel(tau);
    sub_state = zeros(Nt, nr);

    for r = 1:nr
        lambda_prefix = lambda_sub(1:r);
        [xr, ~, MLCache] = local_series_last_node( ...
            lambda_prefix, xj0, tau, Mmax, tol, alpha_eps, DioCache, MLCache);
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


function [xn, info, MLCache] = local_series_last_node( ...
    lambda, x10, t, Mmax, tol, alpha_eps, DioCache, MLCache)

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
        info.prefactor = x10;
        info.n = 1;
        info.alpha_eps = alpha_eps;
        info.cache = [];
        return;
    end

    cache = build_series_cache(lambda, x10, Mmax, DioCache);

    z = -cache.mu .* t;
    xn = zeros(size(t));
    level_contrib = zeros(numel(cache.levels), numel(t));

    for ell = 1:numel(cache.levels)
        level = cache.levels{ell};
        fprintf('   series level %d / %d\n', ell, numel(cache.levels));

        Xm = zeros(size(t));

        for q = 1:numel(level.coeff)
            coeff_q = level.coeff(q);
            beta_q  = level.beta(q);
            gamma_q = level.gamma(q);
            m_q     = level.m;

            Eder = zeros(size(t));
            beta_eff = beta_q + alpha_eps*m_q;
            gamma_eff = m_q + 1;

            for jj = 1:numel(t)
                key = ml_cache_key(z(jj), alpha_eps, beta_eff, gamma_eff);

                if isKey(MLCache, key)
                    Eder(jj) = MLCache(key);
                else
                    val = factorial(m_q) .* ml(z(jj), alpha_eps, beta_eff, gamma_eff);
                    MLCache(key) = val;
                    Eder(jj) = val;
                end
            end

            term = coeff_q .* (t .^ gamma_q) .* Eder;
            Xm = Xm + term;
        end

        level_contrib(ell, :) = cache.prefactor .* Xm;
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
    info.prefactor = cache.prefactor;
    info.n = n;
    info.alpha_eps = alpha_eps;
    info.cache = cache;
    info.level_contrib = level_contrib;
end


function key = ml_cache_key(z, alpha, beta, gamma)
    key = sprintf('%.16e|%.16e|%.16e|%.16e', z, alpha, beta, gamma);
end


function cache = build_series_cache(lambda, x10, Mmax, DioCache)

    lambda = lambda(:).';
    n = numel(lambda);

    a = bateman_poly_coeffs(lambda);
    a0 = a(1);
    a1 = a(2);

    mu = a1 / a0;
    c  = a(end:-1:3) / a0;
    prefactor = x10 * prod(lambda(1:n-1)) / a0;

    w_beta  = (n-1):-1:1;
    w_gamma = n:-1:2;

    levels = cell(Mmax + 1, 1);
    p = n - 1;

    for m = 0:Mmax
        entry = DioCache{p, m + 1};
        K = entry.K;
        log_kfact_vec = entry.log_kfact;

        Nr = size(K, 1);
        beta_vec  = zeros(Nr, 1);
        gamma_vec = zeros(Nr, 1);
        coeff_vec = zeros(Nr, 1);

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
                else
                    coeff_prod = coeff_prod * c(q)^kvec(q);
                end
            end

            if zero_hit
                coeff_vec(r) = 0;
            else
                coeff_vec(r) = (-1)^m * coeff_prod * exp(-log_kfact_vec(r));
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
