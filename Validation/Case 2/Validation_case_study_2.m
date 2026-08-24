% Requeriments
%   - Garrappa's ml(...) is on the MATLAB path
%   - dio_cache_p14_m5.mat
% ================================================================

clear; clc;

chain = ["U-235","U-236","U-237","Np-237","Np-238","Pu-238","Pu-239"];

lambda_eff = [
    1.82437e-8,...
    2.43924e-9,...
    1.20615e-6,...
    1.03533e-8,...
    3.84806e-6,...
    1.22412e-8,...
    6.03196e-8];

b_eff = [
    1.75949e-1,...
    9.65835e-1,...
    9.85383e-1,...
    9.86681e-1,...
    9.84802e-1,...
    9.15377e-1,...
    1.51134e-5];

% To this decay chain U-235 -> ... -> Pu-239
% we use only the consecutive nodes
branch_eff = b_eff(1:numel(chain)-1);

DAYS = [ ...
    2,4,6,8,10,15,20,25,30,35,40,45,50,100];

%DAYS = [ ...
 %   0.00000E+00 2.50000E+00 1.25000E+01 2.50000E+01];

Time_vector = DAYS * 24 * 3600;

x10       = 6.89185e-4;
dt_block  =  3600;    % 1 h per block
Mmax      = 5;
tol       = 1e-16;
alpha_eps = 0.999999999999;
verbose   = false;

load('dio_cache_p14_m5.mat','DioCache');

MLCache = containers.Map('KeyType','char','ValueType','any');
SeriesCache = containers.Map('KeyType','char','ValueType','any');
FactCache = factorial(0:Mmax);

tic
[X_all, history, MLCache, SeriesCache] = bateman_superposition_solver( ...
    lambda_eff, branch_eff, x10, Time_vector, dt_block, Mmax, tol, alpha_eps, ...
    DioCache, MLCache, SeriesCache, FactCache, verbose);

T = array2table([DAYS(:), X_all], ...
    'VariableNames', [{'days'}, cellstr(chain)]);

disp(T)

filename = 'Bateman_chain_with_branches_results_20.txt';
writetable(T, filename, 'Delimiter', '\t');

fprintf('Results written to %s\n', filename);
fprintf('ML cache size     = %d entries\n', MLCache.Count);
fprintf('Series cache size = %d entries\n', SeriesCache.Count);
toc


function [X_out, history, MLCache, SeriesCache] = bateman_superposition_solver( ...
    lambda, branch, x10, Time_vector, dt_block, Mmax, tol, alpha_eps, ...
    DioCache, MLCache, SeriesCache, FactCache, verbose)

    lambda = lambda(:).';
    branch = branch(:).';
    Time_vector = Time_vector(:).';
    n = numel(lambda);

    if numel(branch) ~= n-1
        error('branch must have length n-1.');
    end

    tmax = max(Time_vector);
    block_edges = 0:dt_block:tmax;
    if block_edges(end) < tmax
        block_edges = [block_edges, tmax];
    end

    state0 = zeros(1, n);
    state0(1) = x10;
    X_out = zeros(numel(Time_vector), n);

    % t = 0
    zero_mask = (Time_vector == 0);
    X_out(zero_mask, :) = repmat(state0, sum(zero_mask), 1);

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
            lambda, branch, state0, tau_req, tau_end, Mmax, tol, alpha_eps, ...
            DioCache, MLCache, SeriesCache, FactCache, verbose);

        if ~isempty(t_req)
            X_out(mask, :) = state_req;
        end

        history.block_states{b,2} = state_end;
        state0 = state_end;
    end
end


function [state_req, state_end, MLCache, SeriesCache] = propagate_block_by_superposition( ...
    lambda, branch, state0, tau_req, tau_end, Mmax, tol, alpha_eps, ...
    DioCache, MLCache, SeriesCache, FactCache, verbose)

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
        branch_sub = branch(j:end);

        if ~isempty(tau_req)
            [sub_req, MLCache, SeriesCache] = local_subchain_solution( ...
                lambda_sub, branch_sub, xj0, tau_req, Mmax, tol, alpha_eps, ...
                DioCache, MLCache, SeriesCache, FactCache, verbose);

            state_req(:, j:end) = state_req(:, j:end) + sub_req;
        end

        [sub_end, MLCache, SeriesCache] = local_subchain_solution( ...
            lambda_sub, branch_sub, xj0, tau_end, Mmax, tol, alpha_eps, ...
            DioCache, MLCache, SeriesCache, FactCache, verbose);

        state_end(j:end) = state_end(j:end) + sub_end;
    end
end


function [sub_state, MLCache, SeriesCache] = local_subchain_solution( ...
    lambda_sub, branch_sub, xj0, tau, Mmax, tol, alpha_eps, ...
    DioCache, MLCache, SeriesCache, FactCache, verbose)

    lambda_sub = lambda_sub(:).';
    branch_sub = branch_sub(:).';
    tau = tau(:).';

    nr = numel(lambda_sub);
    Nt = numel(tau);
    sub_state = zeros(Nt, nr);

    for r = 1:nr
        lambda_prefix = lambda_sub(1:r);

        if r == 1
            branch_prefix = [];
        else
            branch_prefix = branch_sub(1:r-1);
        end

        [xr, ~, MLCache, SeriesCache] = local_series_last_node( ...
            lambda_prefix, branch_prefix, xj0, tau, Mmax, tol, alpha_eps, ...
            DioCache, MLCache, SeriesCache, FactCache, verbose);

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
    lambda, branch, x10, t, Mmax, tol, alpha_eps, ...
    DioCache, MLCache, SeriesCache, FactCache, verbose)

    lambda = lambda(:).';
    branch = branch(:).';
    t = t(:).';
    n = numel(lambda);

    if n == 1
        xn = x10 * exp(-lambda(1) * t);
        info = struct();
        return;
    end

    [cache, SeriesCache] = get_or_build_series_cache(lambda, branch, Mmax, DioCache, SeriesCache);

    z = -cache.mu .* t;
    xn = zeros(size(t));
    level_contrib = zeros(numel(cache.levels), numel(t));
    prefactor = x10 * cache.base_prefactor;

    for ell = 1:numel(cache.levels)

        level = cache.levels{ell};

        if isempty(level.coeff)
            continue;
        end

        m_q = level.m;
        fact_m = FactCache(m_q + 1);
        beta_eff_vec = level.beta + alpha_eps * m_q;

        [ubeta, ~, ibeta] = unique(beta_eff_vec);
        Ebeta = cell(numel(ubeta), 1);

        for ib = 1:numel(ubeta)
            [Ebeta{ib}, MLCache] = get_cached_ml_eval( ...
                z, alpha_eps, ubeta(ib), m_q + 1, fact_m, MLCache);
        end

        [ugamma, ~, igamma] = unique(level.gamma);
        Tgamma = cell(numel(ugamma), 1);

        for ig = 1:numel(ugamma)
            Tgamma{ig} = t .^ ugamma(ig);
        end

        Xm = zeros(size(t));

        for q = 1:numel(level.coeff)
            Xm = Xm + level.coeff(q) .* Tgamma{igamma(q)} .* Ebeta{ibeta(q)};
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


function [cache, SeriesCache] = get_or_build_series_cache(lambda, branch, Mmax, DioCache, SeriesCache)

    key = lambda_branch_key(lambda, branch);

    if isKey(SeriesCache, key)
        cache = SeriesCache(key);
    else
        cache = build_series_cache(lambda, branch, Mmax, DioCache);
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


function key = lambda_branch_key(lambda, branch)
    key = ['L|', sprintf('%.16e,', lambda), '|B|', sprintf('%.16e,', branch)];
end


function key = ml_vector_cache_key(z, alpha, beta, gamma)

    if isscalar(z)
        key = sprintf('S|%.16e|%.16e|%.16e|%.16e', z, alpha, beta, gamma);
    else
        key = sprintf('V|%.16e|%.16e|%.16e|%d|', alpha, beta, gamma, numel(z));
        key = [key, sprintf('%.16e,', z)];
    end
end


function cache = build_series_cache(lambda, branch, Mmax, DioCache)

    lambda = lambda(:).';
    branch = branch(:).';
    n = numel(lambda);

    a = bateman_poly_coeffs(lambda);
    a0 = a(1);
    a1 = a(2);

    mu = a1 / a0;
    c  = a(end:-1:3) / a0;

    base_prefactor = prod(branch(1:n-1) .* lambda(1:n-1)) / a0;

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

        BG = [beta_vec, gamma_vec];
        [BGuniq, ~, ic] = unique(BG, 'rows');
        coeff_aggr = accumarray(ic, coeff_vec, [], @sum);

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
