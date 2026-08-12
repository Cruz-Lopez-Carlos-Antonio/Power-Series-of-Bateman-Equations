% ================================================================
% Superposition-based time stepping for the integer-order
% Bateman chain using the optimized recursive power-series solver
% ================================================================

clear; clc;

%% ---------------------------------------------------------------
% INPUT DATA
%% ---------------------------------------------------------------
half_lives = [2,2,3,3,3,4];
lambda = log(2) ./ half_lives;
Time_vector = [0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.010];
%Time_vector = [1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100];
x10 = 6.023e23;

Mmax = 5;
tol  = 1e-10;

% Maximum internal substep size for the superposition stepping
dt_step_max = 1.0;

%% ---------------------------------------------------------------
% INITIAL CONDITION FOR THE FULL CHAIN
%% ---------------------------------------------------------------
n = numel(lambda);
x_init = zeros(n,1);
x_init(1) = x10;

%% ---------------------------------------------------------------
% MAIN CALL
%% ---------------------------------------------------------------
tic
[Xhist, Xn_series, info] = BatemanSuperpositionStepperOptimized( ...
    lambda, Time_vector, x_init, dt_step_max, Mmax, tol);
elapsed_time = toc;

%% ---------------------------------------------------------------
% OPTIONAL DISPLAY (comment out for cleaner benchmarks)
%% ---------------------------------------------------------------
% T = table(Time_vector(:), Xn_series(:), ...
%     'VariableNames', {'t','Xn_series'});
% disp(T);

fprintf('\nChain length n          = %d\n', n);
fprintf('Maximum local terms     = %d\n', info.max_terms_used);
fprintf('Total kernel calls      = %d\n', info.total_kernel_calls);
fprintf('Total internal substeps = %d\n', info.total_substeps);
fprintf('Elapsed solver time     = %.6f s\n', elapsed_time);

%% ---------------------------------------------------------------
% EXPORT RESULTS
%% ---------------------------------------------------------------
filename = 'Bateman_superposition_recursive_results_optimized_2.txt';
fid = fopen(filename,'w');

fprintf(fid,'t\tXn_series\n');
for i = 1:length(Time_vector)
    fprintf(fid,'%d\t%.16e\n', Time_vector(i), Xn_series(i));
end

fclose(fid);

fprintf('Results written to %s\n', filename);

% ================================================================
% LOCAL FUNCTIONS
% ================================================================

function [Xhist, Xn_series, info] = BatemanSuperpositionStepperOptimized( ...
    lambda, Time_vector, x_init, dt_step_max, Mmax, tol)
% ------------------------------------------------
% Advances the FULL Bateman chain by superposition over substeps.
% Uses an optimized local terminal-series solver for each subchain.
%
% INPUTS
%   lambda       : decay constants [lambda_1,...,lambda_n]
%   Time_vector  : output times
%   x_init       : full initial vector [x1(0),...,xn(0)]^T
%   dt_step_max  : maximum allowed internal substep
%   Mmax         : minimum series terms before convergence test
%   tol          : relative tolerance for the local series
%
% OUTPUTS
%   Xhist        : matrix (#times x n), full chain history
%   Xn_series    : last nuclide history
%   info         : diagnostics
% ------------------------------------------------

    lambda = lambda(:).';
    Time_vector = Time_vector(:);
    x = x_init(:);

    n  = numel(lambda);
    nt = numel(Time_vector);

    if numel(x) ~= n
        error('x_init must have the same length as lambda.');
    end

    Xhist = zeros(nt, n);

    t_current = 0.0;

    max_terms_used    = 0;
    total_kernel_calls = 0;
    total_substeps    = 0;

    for it = 1:nt

        t_target = Time_vector(it);
        Delta = t_target - t_current;

        if Delta < 0
            error('Time_vector must be nondecreasing.');
        end

        if Delta == 0
            Xhist(it,:) = x.';
            continue;
        end

        % Number of internal substeps
        nsub = max(1, ceil(Delta / dt_step_max));
        h = Delta / nsub;

        for s = 1:nsub
            [x, step_info] = AdvanceOneStepSuperpositionOptimized(lambda, x, h, Mmax, tol);

            max_terms_used     = max(max_terms_used, step_info.max_terms_used);
            total_kernel_calls = total_kernel_calls + step_info.kernel_calls;
            total_substeps     = total_substeps + 1;
        end

        Xhist(it,:) = x.';
        t_current = t_target;
    end

    Xn_series = Xhist(:, end);

    info.max_terms_used     = max_terms_used;
    info.total_kernel_calls = total_kernel_calls;
    info.total_substeps     = total_substeps;
end

% ----------------------------------------------------------------
function [x_next, info] = AdvanceOneStepSuperpositionOptimized(lambda, x_current, dt, Mmax, tol)
% ------------------------------------------------
% Performs one time step using superposition.
%
% For each j = 1,...,n:
%   - take the subchain lambda(j:end)
%   - assume only x_j is initially nonzero
%   - compute contributions to x_j, x_{j+1}, ..., x_n
%   - accumulate them
% ------------------------------------------------

    n = numel(lambda);
    x_next = zeros(n,1);

    max_terms_used = 0;
    kernel_calls = 0;

    for j = 1:n

        xj0 = x_current(j);

        if xj0 == 0
            continue;
        end

        lambda_sub_full = lambda(j:end);
        Lfull = numel(lambda_sub_full);

        % Contribution to every downstream nuclide
        for ell = 1:Lfull
            lambda_sub = lambda_sub_full(1:ell);

            [contrib, terms_used] = TerminalSeriesSubchainOptimized( ...
                lambda_sub, dt, xj0, Mmax, tol);

            i_global = j + ell - 1;
            x_next(i_global) = x_next(i_global) + contrib;

            if terms_used > max_terms_used
                max_terms_used = terms_used;
            end
            kernel_calls = kernel_calls + 1;
        end
    end

    info.max_terms_used = max_terms_used;
    info.kernel_calls   = kernel_calls;
end

% ----------------------------------------------------------------
function [xn_val, terms_used] = TerminalSeriesSubchainOptimized(lambda_sub, t, x10, Mmax, tol)
% ------------------------------------------------
% Computes the LAST nuclide of a subchain
%
%   x_1 -> x_2 -> ... -> x_L
%
% with initial condition
%
%   x_1(0) = x10,  x_2(0)=...=x_L(0)=0
%
% using the optimized recursive series.
% ------------------------------------------------

    lambda_sub = lambda_sub(:).';
    L = numel(lambda_sub);

    % Trivial case: one nuclide only
    if L == 1
        xn_val = x10 * exp(-lambda_sub(1) * t);
        terms_used = 1;
        return;
    end

    % Prefactor
    prefactor = x10 * prod(lambda_sub(1:L-1));

    % Ratios lambda_r / lambda_{r+1}
    ratios = lambda_sub(1:L-1) ./ lambda_sub(2:L);

    % ratio_powers(r) = (lambda_r/lambda_{r+1})^m
    ratio_powers = ones(1, L-1);   % for m = 0

    % CurrentS(r+1) = S_r^(m), r = 0,...,L-1
    % For m = 0: S_r^(0)=1
    CurrentS = ones(1, L);
    NewS     = zeros(1, L);

    % m = 0 term
    lambda_L_power = 1.0;                 % (-lambda_L)^0
    term_base = t^(L-1) / gamma(L);       % t^(L-1)/Gamma(L)

    sum_val = lambda_L_power * CurrentS(L) * term_base;
    terms_used = 1;

    consecutive_small = 0;
    small_needed = 3;

    max_iter = 100000;

    for m = 1:max_iter

        % Update ratio powers: from m-1 to m
        ratio_powers = ratio_powers .* ratios;

        % Update recurrence S_r^(m)
        NewS(1) = 1.0;
        for r = 1:(L-1)
            NewS(r+1) = CurrentS(r+1) + ratio_powers(r) * NewS(r);
        end
        CurrentS = NewS;

        % Update (-lambda_L)^m recursively
        lambda_L_power = lambda_L_power * (-lambda_sub(L));

        % Update t^(m+L-1)/Gamma(m+L) recursively
        term_base = term_base * t / (m + L - 1);

        % New term
        term = lambda_L_power * CurrentS(L) * term_base;
        sum_val = sum_val + term;
        terms_used = terms_used + 1;

        % Convergence test
        if m >= Mmax
            if abs(term) <= tol * max(1, abs(sum_val))
                consecutive_small = consecutive_small + 1;
            else
                consecutive_small = 0;
            end

            if consecutive_small >= small_needed
                break;
            end
        end
    end

    if m == max_iter
        warning('TerminalSeriesSubchainOptimized reached max_iter without convergence.');
    end

    xn_val = prefactor * sum_val;
end
