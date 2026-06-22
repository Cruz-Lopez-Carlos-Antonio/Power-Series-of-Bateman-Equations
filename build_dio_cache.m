function build_dio_cache(max_p, max_m, filename)
% BUILD_DIO_CACHE
% ---------------------------------------------------------------
% Builds and saves a dictionary of nonnegative Diophantine solutions
% for equations of the form
%
%   x1 + x2 + ... + xp = m,   xi >= 0
%
% for:
%   p = 1,2,...,max_p
%   m = 0,1,...,max_m
%
% The cache is saved in a .mat file and can later be loaded once
% and reused by the Bateman superposition solver.
%
% INPUTS
%   max_p    : maximum number of variables p
%   max_m    : maximum truncation level m
%   filename : output .mat filename
%
% SAVED VARIABLE
%   DioCache : cell array such that
%              DioCache{p, m+1} is a struct with fields:
%                 .p
%                 .m
%                 .K          -> all solutions, one per row
%                 .log_kfact  -> sum(gammaln(K+1),2)
%                 .num_rows   -> number of solutions
%
% EXAMPLE
%   build_dio_cache(14, 5, 'dio_cache_p14_m5.mat');
% ---------------------------------------------------------------

    if nargin < 3 || isempty(filename)
        filename = sprintf('dio_cache_p%d_m%d.mat', max_p, max_m);
    end

    if max_p < 1 || max_m < 0
        error('Require max_p >= 1 and max_m >= 0.');
    end

    DioCache = cell(max_p, max_m + 1);

    fprintf('Building Diophantine cache up to p=%d, m=%d\n', max_p, max_m);
    tic

    for p = 1:max_p
        fprintf('  Processing p = %d / %d\n', p, max_p);

        for m = 0:max_m
            K = compositions_nonnegative_fast(m, p);
            log_kfact = sum(gammaln(K + 1), 2);

            entry = struct();
            entry.p = p;
            entry.m = m;
            entry.K = K;
            entry.log_kfact = log_kfact;
            entry.num_rows = size(K, 1);

            DioCache{p, m + 1} = entry;
        end
    end

    save(filename, 'DioCache', 'max_p', 'max_m', '-v7.3');

    elapsed = toc;
    fprintf('Cache saved to %s\n', filename);
    fprintf('Elapsed time: %.6f seconds\n', elapsed);
end


function K = compositions_nonnegative_fast(m, p)
% COMPOSITIONS_NONNEGATIVE_FAST
% ---------------------------------------------------------------
% Returns all nonnegative integer p-tuples [x1,...,xp] satisfying
%
%   x1 + ... + xp = m
%
% One solution per row.
%
% Uses the stars-and-bars construction.
% ---------------------------------------------------------------

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