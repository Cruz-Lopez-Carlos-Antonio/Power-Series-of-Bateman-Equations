clear; clc;

% ========================================================================
% Generalized Bateman solution in standard double precision
% ========================================================================
%
% This script evaluates a generalized Bateman-type closed-form solution for
% a linear transmutation/decay chain when repeated decay constants are
% present. The implementation follows a grouped representation in which
% identical decay constants are collected into a single value lambda_i with
% multiplicity mu_i + 1.
%
% In the notation used in the associated manuscript, the implemented
% structure is essentially
%
%   X_N(t) = X_1(0)/lambda_N * prod_k lambda_k^(mu_k+1)
%            * sum_i exp(-lambda_i t)
%              / prod_{j ~= i} (lambda_j - lambda_i)^(mu_j+1)
%              * sum_{ell=0}^{mu_i} t^ell/ell! * chi_{i,mu_i-ell},
%
% where
%
%   chi_{i,r} = sum_{h_1+...+h_{i-1}+h_{i+1}+...+h_n = r}
%               prod_{k ~= i} binomial(h_k + mu_k, mu_k)
%               / (lambda_i - lambda_k)^h_k.
%
% The restriction j ~= i, or k ~= i, is implemented by removing the i-th
% entry from the corresponding MATLAB arrays. This is done, for example, by
% using the index vector [1:i-1, i+1:end].
%
% Important numerical note:
% This version intentionally uses standard IEEE double precision. For small
% times, the formula may involve severe cancellation among large terms, so
% the output can display the expected loss of numerical accuracy.
% ========================================================================

% ------------------------------------------------------------------------
% Input data
% ------------------------------------------------------------------------
% Half-lives defining the decay constants lambda = log(2)/T. Repeated
% half-lives lead to repeated decay constants and are therefore represented
% through multiplicities in the grouped formula.
Half_lifes = [2, 2, 3, 3, 3, 4];
D = log(2) ./ Half_lifes;

% Two typical time grids are useful for testing:
%   1) Moderate/large times, where the computation is less affected by
%      cancellation.
%   2) Very small times, where cancellation and round-off errors become
%      visible in standard double precision.
%Time_vector = [1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100];
Time_vector = [0.001, 0.002, 0.003, 0.004, 0.005, ...
               0.006, 0.007, 0.008, 0.009, 0.010];

% Initial number of atoms/nuclei in the first member of the chain.
X0 = 6.023e23;

% Preallocate the output vector.
Xn_GPS = zeros(size(Time_vector));

% ------------------------------------------------------------------------
% Evaluation over the selected time grid
% ------------------------------------------------------------------------
tic
for k = 1:length(Time_vector)
    Xn_GPS(k) = GPS(X0, D, Time_vector(k));
end

fprintf('Computed Bateman values using GPS formula:\n')
disp(table(Time_vector(:), Xn_GPS(:), ...
    'VariableNames', {'t', 'Xn_GPS'}))

% Save the computed values in a plain text file, using high-precision
% scientific notation to facilitate comparisons with Python or other codes.
filename = 'Bateman_GPS_results_double.txt';
fid = fopen(filename, 'w');
fprintf(fid, 't\tXn_GPS\n');
for k = 1:length(Time_vector)
    fprintf(fid, '%.16e\t%.16e\n', Time_vector(k), Xn_GPS(k));
end
fclose(fid);

fprintf('Results written to %s\n', filename);
toc

% ========================================================================
% Local functions
% ========================================================================

function solution = GPS(X0, DC, t)
    % GPS evaluates the generalized closed-form Bateman expression.
    %
    % Inputs:
    %   X0 : initial amount in the first nuclide/member of the chain.
    %   DC : vector of decay constants. Repeated entries are allowed.
    %   t  : time at which X_N(t) is evaluated.
    %
    % Output:
    %   solution : value of the final chain member at time t, computed in
    %              standard double precision.

    % First prefactor, X_1(0)/lambda_N.
    term1 = X0 / DC(end);

    % Group equal decay constants while preserving their first-appearance
    % order. These grouped values play the role of lambda_i in the formula.
    Au1 = unique(DC, 'stable');

    % Multiplicity parameter mu_i. If a grouped decay constant appears q
    % times in the original chain, then mu_i = q - 1.
    Mu = zeros(size(Au1));
    for q = 1:length(Au1)
        Mu(q) = sum(DC == Au1(q)) - 1;
    end

    % Product prefactor prod_k lambda_k^(mu_k+1).
    term2 = 1;
    for z = 1:length(Au1)
        term2 = term2 * Au1(z)^(Mu(z) + 1);
    end

    % Main summation over the grouped decay constants lambda_i.
    s1 = 0;

    for i = 1:length(Au1)
        % Exponential factor exp(-lambda_i t).
        p1 = exp(-Au1(i) * t);

        % Implement the condition j ~= i by removing the i-th entry from
        % the lists of grouped decay constants and multiplicities.
        %
        % This is the MATLAB counterpart of the mathematical products over
        % all indices except i.
        Au2  = Au1([1:i-1, i+1:end]);
        Mu_m = Mu([1:i-1, i+1:end]);

        % Denominator product over j ~= i. The code keeps the same order of
        % multiplication as the corresponding Python implementation, which
        % is relevant when comparing double-precision round-off effects.
        p2 = 1;
        for j = 1:length(Au2)
            p2 = p2 * (1 / ((Au2(j) - Au1(i))^(Mu_m(j) + 1)));
        end

        % Polynomial contribution associated with repeated roots. For each
        % ell, the coefficient chi_{i,mu_i-ell} is computed by enumerating
        % all restricted integer partitions of mu_i - ell.
        s2 = 0;
        for ell = 0:Mu(i)
            z1 = t^ell / factorial(ell);

            % Rows of L contain the integer vectors h_k satisfying
            % sum_{k ~= i} h_k = mu_i - ell.
            L = partitions_restricted(Mu(i) - ell, length(Au2));

            % Compute chi_{i,mu_i-ell}.
            z2 = chi(i, Mu(i) - ell, Mu, Au1, L);

            s2 = s2 + z1 * z2;
        end

        s1 = s1 + p1 * p2 * s2;
    end

    % Final value of X_N(t).
    solution = term1 * term2 * s1;
end

function L = partitions_restricted(m, n)
    % partitions_restricted generates all n-tuples of nonnegative integers
    % whose entries sum to m.
    %
    % In the expression for chi_{i,r}, these tuples represent the indices
    % (h_1,...,h_{i-1},h_{i+1},...,h_n) subject to
    %
    %   h_1 + ... + h_{i-1} + h_{i+1} + ... + h_n = r.
    %
    % The case n = 0 is included for completeness: the only valid empty
    % tuple occurs when m = 0.

    if n == 0
        if m == 0
            L = zeros(1, 0);
        else
            L = zeros(0, 0);
        end
        return;
    end

    L = partitions_rec(m, n);
end

function L = partitions_rec(m, n)
    % Recursive helper for partitions_restricted.
    %
    % The first component is chosen as v = 0,...,m. The remaining n-1
    % components are then generated recursively so that the total sum is m.

    if n == 1
        L = m;
        return;
    end

    L = [];
    for v = 0:m
        Tail = partitions_rec(m - v, n - 1);
        if isempty(Tail)
            L = [L; v]; %#ok<AGROW>
        else
            L = [L; [v * ones(size(Tail, 1), 1), Tail]]; %#ok<AGROW>
        end
    end
end

function c = chi(i, ~, Mu, Lambd, L)
    % chi evaluates the coefficient chi_{i,r}.
    %
    % Inputs:
    %   i     : selected grouped decay constant lambda_i.
    %   ~     : integer order r in chi_{i,r}. This value is already encoded
    %           in the rows of L, so MATLAB does not need to use it directly.
    %           It is kept in the call for consistency with the formula.
    %   Mu    : vector of multiplicity parameters mu_k.
    %   Lambd : vector of grouped decay constants lambda_k.
    %   L     : matrix whose rows contain the admissible vectors h_k.
    %
    % Output:
    %   c     : value of chi_{i,r} in standard double precision.

    c = 0;

    % Implement k ~= i by removing the i-th entry. Aux_L contains the
    % multiplicities mu_k for k ~= i, and Aux2 contains the corresponding
    % grouped decay constants lambda_k.
    Aux_L = Mu([1:i-1, i+1:end]);
    Aux2  = Lambd([1:i-1, i+1:end]);

    % Sum over all restricted integer partitions. Each row of L is one
    % admissible vector (h_k)_{k ~= i}.
    for row = 1:size(L, 1)
        u = L(row, :);
        a = 1;

        % Product over k ~= i of
        % binomial(h_k + mu_k, mu_k) / (lambda_i - lambda_k)^h_k.
        for k = 1:length(Aux_L)
            b_f = nchoosek(Aux_L(k) + u(k), Aux_L(k));
            dif = 1 / ((Lambd(i) - Aux2(k))^u(k));
            a = a * b_f * dif;
        end

        c = c + a;
    end
end
