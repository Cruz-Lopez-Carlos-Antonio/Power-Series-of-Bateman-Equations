clear; clc;

% ========================================================================
% Generalized Bateman solution with selective multiprecision
% ========================================================================
%
% This script evaluates a generalized Bateman-type closed-form solution for
% a linear transmutation/decay chain when repeated decay constants are
% present. The implementation follows a grouped representation in which
% identical decay constants are collected into a single value lambda_i with
% multiplicity mu_i + 1. (see Cruz-López et al., 2024, Annals of Nuclear Energy)
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
% This version uses selective variable-precision arithmetic. The integer
% combinatorial structure is kept in standard MATLAB arithmetic, while the
% numerically delicate products, sums, and exponential evaluations are
% carried out with vpa arithmetic.
% ========================================================================

% ------------------------------------------------------------------------
% Input data
% ------------------------------------------------------------------------
% Half-lives defining the decay constants lambda = log(2)/T. Repeated
% half-lives lead to repeated decay constants and are therefore represented
% through multiplicities in the grouped formula.

Half_lifes = [2, 2, 3, 3, 3, 4, 4, 4, 4];
D = log(2) ./ Half_lifes;

% Number of significant decimal digits used by MATLAB symbolic arithmetic.
% Increase this value if stronger cancellation is expected.
mp_digits = 16;
digits(mp_digits);

% Two typical time grids are useful for testing:
%   1) Moderate/large times, where the computation is less affected by
%      cancellation.
%   2) Very small times, where cancellation and round-off errors become
%      visible in standard double precision.

%Time_vector = [1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100];
Time_vector = [0.001, 0.002, 0.003, 0.004, 0.005, ...
               0.006, 0.007, 0.008, 0.009, 0.010];

% Initial number of atoms/nuclei in the first member of the chain.
%proposed originally by Dreher (Annals of Nuclear Energy 53 (2013))
X0 = 6.023e23;

% Preallocate the output vector. Since the final result is produced by
% variable-precision arithmetic, the output is stored as a symbolic vector.
Xn_GPS = sym(zeros(size(Time_vector)));

% ------------------------------------------------------------------------
% Evaluation over the selected time grid
% ------------------------------------------------------------------------
tic
for k = 1:length(Time_vector)
    Xn_GPS(k) = GPS(X0, D, Time_vector(k), mp_digits);
end

fprintf('Computed Bateman values using GPS formula:\n')
disp(table(Time_vector(:), string(vpa(Xn_GPS(:), 25)), ...
    'VariableNames', {'t', 'Xn_GPS'}))

% Save the computed values in a plain text file, using high-precision
% scientific notation to facilitate comparisons with Python or other codes.

filename = 'Bateman_GPS_results_multiprecision.txt';
fid = fopen(filename, 'w');
fprintf(fid, 't\tXn_GPS\n');
for k = 1:length(Time_vector)
    fprintf(fid, '%.16e\t%s\n', Time_vector(k), char(vpa(Xn_GPS(k), mp_digits)));
end
fclose(fid);

fprintf('Results written to %s\n', filename);
toc

% ========================================================================
% Local functions
% ========================================================================

function solution = GPS(X0, DC, t, mp_digits)
    % GPS evaluates the generalized closed-form Bateman expression.
    %
    % Inputs:
    %   X0 : initial amount in the first nuclide/member of the chain.
    %   DC : vector of decay constants. Repeated entries are allowed.
    %   t  : time at which X_N(t) is evaluated.
    %
    % Output:
    %   solution : value of the final chain member at time t, computed in
    %              selective variable precision.

    % Convert the continuous numerical quantities to variable precision.
    % The grouped structure itself is still identified from the original
    % double-precision vector DC, because exact equality of repeated decay
    % constants is already encoded there.
    DC_mp = vpa(DC, mp_digits);
    t_mp  = vpa(t, mp_digits);
    X0_mp = vpa(X0, mp_digits);

    % First prefactor, X_1(0)/lambda_N.
    term1 = X0_mp / DC_mp(end);

    % Group equal decay constants while preserving their first-appearance
    % order. These grouped values play the role of lambda_i in the formula.
    Au1_double = unique(DC, 'stable');
    Au1 = vpa(Au1_double, mp_digits);

    % Multiplicity parameter mu_i. If a grouped decay constant appears q
    % times in the original chain, then mu_i = q - 1.
    Mu = zeros(size(Au1_double));
    for q = 1:length(Au1_double)
        Mu(q) = sum(DC == Au1_double(q)) - 1;
    end

    % Product prefactor prod_k lambda_k^(mu_k+1).
    term2 = vpa(1, mp_digits);
    for z = 1:length(Au1)
        term2 = term2 * Au1(z)^(Mu(z) + 1);
    end

    % Main summation over the grouped decay constants lambda_i.
    s1 = vpa(0, mp_digits);

    for i = 1:length(Au1)
        % Exponential factor exp(-lambda_i t).
        p1 = exp(-Au1(i) * t_mp);

        % Implement the condition j ~= i by removing the i-th entry from
        % the lists of grouped decay constants and multiplicities.
        %
        % This is the MATLAB counterpart of the mathematical products over
        % all indices except i used in the original Python implementation
        Au2  = Au1([1:i-1, i+1:end]);
        Mu_m = Mu([1:i-1, i+1:end]);

        % Denominator product over j ~= i. The code keeps the same order of
        % multiplication as the corresponding Python implementation, which
        % is relevant when comparing double-precision round-off effects.
        p2 = vpa(1, mp_digits);
        for j = 1:length(Au2)
            p2 = p2 * (vpa(1, mp_digits) / ((Au2(j) - Au1(i))^(Mu_m(j) + 1)));
        end

        % Polynomial contribution associated with repeated roots. For each
        % ell, the coefficient chi_{i,mu_i-ell} is computed by enumerating
        % all restricted integer partitions of mu_i - ell.
        s2 = vpa(0, mp_digits);
        for ell = 0:Mu(i)
            z1 = t_mp^ell / vpa(factorial(ell), mp_digits);

            % Rows of L contain the integer vectors h_k satisfying
            % sum_{k ~= i} h_k = mu_i - ell.
            L = partitions_restricted(Mu(i) - ell, length(Au2));

            % Compute chi_{i,mu_i-ell}.
            z2 = chi(i, Mu, Au1, L, mp_digits);

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

function c = chi(i, Mu, Lambd, L, mp_digits)
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
    %   c     : value of chi_{i,r} in selective variable precision.

    c = vpa(0, mp_digits);

    % Implement k ~= i by removing the i-th entry. Aux_L contains the
    % multiplicities mu_k for k ~= i, and Aux2 contains the corresponding
    % grouped decay constants lambda_k.
    Aux_L = Mu([1:i-1, i+1:end]);
    Aux2  = Lambd([1:i-1, i+1:end]);

    % Sum over all restricted integer partitions. Each row of L is one
    % admissible vector (h_k)_{k ~= i}.
    for row = 1:size(L, 1)
        u = L(row, :);
        a = vpa(1, mp_digits);

        % Product over k ~= i of
        % binomial(h_k + mu_k, mu_k) / (lambda_i - lambda_k)^h_k.
        for k = 1:length(Aux_L)
            b_f = vpa(nchoosek(Aux_L(k) + u(k), Aux_L(k)), mp_digits);
            dif = vpa(1, mp_digits) / ((Lambd(i) - Aux2(k))^u(k));
            a = a * b_f * dif;
        end

        c = c + a;
    end
end
