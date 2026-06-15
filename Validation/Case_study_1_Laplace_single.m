clear; clc;

% ========================================================================
% Generalized Bateman solution in standard double precision
% with effective decay constants and effective branching ratios
% ========================================================================

chain = ["U-235", ...
         "U-236", ...
         "U-237", ...
         "Np-237", ...
         "Np-238", ...
         "Pu-238", ...
         "Pu-239"];

lambda_eff = [ ...
    1.82437e-8, ...
    2.43924e-9, ...
    1.20615e-6, ...
    1.03533e-8, ...
    3.84806e-6, ...
    1.22412e-8, ...
    6.03196e-8];

b_eff = [ ...
    1.75949e-1, ...
    9.65835e-1, ...
    9.85383e-1, ...
    9.86681e-1, ...
    9.84802e-1, ...
    9.15377e-1, ...
    1.51134e-5];

% For the broken linear chain U-235 -> ... -> Pu-239,
% only the first six branches are used.
branch_eff = b_eff(1:numel(chain)-1);

% The last branch would close the cycle Pu-239 -> U-235.
branch_cycle_closure = b_eff(end); %#ok<NASGU>

DAYS = [ ...
    0.00000E+00 2.50000E+00 1.25000E+01 2.50000E+01 ...
    1.25000E+02 2.50000E+02 3.75000E+02 5.00000E+02 ...
    6.25000E+02 7.50000E+02 8.75000E+02 1.00000E+03];

Time_vector = DAYS * 24 * 3600;

X0 = 6.89185e-4;

n = numel(lambda_eff);
X_all = zeros(numel(Time_vector), n);

tic

for tt = 1:numel(Time_vector)
    t = Time_vector(tt);

    for r = 1:n
        DC_prefix = lambda_eff(1:r);

        if r == 1
            branch_prefix = [];
        else
            branch_prefix = branch_eff(1:r-1);
        end

        X_all(tt, r) = GPS(X0, DC_prefix, branch_prefix, t);
    end
end

T = array2table([DAYS(:), X_all], ...
    'VariableNames', [{'days'}, cellstr(chain)]);

disp(T)

filename = 'Bateman_GPS_effective_branches_double.txt';
fid = fopen(filename, 'w');

fprintf(fid, 'days');
for r = 1:n
    fprintf(fid, '\t%s', chain(r));
end
fprintf(fid, '\n');

for tt = 1:numel(DAYS)
    fprintf(fid, '%.16e', DAYS(tt));
    for r = 1:n
        fprintf(fid, '\t%.16e', X_all(tt, r));
    end
    fprintf(fid, '\n');
end

fclose(fid);

fprintf('Results written to %s\n', filename);
toc

% ========================================================================
% Local functions
% ========================================================================

function solution = GPS(X0, DC, branch, t)

    DC = DC(:).';
    branch = branch(:).';

    if numel(branch) ~= numel(DC) - 1
        error('The branch vector must have length numel(DC)-1.');
    end

    % Effective branching prefactor.
    Bn = prod(branch);

    % First prefactor, X_1(0)/lambda_N.
    term1 = X0 / DC(end);

    % Group equal effective decay constants.
    Au1 = unique(DC, 'stable');

    % Multiplicity parameter mu_i.
    Mu = zeros(size(Au1));
    for q = 1:length(Au1)
        Mu(q) = sum(DC == Au1(q)) - 1;
    end

    % Product prefactor prod_k lambda_k^(mu_k+1).
    term2 = 1;
    for z = 1:length(Au1)
        term2 = term2 * Au1(z)^(Mu(z) + 1);
    end

    s1 = 0;

    for i = 1:length(Au1)

        p1 = exp(-Au1(i) * t);

        Au2  = Au1([1:i-1, i+1:end]);
        Mu_m = Mu([1:i-1, i+1:end]);

        p2 = 1;
        for j = 1:length(Au2)
            p2 = p2 * (1 / ((Au2(j) - Au1(i))^(Mu_m(j) + 1)));
        end

        s2 = 0;

        for ell = 0:Mu(i)
            z1 = t^ell / factorial(ell);

            L = partitions_restricted(Mu(i) - ell, length(Au2));

            z2 = chi(i, Mu, Au1, L);

            s2 = s2 + z1 * z2;
        end

        s1 = s1 + p1 * p2 * s2;
    end

    solution = Bn * term1 * term2 * s1;
end

function L = partitions_restricted(m, n)

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

function c = chi(i, Mu, Lambd, L)

    c = 0;

    Aux_L = Mu([1:i-1, i+1:end]);
    Aux2  = Lambd([1:i-1, i+1:end]);

    for row = 1:size(L, 1)

        u = L(row, :);
        a = 1;

        for k = 1:length(Aux_L)
            b_f = nchoosek(Aux_L(k) + u(k), Aux_L(k));
            dif = 1 / ((Lambd(i) - Aux2(k))^u(k));
            a = a * b_f * dif;
        end

        c = c + a;
    end
end
