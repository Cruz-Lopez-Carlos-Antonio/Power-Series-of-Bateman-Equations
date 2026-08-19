clear; clc;

% ========================================================================
% Generalized Bateman solution with selective multiprecision
% using effective decay constants and effective branching ratios
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

% The last branch closes the cycle Pu-239 -> U-235.
branch_cycle_closure = b_eff(end); %#ok<NASGU>

DAYS = [ ...
    2,4,6,8,10,15,20,25,30,35,40,45,50,100];


Time_vector = DAYS * 24 * 3600;

X0 = 6.89185e-4;

mp_digits = 200;
digits(mp_digits);

n = numel(lambda_eff);
X_all = sym(zeros(numel(Time_vector), n));

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

        X_all(tt, r) = GPS(X0, DC_prefix, branch_prefix, t, mp_digits);

    end
end

fprintf('Computed Bateman values using multiprecision GPS formula:\n')

T = table(DAYS(:), ...
    string(vpa(X_all(:,1), 25)), ...
    string(vpa(X_all(:,2), 25)), ...
    string(vpa(X_all(:,3), 25)), ...
    string(vpa(X_all(:,4), 25)), ...
    string(vpa(X_all(:,5), 25)), ...
    string(vpa(X_all(:,6), 25)), ...
    string(vpa(X_all(:,7), 25)), ...
    'VariableNames', [{'days'}, cellstr(chain)]);

disp(T)

filename = 'Bateman_GPS_effective_branches_multiprecision.txt';
fid = fopen(filename, 'w');

fprintf(fid, 'days');
for r = 1:n
    fprintf(fid, '\t%s', chain(r));
end
fprintf(fid, '\n');

for tt = 1:numel(DAYS)
    fprintf(fid, '%.16e', DAYS(tt));
    for r = 1:n
        fprintf(fid, '\t%s', char(vpa(X_all(tt, r), mp_digits)));
    end
    fprintf(fid, '\n');
end

fclose(fid);

fprintf('Results written to %s\n', filename);
toc

% ========================================================================
% Local functions
% ========================================================================

function solution = GPS(X0, DC, branch, t, mp_digits)

    DC = DC(:).';
    branch = branch(:).';

    if numel(branch) ~= numel(DC) - 1
        error('The branch vector must have length numel(DC)-1.');
    end

    DC_mp = vpa(DC, mp_digits);
    t_mp  = vpa(t, mp_digits);
    X0_mp = vpa(X0, mp_digits);

    % Effective branching prefactor.
    Bn = vpa(prod(branch), mp_digits);

    % First prefactor, X_1(0)/lambda_N.
    term1 = X0_mp / DC_mp(end);

    % Group equal effective decay constants.
    Au1_double = unique(DC, 'stable');
    Au1 = vpa(Au1_double, mp_digits);

    % Multiplicity parameter mu_i.
    Mu = zeros(size(Au1_double));
    for q = 1:length(Au1_double)
        Mu(q) = sum(DC == Au1_double(q)) - 1;
    end

    % Product prefactor prod_k lambda_k^(mu_k+1).
    term2 = vpa(1, mp_digits);
    for z = 1:length(Au1)
        term2 = term2 * Au1(z)^(Mu(z) + 1);
    end

    s1 = vpa(0, mp_digits);

    for i = 1:length(Au1)

        p1 = exp(-Au1(i) * t_mp);

        Au2  = Au1([1:i-1, i+1:end]);
        Mu_m = Mu([1:i-1, i+1:end]);

        p2 = vpa(1, mp_digits);
        for j = 1:length(Au2)
            p2 = p2 * (vpa(1, mp_digits) / ...
                ((Au2(j) - Au1(i))^(Mu_m(j) + 1)));
        end

        s2 = vpa(0, mp_digits);

        for ell = 0:Mu(i)

            z1 = t_mp^ell / vpa(factorial(ell), mp_digits);

            L = partitions_restricted(Mu(i) - ell, length(Au2));

            z2 = chi(i, Mu, Au1, L, mp_digits);

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

function c = chi(i, Mu, Lambd, L, mp_digits)

    c = vpa(0, mp_digits);

    Aux_L = Mu([1:i-1, i+1:end]);
    Aux2  = Lambd([1:i-1, i+1:end]);

    for row = 1:size(L, 1)

        u = L(row, :);
        a = vpa(1, mp_digits);

        for k = 1:length(Aux_L)
            b_f = vpa(nchoosek(Aux_L(k) + u(k), Aux_L(k)), mp_digits);
            dif = vpa(1, mp_digits) / ((Lambd(i) - Aux2(k))^u(k));
            a = a * b_f * dif;
        end

        c = c + a;
    end
end
