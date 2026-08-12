clear; clc;

% =========================
% Input data, GPS-style
% =========================

X0 = 6.023e23;

Half_lifes = [2, 2, 3, 3, 3, 4];

DC = log(2) ./ Half_lifes;

L  = DC;
Nu = DC(1:end-1);

N0 = zeros(1, length(DC));
N0(1) = X0;

Time_vector = [0.001, 0.002, 0.003, 0.004, 0.005, ...
               0.006, 0.007, 0.008, 0.009, 0.010];

% =========================
% Output similar to Python
% =========================

for idx = 1:length(Time_vector)
    t = Time_vector(idx);
    solution = Nn(L, Nu, N0, t);
    fprintf('%.12e\n', solution);
end


% ============================================================
% Local functions
% ============================================================

function p = myprod(values)
    p = 1.0;
    for i = 1:length(values)
        p = p * values(i);
    end
end


function val = gkl(L, beta_k, l)
    val = sum(L(l:end) == beta_k);
end


function indices = Setk(L, beta_k)
    indices = find(L == beta_k);
end


function val = Kkm(L, beta_k, m)
    indices = Setk(L, beta_k);
    val = indices(end - m);
end


function S = Skl(L, beta_k, l)
    n = length(L);
    s = l:n;
    S = setdiff(s, Setk(L, beta_k));
end


function C = combinations_with_replacement(values, r)
    if r == 0
        C = zeros(1, 0);
        return;
    end

    if isempty(values)
        C = [];
        return;
    end

    numVals = length(values);
    indexComb = nchoosek(1:(numVals + r - 1), r);
    C = zeros(size(indexComb));

    for row = 1:size(indexComb, 1)
        idx = indexComb(row, :) - (0:r-1);
        C(row, :) = values(idx);
    end
end


function C = gklmcomb(L, beta_k, l, m)
    r = gkl(L, beta_k, l) - 1 - m;
    S = Skl(L, beta_k, l);
    C = combinations_with_replacement(S, r);
end


function total = degeneratefactor(L, beta_k, l, m)
    if gkl(L, beta_k, l) == 1 + m
        total = 1.0;
        return;
    end

    S = Skl(L, beta_k, l);

    if isempty(S)
        total = 0.0;
        return;
    end

    Comb = gklmcomb(L, beta_k, l, m);

    total = 0.0;

    for row = 1:size(Comb, 1)
        factor = 1.0;

        for j = 1:size(Comb, 2)
            elem = Comb(row, j);
            factor = factor / (beta_k - L(elem));
        end

        total = total + factor;
    end
end


function total = Akm(L, Nu, N0, beta_k, m)
    total = 0.0;

    upper = Kkm(L, beta_k, m);

    for l = 1:upper
        factor = N0(l) * myprod(Nu(l:end));

        S = Skl(L, beta_k, l);

        if ~isempty(S)
            for j = 1:length(S)
                ind = S(j);
                factor = factor / (L(ind) - beta_k);
            end
        end

        factor = factor * degeneratefactor(L, beta_k, l, m);
        total = total + factor;
    end

    total = total / factorial(m);
end


function Nnt = Nn(L, Nu, N0, t)
    Betas = unique(L, 'stable');
    lambda_n = L(end);

    Nnt = N0(end) * exp(-lambda_n * t);

    for i = 1:length(Betas)
        beta_k = Betas(i);

        if beta_k ~= lambda_n
            Nnt = Nnt + Akm(L, Nu, N0, beta_k, 0) * ...
                (exp(-beta_k * t) - exp(-lambda_n * t));
        end
    end

    for i = 1:length(Betas)
        beta_k = Betas(i);
        mu_k = gkl(L, beta_k, 1) - 1;

        if mu_k > 0
            polynomial_part = 0.0;

            for m = 1:mu_k
                polynomial_part = polynomial_part + ...
                    (t^m) * Akm(L, Nu, N0, beta_k, m);
            end

            Nnt = Nnt + polynomial_part * exp(-beta_k * t);
        end
    end
end
