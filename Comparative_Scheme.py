import math
import time
import numpy as np
from scipy.integrate import quad
import numpy.linalg as npl
from scipy.special import gamma, gammaincc  # Palma

# ------------------------------------------------------------
# Parámetros físicos (comunes a todos los métodos)
# ------------------------------------------------------------
a       = 1.0e-4     # pendiente de la rampa, rho(t) = a t + b
beta    = 0.0075     # fracción de retrasados
lam     = 0.001      # lambda (decaimiento de precursores)
Lambda1 = 0.0015     # tiempo de generación de neutrones
q       = 1.0e8      # fuente externa
b       = -6.0e-5    # término constante de la rampa

# Condición inicial típica de arranque
n0  = q * Lambda1 / abs(b)
dn0 = 0.0

# Notación interna para la solución analítica (Smets)
gamma_1  = a
lambda_1 = lam
Lambda_1 = Lambda1
rho_0    = b

# Cuadratura para tus integrales (Sección 6 de tu artículo)
_QKWARGS = dict(epsabs=0.0, epsrel=1e-9, limit=200)

# ============================================================
# 1) Solución analítica tipo Smets (tu Eq. (69))
# ============================================================

def z(t, rho_0, beta, Lambda_1, gamma_1, lambda_1):
    a1 = math.sqrt(gamma_1 / Lambda_1)
    return a1 * (t + (rho_0 - beta)/gamma_1 + Lambda_1*lambda_1/gamma_1)

def I_1(mu, zval):
    def integrand(x):
        return x**mu * np.exp(-x**2/2 + x*zval)
    val, _ = quad(integrand, 0.0, np.inf, **_QKWARGS)
    return float(val)

def I_2(mu, zval):
    def integrand(x):
        return x**mu * np.exp(-x**2/2 - x*zval)
    val, _ = quad(integrand, 0.0, np.inf, **_QKWARGS)
    return float(val)

def I_3(mu, zval):
    def integrand(x):
        return x**(mu+1) * np.exp(-x**2/2 + x*zval)
    val, _ = quad(integrand, 0.0, np.inf, **_QKWARGS)
    return float(val)

def I_4(mu, zval):
    def integrand(x):
        return x**(mu+1) * np.exp(-x**2/2 - x*zval)
    val, _ = quad(integrand, 0.0, np.inf, **_QKWARGS)
    return float(val)

def I_5(t, rho_0, lambda_1, beta, gamma_1, Lambda_1):
    mu = lambda_1*beta/gamma_1
    A = Lambda_1/(2*gamma_1)
    B = (rho_0 - beta)/gamma_1 + t
    def integrand(x):
        return (x + lambda_1)**mu * np.exp(-(A*x**2 - B*x))
    val, _ = quad(integrand, 0.0, np.inf, **_QKWARGS)
    return float(val)

def I_6(t, rho_0, lambda_1, beta, gamma_1, Lambda_1):
    mu = lambda_1*beta/gamma_1
    A = Lambda_1/(2*gamma_1)
    B = (rho_0 - beta)/gamma_1 + t
    def integrand(x):
        return x*(x + lambda_1)**mu * np.exp(-(A*x**2 - B*x))
    val, _ = quad(integrand, 0.0, np.inf, **_QKWARGS)
    return float(val)

def prefactor_F(lambda_1, beta, gamma_1, Lambda_1):
    mu = lambda_1*beta/gamma_1
    return (Lambda_1/gamma_1) * (lambda_1**(-mu))

def Constants_Ini_con(lambda_1, beta, Lambda_1, gamma_1, rho_0, n_0, dn_0, q):
    zeta = z(0.0, rho_0, beta, Lambda_1, gamma_1, lambda_1)
    mu   = lambda_1*beta/gamma_1

    Int1 = I_1(mu, zeta)
    Int2 = I_2(mu, zeta)
    Int3 = I_3(mu, zeta)
    Int4 = I_4(mu, zeta)
    Int5 = I_5(0.0, rho_0, lambda_1, beta, gamma_1, Lambda_1)
    Int6 = I_6(0.0, rho_0, lambda_1, beta, gamma_1, Lambda_1)

    F = prefactor_F(lambda_1, beta, gamma_1, Lambda_1)

    l1 = np.array([
        [Int1,                                  Int2],
        [math.sqrt(gamma_1/Lambda_1)*Int3, -math.sqrt(gamma_1/Lambda_1)*Int4]
    ], dtype=float)

    rhs1 = n_0 - q * F * Int5
    rhs2 = dn_0 + lambda_1 * (n_0 - q * F * Int5) - q * F * Int6

    l2 = np.array([rhs1, rhs2], dtype=float)

    col_norms = npl.norm(l1, axis=0)
    l1_scaled = l1 / col_norms

    max_rhs = np.max(np.abs(l2))
    if max_rhs == 0.0:
        max_rhs = 1.0
    l2_scaled = l2 / max_rhs

    A_scaled, *_ = npl.lstsq(l1_scaled, l2_scaled, rcond=None)
    A1, A2 = (A_scaled / col_norms) * max_rhs
    return float(A1), float(A2)

def Analytic_n(t, rho_0, beta, Lambda_1, gamma_1, lambda_1, n_0, dn_0, q):
    A1, A2 = Constants_Ini_con(lambda_1, beta, Lambda_1, gamma_1, rho_0, n_0, dn_0, q)
    mu   = lambda_1*beta/gamma_1
    zeta = z(t, rho_0, beta, Lambda_1, gamma_1, lambda_1)
    F    = prefactor_F(lambda_1, beta, gamma_1, Lambda_1)

    first  = A1 * math.exp(-lambda_1 * t) * I_1(mu, zeta)
    second = A2 * math.exp(-lambda_1 * t) * I_2(mu, zeta)
    third  = q  * F * I_5(t, rho_0, lambda_1, beta, gamma_1, Lambda_1)

    return float(first + second + third)

# ============================================================
# 2) Zhang
# ============================================================

def n_zhang(t, a, beta, Lambda1, q, b, n0):
    numerator   = beta * n0 + q * Lambda1
    denominator = beta - (a * t + b)
    return numerator / denominator

# ============================================================
# 3) Palma et al. (a partir de tu Corrected Palma.py)
# ============================================================

# Mapeo de notación:
q0  = b
r   = a
k   = lam
ell = Lambda1

k1 = (k * q0 + r) / r
k2 = (beta - q0) / r      # = A3
k3 = k * q * ell / r
s  = 1.0 + k * beta / r   # orden del exponente en Palma

A3 = k2

def upper_incomplete_gamma(s_param, z):
    # Γ(s, z) = Γ(s) * gammaincc(s, z)
    return gamma(s_param) * gammaincc(s_param, z)

A2 = -upper_incomplete_gamma(s, k * A3) \
     + (r * n0 * (k * A3)**s) / (k * q * ell * np.exp(k * A3))

B_palma = n0 * (A3 ** s) / (upper_incomplete_gamma(s, k * A3) + A2)

def n_palma(t):
    t = float(t)
    x = A3 - t
    val = B_palma * math.exp(-k * t) / (x**s) * (upper_incomplete_gamma(s, k * x) + A2)
    return float(val)

# ============================================================
# 4) Runge–Kutta 4° orden (1 grupo de retardados, rampa)
# ============================================================

def rho_t(t, a, b):
    return a * t + b

def rhs_point_kinetics(t, y, a, b, beta, lam, Lambda1, q):
    """
    y = [n, C]
    """
    n, C = y
    rho = rho_t(t, a, b)
    dn = (rho - beta)/Lambda1 * n + lam * C + q
    dC = beta/Lambda1 * n - lam * C
    return np.array([dn, dC], dtype=float)

def rk4_step(f, t, y, h, args):
    k1 = f(t,          y,              *args)
    k2 = f(t + 0.5*h,  y + 0.5*h*k1,   *args)
    k3 = f(t + 0.5*h,  y + 0.5*h*k2,   *args)
    k4 = f(t + h,      y + h*k3,       *args)
    return y + (h/6.0)*(k1 + 2*k2 + 2*k3 + k4)

def rk4_solve_and_sample(times, a, b, beta, lam, Lambda1, q, n0):
    """
    Integra una vez el sistema entre t=0 y t_max usando RK4
    con paso fijo h, y luego devuelve n(t) en los tiempos 'times'.
    """
    times = np.asarray(times, dtype=float)
    t_max = float(times.max())
    # Paso de integración: ajústalo para la precisión que quieras.
    h = 1.0e-3  # por ejemplo, 0.001 s -> 20,000 pasos hasta t=20
    N = int(round(t_max / h))

    # Condiciones iniciales
    C0 = beta * n0 / (lam * Lambda1)
    y = np.array([n0, C0], dtype=float)

    # Arrays para guardar la solución completa (solo para muestreo)
    n_full = np.empty(N+1, dtype=float)
    t_full = np.empty(N+1, dtype=float)

    n_full[0] = n0
    t_full[0] = 0.0

    args = (a, b, beta, lam, Lambda1, q)

    t = 0.0
    for k in range(1, N+1):
        y = rk4_step(rhs_point_kinetics, t, y, h, args)
        t = k * h
        n_full[k] = y[0]
        t_full[k] = t

    # Muestreamos n(t) en los tiempos requeridos (t son múltiplos de h)
    n_out = np.empty_like(times)
    inv_h = 1.0/h
    for i, tt in enumerate(times):
        idx = int(round(tt * inv_h))
        n_out[i] = n_full[idx]

    return n_out

# ============================================================
# Rutina genérica de benchmark con promedio (para métodos n(t))
# ============================================================

def benchmark_method(method, times, args=(), n_runs=5, warmup=True):
    """
    method(t, *args) escalar: se evalúa punto a punto.
    """
    times = np.asarray(times, dtype=float)
    if warmup:
        _ = [method(t, *args) for t in times]

    elapsed_list = []
    for _ in range(n_runs):
        t0 = time.perf_counter()
        _ = [method(t, *args) for t in times]
        t1 = time.perf_counter()
        elapsed_list.append(t1 - t0)

    elapsed_array = np.array(elapsed_list)
    mean_elapsed = elapsed_array.mean()
    std_elapsed  = elapsed_array.std(ddof=1) if n_runs > 1 else 0.0
    time_per_eval = mean_elapsed / len(times)
    return mean_elapsed, std_elapsed, time_per_eval

def benchmark_rk4(times, a, b, beta, lam, Lambda1, q, n0, n_runs=5, warmup=True):
    """
    Benchmark específico para RK4: integra una vez todo el intervalo
    y muestrea en 'times'.
    """
    times = np.asarray(times, dtype=float)
    if warmup:
        _ = rk4_solve_and_sample(times, a, b, beta, lam, Lambda1, q, n0)

    elapsed_list = []
    for _ in range(n_runs):
        t0 = time.perf_counter()
        _ = rk4_solve_and_sample(times, a, b, beta, lam, Lambda1, q, n0)
        t1 = time.perf_counter()
        elapsed_list.append(t1 - t0)

    elapsed_array = np.array(elapsed_list)
    mean_elapsed = elapsed_array.mean()
    std_elapsed  = elapsed_array.std(ddof=1) if n_runs > 1 else 0.0
    time_per_eval = mean_elapsed / len(times)
    return mean_elapsed, std_elapsed, time_per_eval

# ============================================================
# Programa principal: comparación en t = 0,1,...,20
# ============================================================

if __name__ == "__main__":

    times = np.arange(0, 21, 1, dtype=float)

    analytic_args = (rho_0, beta, Lambda_1, gamma_1, lambda_1, n0, dn0, q)
    zhang_args    = (a, beta, Lambda1, q, b, n0)

    n_runs = 5

    print("=== Benchmark en t = 0,1,2,...,20 ===")
    print(f"Parámetros: a={a}, beta={beta}, lambda={lam}, Lambda={Lambda1}, q={q}, b={b}")
    print(f"n0 inicial = {n0:.6e}")
    print()

    # --- Smets (tu método analítico) ---
    mean_A, std_A, per_eval_A = benchmark_method(
        Analytic_n,
        times,
        args=analytic_args,
        n_runs=n_runs,
        warmup=True
    )
    print("Solución analítica (Eq. (69), método de Smets):")
    print(f"  Tiempo medio total: {mean_A:.6f} s ± {std_A:.6f} s")
    print(f"  Tiempo medio por evaluación: {per_eval_A:.6e} s\n")

    # --- Zhang ---
    mean_Z, std_Z, per_eval_Z = benchmark_method(
        n_zhang,
        times,
        args=zhang_args,
        n_runs=n_runs,
        warmup=True
    )
    print("Solución aproximada de Zhang et al.:")
    print(f"  Tiempo medio total: {mean_Z:.6f} s ± {std_Z:.6f} s")
    print(f"  Tiempo medio por evaluación: {per_eval_Z:.6e} s\n")

    # --- Palma ---
    mean_P, std_P, per_eval_P = benchmark_method(
        n_palma,
        times,
        args=(),
        n_runs=n_runs,
        warmup=True
    )
    print("Solución analítica de Palma et al.:")
    print(f"  Tiempo medio total: {mean_P:.6f} s ± {std_P:.6f} s")
    print(f"  Tiempo medio por evaluación: {per_eval_P:.6e} s\n")

    # --- RK4 ---
    mean_R, std_R, per_eval_R = benchmark_rk4(
        times,
        a, b, beta, lam, Lambda1, q, n0,
        n_runs=n_runs,
        warmup=True
    )
    print("Solución numérica Runge–Kutta 4° orden:")
    print(f"  Tiempo medio total: {mean_R:.6f} s ± {std_R:.6f} s")
    print(f"  Tiempo medio por evaluación (sobre estos 21 puntos): {per_eval_R:.6e} s\n")

    # (Opcional) Mostrar algunos valores comparados
    print("Comparación de n(t) en t = 0, 5, 10, 15, 20:")
    nA_list = [Analytic_n(t, *analytic_args) for t in [0,5,10,15,20]]
    nZ_list = [n_zhang(t, *zhang_args)       for t in [0,5,10,15,20]]
    nP_list = [n_palma(t)                    for t in [0,5,10,15,20]]
    nR_full = rk4_solve_and_sample([0,5,10,15,20], a, b, beta, lam, Lambda1, q, n0)

    for t, nA, nZ, nP, nR in zip([0,5,10,15,20], nA_list, nZ_list, nP_list, nR_full):
        print(f"  t = {t:4.1f} s -> Smets = {nA:.6e}, Zhang = {nZ:.6e}, Palma = {nP:.6e}, RK4 = {nR:.6e}")