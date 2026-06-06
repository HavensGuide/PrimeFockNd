import numpy as np

from prime_fock_nd import ValuationArithmetic, pairwise_valuation_kernel


def test_valuation_arithmetic_identities():
    logs = np.log(np.array([2, 3, 5], dtype=float))
    a = np.array([2, 0, 1], dtype=float)  # 2^2 * 5
    b = np.array([1, 3, 0], dtype=float)  # 2 * 3^3

    assert np.allclose(ValuationArithmetic.multiply(a, b), np.array([3, 3, 1]))
    assert np.allclose(ValuationArithmetic.gcd(a, b), np.array([1, 0, 0]))
    assert np.allclose(ValuationArithmetic.lcm(a, b), np.array([2, 3, 1]))
    assert ValuationArithmetic.divides(np.array([1, 0, 0]), a)
    assert not ValuationArithmetic.divides(np.array([3, 0, 0]), a)

    log_gcd = ValuationArithmetic.log_gcd(a, b, logs)
    assert np.isclose(log_gcd, np.log(2.0))

    dist = ValuationArithmetic.arithmetic_distance(a, b, logs)
    expected = np.abs(a - b) @ logs
    assert np.isclose(dist, expected)


def test_pairwise_kernels():
    logs = np.log(np.array([2, 3, 5], dtype=float))
    V = np.array([
        [1, 0, 1],
        [1, 2, 0],
        [0, 0, 1],
    ], dtype=float)

    kg = pairwise_valuation_kernel(V, logs, kernel="log_gcd").matrix
    kl = pairwise_valuation_kernel(V, logs, kernel="log_lcm").matrix
    kd = pairwise_valuation_kernel(V, logs, kernel="arithmetic_distance").matrix
    kj = pairwise_valuation_kernel(V, logs, kernel="jaccard_gcd_lcm").matrix
    kc = pairwise_valuation_kernel(V, logs, kernel="coprime_indicator").matrix

    assert np.isclose(kg[0, 1], np.log(2.0))
    assert np.isclose(kg[0, 2], np.log(5.0))
    assert kl[0, 1] > kg[0, 1]
    assert kd[0, 1] > 0
    assert 0 <= kj[0, 1] <= 1
    assert kc[1, 2] == 1.0
