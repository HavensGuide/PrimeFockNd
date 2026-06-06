from prime_fock_nd.primes import first_primes


def test_many_primes_without_growth_issue():
    primes = first_primes(1000)
    assert len(primes) == 1000
    assert primes[0] == 2
    assert primes[-1] == 7919
