from prime_fock_nd import first_primes


def test_first_primes():
    assert first_primes(10) == [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
