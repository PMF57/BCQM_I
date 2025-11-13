import numpy as np
from bcqm_tools.ramsey import compute_D, W_from_D

def test_constant_gamma_D():
    t = np.linspace(0, 10.0, 101)
    gamma = 0.3
    D = compute_D(t, gamma=gamma)
    assert np.allclose(D, np.exp(-gamma*t))

def test_W_bound_constant_gamma():
    gamma = 2.0
    fstar = 0.9
    t = np.linspace(0, 5.0, 5001)
    D = compute_D(t, gamma=gamma)
    W, idx = W_from_D(t, D, fstar)
    W_theory = (1.0/gamma)*np.log(1.0/(2.0*fstar - 1.0))
    assert W is not None
    assert abs(W - W_theory) <= (t[1]-t[0])
