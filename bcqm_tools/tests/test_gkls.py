import numpy as np
from bcqm_tools.gkls import simulate_gkls

def test_pure_dephasing_rx_decay():
    t = np.linspace(0, 10.0, 1001)
    gamma_phi = 0.2
    r = simulate_gkls(t, r0=np.array([1.0,0.0,0.0]), gamma_phi=gamma_phi, gamma_relax=0.0)
    rx = r[:,0]
    assert np.allclose(rx, np.exp(-gamma_phi*t), atol=5e-3)
