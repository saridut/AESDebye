import numpy as np
import pytest
import aesdebye


@pytest.fixture
def calculator():
    """Return a standard DebyeCalculator instance."""
    return aesdebye.DebyeCalculator()


def brute_force_debye_partial(coords_i, coords_j, el_i, el_j, q_vec, calc=None, cross_symmetric=False):
    """
    Compute Debye scattering contribution between position set I and set J using
    the exact brute force double sum.
    
    If cross_symmetric is True, multiplies by 2 to account for both I -> J and J -> I directions.
    """
    if calc is not None and el_i not in ("None", "", None):
        f_i = np.array(calc.calculateASFProfile(list(q_vec), el_i))
    else:
        f_i = np.ones(len(q_vec))

    if calc is not None and el_j not in ("None", "", None):
        f_j = np.array(calc.calculateASFProfile(list(q_vec), el_j))
    else:
        f_j = np.ones(len(q_vec))

    diff = coords_i[:, None, :] - coords_j[None, :, :]
    dists = np.linalg.norm(diff, axis=-1).flatten()

    intensity = np.zeros(len(q_vec))
    for r in dists:
        if r < 1e-9:
            sinc = np.ones_like(q_vec)
        else:
            qr = q_vec * r
            sinc = np.where(qr > 1e-7, np.sin(qr) / qr, 1.0)
        intensity += sinc

    if cross_symmetric:
        intensity *= 2.0

    return intensity * f_i * f_j
