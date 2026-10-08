import numpy as np
import pytest
import aesdebye
from conftest import brute_force_debye_partial


class TestBruteForceValidation:
    """Validate Debye calculation against pure Python brute force implementation."""

    def test_single_element_against_brute_force(self, calculator):
        """Validate single-element Debye profile against brute force double sum."""
        coords = np.array([
            [0.0, 0.0, 0.0],
            [2.77, 0.0, 0.0],
            [0.0, 2.77, 0.0],
            [2.77, 2.77, 0.0],
        ])
        pos = aesdebye.Positions(
            chemicalSymbols=["Pt"] * len(coords), coordinates=coords
        ).filterByElement("Pt")

        q_start, q_end, steps = 0.5, 8.0, 50
        pdf = calculator.calculatePDF(pos)
        profile = calculator.calculateIntensity(pdf, q_start, q_end, steps, False, 0.4)

        q_vec = np.array(profile.q)
        bf_intensity = brute_force_debye_partial(
            coords, coords, "Pt", "Pt", q_vec, calc=calculator, cross_symmetric=False
        )

        aes_intensity = np.array(profile.intensity)
        # Should match within numerical binning error (< 1e-4 relative tolerance)
        np.testing.assert_allclose(aes_intensity, bf_intensity, rtol=1e-4, atol=1e-6)

    def test_cross_elements_against_brute_force(self, calculator):
        """
        Validate multi-element system (Pt and Au) against brute force,
        checking partials (Pt-Pt, Au-Au, Pt-Au) and the total profile.
        """
        coords_pt = np.array([
            [0.0, 0.0, 0.0],
            [2.77, 0.0, 0.0],
            [0.0, 2.77, 0.0],
        ])
        coords_au = np.array([
            [2.77, 2.77, 0.0],
            [1.38, 1.38, 2.0],
        ])

        all_coords = np.vstack([coords_pt, coords_au])
        all_symbols = ["Pt"] * len(coords_pt) + ["Au"] * len(coords_au)
        pos_all = aesdebye.Positions(chemicalSymbols=all_symbols, coordinates=all_coords)

        q_start, q_end, steps = 1.0, 10.0, 60
        results = calculator.calculateProfile(pos_all, q_start, q_end, steps, False, 0.4, "")

        q_vec = np.array(results["total"][1].q)

        # 1. Brute force partials
        bf_pt_pt = brute_force_debye_partial(
            coords_pt, coords_pt, "Pt", "Pt", q_vec, calc=calculator, cross_symmetric=False
        )
        bf_au_au = brute_force_debye_partial(
            coords_au, coords_au, "Au", "Au", q_vec, calc=calculator, cross_symmetric=False
        )
        # Cross partial must be doubled (accounts for both Pt -> Au and Au -> Pt)
        bf_pt_au = brute_force_debye_partial(
            coords_pt, coords_au, "Pt", "Au", q_vec, calc=calculator, cross_symmetric=True
        )
        bf_total = bf_pt_pt + bf_au_au + bf_pt_au

        # 2. aesdebye partials & total
        aes_pt_pt = np.array(results["Pt-Pt"][1].intensity)
        aes_au_au = np.array(results["Au-Au"][1].intensity)
        aes_pt_au = np.array(results["Pt-Au"][1].intensity)
        aes_total = np.array(results["total"][1].intensity)

        # 3. Assertions
        np.testing.assert_allclose(aes_pt_pt, bf_pt_pt, rtol=1e-4, atol=1e-6)
        np.testing.assert_allclose(aes_au_au, bf_au_au, rtol=1e-4, atol=1e-6)
        np.testing.assert_allclose(aes_pt_au, bf_pt_au, rtol=1e-4, atol=1e-6)
        np.testing.assert_allclose(aes_total, bf_total, rtol=1e-4, atol=1e-6)

        # Verify that total profile is the exact sum of partial profiles
        np.testing.assert_allclose(aes_total, aes_pt_pt + aes_pt_au + aes_au_au, rtol=1e-9)

    def test_pure_geometric_cross_against_brute_force(self, calculator):
        """
        Validate cross calculation with element='None' (pure geometric Debye scattering,
        form factors f(q) = 1).
        """
        coords_i = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        coords_j = np.array([[0.0, 2.0, 0.0], [0.0, 4.0, 0.0], [1.0, 2.0, 0.0]])

        pos_i = aesdebye.Positions(chemicalSymbols=["None"] * len(coords_i), coordinates=coords_i)
        pos_j = aesdebye.Positions(chemicalSymbols=["None"] * len(coords_j), coordinates=coords_j)

        pdf_cross = calculator.calculatePDF(pos_i, pos_j)
        profile_cross = calculator.calculateIntensity(pdf_cross, 0.5, 6.0, 40, False, 0.4)

        q_vec = np.array(profile_cross.q)
        bf_cross = brute_force_debye_partial(
            coords_i, coords_j, "None", "None", q_vec, calc=None, cross_symmetric=True
        )

        np.testing.assert_allclose(np.array(profile_cross.intensity), bf_cross, rtol=1e-4, atol=1e-6)
