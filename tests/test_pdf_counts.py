import numpy as np
import pytest
import aesdebye


class TestPDFCounts:
    """Sanity tests for PDF pair counts and normalization."""

    def test_same_positions_counts(self, calculator):
        """Verify that samePositions produces N^2 total counts including self-pairs."""
        coords = np.array([
            [0.0, 0.0, 0.0],
            [1.5, 0.0, 0.0],
            [0.0, 2.0, 0.0],
            [1.5, 2.0, 0.0],
        ])
        n_atoms = len(coords)
        pos = aesdebye.Positions(
            chemicalSymbols=["Pt"] * n_atoms,
            coordinates=coords,
        )

        pdf = calculator.calculatePDF(pos)
        assert pdf.testPassed, "pdf.test() must pass for samePositions"
        assert np.isclose(sum(pdf.counts), n_atoms ** 2), (
            f"Expected {n_atoms**2} counts, got {sum(pdf.counts)}"
        )
        assert pdf.counts[0] >= n_atoms, "Bin 0 must include self-pairs"

    def test_cross_positions_counts(self, calculator):
        """Verify that cross positions produce 2 * N_I * N_J total counts."""
        coords_i = np.array([
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
        ])
        coords_j = np.array([
            [0.0, 3.0, 0.0],
            [0.0, 5.0, 0.0],
            [0.0, 7.0, 0.0],
        ])
        n_i = len(coords_i)
        n_j = len(coords_j)

        pos_i = aesdebye.Positions(chemicalSymbols=["Pt"] * n_i, coordinates=coords_i)
        pos_j = aesdebye.Positions(chemicalSymbols=["Au"] * n_j, coordinates=coords_j)

        pdf = calculator.calculatePDF(pos_i, pos_j)
        assert pdf.testPassed, "pdf.test() must pass for cross positions"
        assert np.isclose(sum(pdf.counts), 2 * n_i * n_j), (
            f"Expected {2 * n_i * n_j} counts for bidirectional cross pairs, got {sum(pdf.counts)}"
        )

    def test_multi_element_profile_counts_expansion(self, calculator):
        """
        Verify that a multi-element system satisfies:
        total = a^2 + 2ab + b^2 = (a + b)^2
        """
        n_pt = 3
        n_au = 4
        coords_pt = np.array([[i * 1.5, 0.0, 0.0] for i in range(n_pt)])
        coords_au = np.array([[0.0, j * 2.0, 1.0] for j in range(n_au)])

        coords_all = np.vstack([coords_pt, coords_au])
        symbols_all = ["Pt"] * n_pt + ["Au"] * n_au
        pos_all = aesdebye.Positions(chemicalSymbols=symbols_all, coordinates=coords_all)

        results = calculator.calculateProfile(pos_all, 0.5, 5.0, 20, False, 0.4, "")

        pdf_pt_pt = results["Pt-Pt"][0]
        pdf_au_au = results["Au-Au"][0]
        pdf_pt_au = results["Pt-Au"][0]
        pdf_total = results["total"][0]

        assert pdf_pt_pt.testPassed
        assert pdf_au_au.testPassed
        assert pdf_pt_au.testPassed

        assert np.isclose(sum(pdf_pt_pt.counts), n_pt ** 2)
        assert np.isclose(sum(pdf_au_au.counts), n_au ** 2)
        assert np.isclose(sum(pdf_pt_au.counts), 2 * n_pt * n_au)

        expected_total = (n_pt + n_au) ** 2
        assert np.isclose(sum(pdf_total.counts), expected_total)
        assert np.isclose(
            sum(pdf_pt_pt.counts) + sum(pdf_pt_au.counts) + sum(pdf_au_au.counts),
            sum(pdf_total.counts),
        )

    def test_three_element_system_counts(self, calculator):
        """
        Verify that a 3-element system satisfies:
        (a + b + c)^2 = a^2 + b^2 + c^2 + 2ab + 2bc + 2ac
        """
        n_a, n_b, n_c = 2, 3, 2
        coords_a = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        coords_b = np.array([[0.0, 2.0, 0.0], [1.0, 2.0, 0.0], [2.0, 2.0, 0.0]])
        coords_c = np.array([[0.0, 0.0, 3.0], [1.0, 0.0, 3.0]])

        coords_all = np.vstack([coords_a, coords_b, coords_c])
        symbols_all = ["Pt"] * n_a + ["Au"] * n_b + ["Ag"] * n_c
        pos_all = aesdebye.Positions(chemicalSymbols=symbols_all, coordinates=coords_all)

        results = calculator.calculateProfile(pos_all, 0.5, 5.0, 20, False, 0.4, "")

        expected_total = (n_a + n_b + n_c) ** 2
        total_pdf = results["total"][0]
        assert np.isclose(sum(total_pdf.counts), expected_total)
