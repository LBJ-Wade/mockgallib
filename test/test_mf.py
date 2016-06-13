import unittest
import mockgallib as mock
import math

class TestSigma(unittest.TestCase):
    def setUp(self):
        mock.set_loglevel(3)
        mock.cosmology_set(0.31)
        self.ps = mock.PowerSpectrum('../data/planck1_matterpower.dat')
        self.s = mock.Sigma(self.ps)
        self.mf = mock.MassFunction(self.s)

    def test_mf(self):
        M = [1.0e12, 1.0e13, 1.0e14, 1.0e15]
        expected = [0.034978568962731604,
                    0.0025380630126614636,
                    0.00011791330704897034,
                    8.625465285231354e-07]
        z = 0.5
        for MM, ans in zip(M, expected):
            dndlnM = self.mf.dndlnM(MM, z)
            self.assertAlmostEqual(dndlnM, ans)

if __name__ == '__main__':
    unittest.main()
