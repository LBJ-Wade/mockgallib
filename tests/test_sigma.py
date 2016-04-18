import unittest
import mockgallib as mock
import math

class TestSigma(unittest.TestCase):
    def setUp(self):
        mock.set_loglevel(3)
        mock.cosmology_set(0.31)
        self.ps = mock.PowerSpectrum('../data/planck1_matterpower.dat')
        self.sigma = mock.Sigma(self.ps)

    def test_len(self):
        self.assertEqual(len(self.sigma), 1001)

    def test_sigma(self):
        # regression test (value computed with this module)
        self.assertAlmostEqual(self.sigma(1.0e10), 3.7489810407407171)
        self.assertAlmostEqual(self.sigma(1.0e12), 2.1281353676997923)
        self.assertAlmostEqual(self.sigma(1.0e15), 0.5396188320554273)

    def test_inv(self):
        # regression test (value computed with this module)
        self.assertAlmostEqual(math.log10(self.sigma.inv(3.748981040740717)), 10.0)
        self.assertAlmostEqual(math.log10(self.sigma.inv(0.2707947263660215)), 16.0)
        self.assertAlmostEqual(math.log10(self.sigma.inv(1.0)), math.log10(7.478911479313975e13))


if __name__ == '__main__':
    unittest.main()
