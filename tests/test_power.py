import unittest
import mockgallib as mock

class TestLightCone(unittest.TestCase):
    def setUp(self):
        mock.set_loglevel(3)
        self.ps= mock.PowerSpectrum('../data/planck1_matterpower.dat')

    def test_len(self):
        self.assertEqual(len(self.ps), 907)

    def test_Pk(self):
        self.assertEqual(self.ps.k[0], 0.10000e-03)
        self.assertEqual(self.ps.P[0], 0.45197e+03)

        self.assertEqual(self.ps.k[-1], 0.74031e+04)
        self.assertEqual(self.ps.P[-1], 0.79671e-09)


if __name__ == '__main__':
    unittest.main()
