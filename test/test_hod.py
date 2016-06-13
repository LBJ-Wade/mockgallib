import unittest
import mockgallib as mock

class TestGrowth(unittest.TestCase):
    def setUp(self):
        mock.set_loglevel(3)
        mock.cosmology_set(0.31)
        self.hod = mock.Hod()
        self.n = len(self.hod)
        self.c = [12.0, 0.5, 0.1, 0.05, 0.95, 0.12, 13.0, 1.0, 1.5, 0.15]
        for i,c in enumerate(self.c):
            self.hod[i] = c

    def test_n(self):
        self.assertEqual(self.n, 10)
        
    def test_get_coef(self):
        for i,c in enumerate(self.c):
            self.assertAlmostEqual(self.hod[i], c)
            self.assertAlmostEqual(self.hod.get_coef()[i], c)

    def test_set_coef(self):
        c2 = [12.5, 1.5, 0.15, 0.15, 1.00, 0.15, 14.0, 1.5, 2.0, 0.20]
        self.hod.set_coef(c2)
        for i,cc in enumerate(c2):
            self.assertAlmostEqual(self.hod[i], cc)

    def test_hod_param(self):
        z = 0.7
        self.assertAlmostEqual(self.hod.logMmin(z), 12.1044)
        self.assertAlmostEqual(self.hod.logM1(z), 13.2)
        self.assertAlmostEqual(self.hod.alpha(z), 1.53)

    def test_ngal1(self):
        z = 0.7
        M = 10**12.0
        self.assertAlmostEqual(self.hod.ncen(M, z), 0.43975709)
        self.assertAlmostEqual(self.hod.nsat(M, z), 0.0)

    def test_ngal2(self):
        z = 1.2
        M = 10**14.0
        self.assertAlmostEqual(self.hod.ncen(M, z), 0.98485394)
        self.assertAlmostEqual(self.hod.nsat(M, z), 2.90460207)

if __name__ == '__main__':
    unittest.main()
