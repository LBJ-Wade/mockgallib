import unittest
import mockgallib as mock
import math

class TestGrowth(unittest.TestCase):
    def setUp(self):
        mock.set_loglevel(3)
        mock.cosmology_set(0.31)

    def test_growth(self):
        self.assertAlmostEqual(mock.growth_D(1.0), 1.0)
        self.assertAlmostEqual(mock.growth_D(0.0), 0.0)
        
        # regression test (value computed with this module)
        self.assertAlmostEqual(mock.growth_D(0.5), 0.6084071313593205)


if __name__ == '__main__':
    unittest.main()
