import unittest
import mockgallib as mock

class TestLightcones(unittest.TestCase):
    def setUp(self):
        mock.set_loglevel(3)
        self.lightcones = mock.LightCones()

    def test_load_filename(self):
        self.lightcones.load_h5(['../data/halo_lightcone_mini.h5'])
        self.assertEqual(len(self.lightcones), 1)


if __name__ == '__main__':
    unittest.main()
        
