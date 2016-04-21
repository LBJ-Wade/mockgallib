import unittest
import mockgallib as mock

class TestLightcones(unittest.TestCase):
    def setUp(self):
        self.lightcones = mock.LightCones()

    def test_load_filename(self):
        self.lightcones.load('../data/lightcone_00001.h5')
        self.assertEqual(len(self.lightcones), 1)
        self.lightcones.load('../data/lightcone_00002.h5')
        self.assertEqual(len(self.lightcones), 2)

    def test_load_filenames(self):
        filenames = [ "../data/lightcone_%05d.h5" % i for i in range(1,4) ]
        self.lightcones.load(filenames)
        self.assertEqual(len(self.lightcones), 3)


if __name__ == '__main__':
    unittest.main()
        
