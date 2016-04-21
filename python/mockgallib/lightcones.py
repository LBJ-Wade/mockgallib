import mockgallib._mockgallib as c

class LightCones:
    def __init__(self):
        self._lt = c._lightcones_alloc()

    def __len__(self):
        return c._lightcones_len(self._lt)
        
    def load(self, filenames):
        if isinstance(filenames, str):
            c._lightcones_load(self._lt, filenames)
        else:
            for filename in filenames:
                c._lightcones_load(self._lt, filename)
            

    
