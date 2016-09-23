import mockgallib._mockgallib as c

class Snapshots:
    """A collection of mock catalogues"""
    def __init__(self):
        self._snps = c._snapshots_alloc()

    def __len__(self):
        return c._snapshots_len(self._snps)

    def __getitem__(self, i):
        return c._snapshots_get(self._snps, i)
        
    def insert(self, fof_filename, part_filename, halo_mass_filename, a_snp):
        c._snapshots_insert(self._snps,
                            fof_filename,
                            part_filename,
                            halo_mass_filename,
                            a_snp[0], a_snp[1], a_snp[2])

    def clear(self):
        c._snapshots_clear(self._snps)
    
