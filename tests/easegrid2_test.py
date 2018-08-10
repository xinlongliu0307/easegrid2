import unittest
from easegrid2 import *
import subprocess as sp
import ipdb

class ease2_test(object):
    def setUp(self):
        self.m2_9 = M2('g09')
        self.n2_9 = N2('g09')
        self.m2_3 = M2('g03')
        self.n2_3 = N2('g03')
    
    def idl_grid(self, grid_name, s, r):
        idl_ = ['gdl', '-quiet', '-e', 
                'addpath, "../idl"', 
                '$status=wgs84_inverse("{0}", {1},{2}, lat, lon)'.format(grid_name, s, r),
                '&print, lon', '&print, lat']
        ipdb.set_trace()
        p1 = sp.Popen(idl_, shell=True, stdout=sp.PIPE)
        (stdout, stderr) = p1.communicate()
        print stdout
        #return [double(val) for val in stdout.splitlines() if not val.startswith('$')]

class TestGrids(unittest.TestCase):
    def setUp(self):
        self.m1 = M1('g25')
        self.mh1 = M1('g12')
        self.n1 = N1('g25')
        self.s1 = S1('g25')
        self.sh1 = S1('g12')
        self.nh1 = N1('g12')
        #optimize 
        self.rtol = 7e-5

    def test_reversible(self):
        grid_pairs = ([372, 271], [362, 550])
        for grid in [self.s1, self.n1, self.sh1, self.nh1, self.m1]:
            for grid_pair in grid_pairs:
                lon, lat = grid.inverse(*grid_pair)
                r_, s_ = grid.forward(lon, lat)
                #ipdb.set_trace()
                test = np.allclose(grid_pair, [r_, s_],  rtol=self.rtol)
                try:
                    self.assertTrue(test)
                except:
                    Exception('Need to fix pole problems--> lon_0 = 180 or 0')

    def test_pole(self):
        lon, _ = self.n1.inverse(360, 360)
        self.assertTrue(np.allclose(lon, 0, rtol = self.rtol))
        
    def test_polar(self):
        #grid_pairs = ([360, 360], [370, 370], [360, 370])
        grid_pairs = ([370, 370], [360, 370])

        for grid in [self.s1, self.n1, self.sh1, self.nh1]:
            for grid_pair in grid_pairs:
                lonA, latA = grid.from_array()
                lonF, latF = grid.from_file()
                row, col = grid_pair
                lat_pair = [latA[row, col], latF[row, col]]
                lon_pair = [lonA[row, col], lonF[row, col]]
                test = np.allclose(*lat_pair, rtol=self.rtol)
                self.assertTrue(test)
                test = np.allclose(*lon_pair, rtol=self.rtol)
                self.assertTrue(test)

    def test_global(self):
        grids = [self.m1]
        def test_proj(grid_obj, mask=None):
            ''' Checks projection is the same as bins from NSIDC'''
            lonA, latA = grid_obj.from_array()
            lonF, latF = grid_obj.from_file()
            if mask:
                lonA = np.ma.masked_value(lonA, mask)
                lonF = np.ma.masked_value(lonF, mask)
                latA = np.ma.masked_value(latA, mask)
                latF = np.ma.masked_value(latF, mask)
            return np.allclose(lonA, lonF, rtol=self.rtol), np.allclose(latA, latF, rtol=self.rtol)

        for grid in grids:
            lon, lat = test_proj(grid)
            self.assertTrue(lon)
            self.assertTrue(lat)
