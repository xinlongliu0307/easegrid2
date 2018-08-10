#! /usr/bin/python
'''
Author 
~~~~~~
n.c.steiner, 2013
modified by: kat jensen, 2016

Requires: py-yaml, numpy, pyproj

'''
import os
import sys
from collections import namedtuple

import yaml
import numpy as np
from pyproj import Proj

class _Grid(object):
    '''EASEGrid grid object template.

    Keywords:
    grid_name -- Subgrid name (e.g. g09, g03) 

    ''' 
    _src_path = os.path.split(os.path.abspath(__file__))[0]
    _dat_path = os.path.join(_src_path, 'dat')
    _constant = yaml.load(open(os.path.join(_src_path, 'grids.yaml'), 'r').read())
    coord_scale = 1e-05

    def __init__(self, grid_name):
        # get grid name
        try:
            assert grid_name in self.grid_names
        except:
            raise Exception('Valid grid name in {0}'.format(self.grid_names))
        self.grid_name = grid_name
        # set the projection
        self.p = self.p_args
        # set the grid
        _constant = namedtuple('grid_constant', 'size cols rows r0 s0')
        self.constant = _constant(**self._constant[self.name][grid_name])
        # id string unique to grid/proj combinatio
        self.grid_id = self.name + '_' + getattr(self, 'hemi', 'G') + '_' + self.grid_name
        
    '''
    Properties
    ~~~~~~~~~~
    '''

    @property
    def p(self):
        return self._p

    @p.setter
    def p(self, val):
        ''' Proj4 object.'''
        val = dict(zip(val, [getattr(self, v, None) for v in val]))
        self._p = Proj(val)

    @property
    def geotransform(self):
        ''' Affine transform array.'''
        return [self.constant.s0, self.constant.size, 0, self.constant.r0, 0, self.constant.size]

    @property
    def latitude(self):
        if not hasattr(self, '_latitude'):
            self.set_coords()
        return self._latitude
    
    @property
    def longitude(self):
        if not hasattr(self, '_longitude'):
            self.set_coords()
        return self._longitude

    def set_coords(self):
        try:
            lon, lat = self.from_file()
        except:
            sys.stdout.write('\nNo coordfiles ... calculating coords\n')
            lon, lat = self.from_array()
        self._latitude = lat
        self._longitude = lon

    '''
    Methods
    ~~~~~~~
    '''

    FILES = {'g25': 'L', 'g12': 'H'}
    def from_file(self):
        ''' Not available for all files. '''
        assert self.name[1] == '1'
        assert self.grid_name in self.FILES
        lat_file, lon_file = self._coord_file()
        lat = self._read_coord(lat_file)
        lon = self._read_coord(lon_file)
        return lon, lat 

    def _coord_file(self):
        coord_ = self.hemi + self.FILES[self.grid_name]
        return coord_ + 'LATLSB', coord_ + 'LONLSB'

    def _read_coord(self, file_name):
        file_name = os.path.join(self._dat_path, file_name)
        coord = np.fromfile(file_name, dtype='i')
        coord.resize((self.constant.rows, self.constant.cols))
        return coord * self.coord_scale        
                    
    def from_array(self):
        one_=np.array(1)
        s_, r_ = np.meshgrid(np.arange(1, self.constant.cols + 1), np.arange(1, self.constant.rows + 1))
        return self.inverse(s_, r_)

    def forward(self, lon, lat):
        '''Lat/lon forward projected to raster/scan.
        note: one-based index 
        '''
        x_, y_ = self.p(lon, lat)
        assert self.units == 'm'
        dR, dS = x_/self.constant.size, -y_/self.constant.size 
        r_ = self.constant.r0 + dR + 1
        s_ = self.constant.s0 + dS + 1
        return r_, s_
    
    def inverse(self, r_, s_):
        '''Raster/scan inverse projected to lat/lon.
        note: one-based index 
        '''
        assert self.units == 'm'
        x_ = (r_ - self.constant.r0 - 1) * self.constant.size
        y_ = (s_ - self.constant.s0 - 1) * -self.constant.size
        lon, lat = self.p(x_, y_, inverse=True)
        return lon, lat

'''
Global Cylindrical
~~~~~~~~~~~~~~~~~~
'''

class M1(_Grid):
    ''' EASEGrid ver1 Global, an equal area cylindrical projection.
        
    Keywords:
    grid_name -- Name of grid {'l', 'h'}
    '''
    name = 'M1'
    hemi = 'M'
    proj = 'cea'
    lat_0 = 0
    lon_0 = 0
    lat_ts = 30
    a = 6371228.0
    units = 'm'
    p_args = ['proj', 'lat_0', 'lon_0', 'lat_ts', 'a', 'units']
    grid_names = _Grid._constant[name].keys()

class M2(_Grid):
    ''' EASEGrid ver1 Global, an equal area cylindrical projection.
        
    Keywords:
    grid_name -- Name of grid {'l', 'h'}
    '''
    name = 'M2'
    hemi = 'M'
    proj = 'cea'
    lat_0 = 0
    lon_0 = 0
    lat_1 = 30
    x_0 = 0
    y_0 = 0
    ellps = 'WGS84'
    datum = 'WGS84'
    units = 'm'
    p_args = ['proj', 'lat_0', 'lon_0', 'lat_1', 'x_0', 'y_0', 'ellps', 'datum', 'units']
    grid_names = _Grid._constant[name].keys()

'''
Polar Azimuthal
~~~~~~~~~~~~~~~
'''

class P2(_Grid):
    name = 'P2'
    proj = 'laea'
    lon_0 = 0
    x_0 = 0
    y_0 = 0
    ellps = 'WGS84'
    datum = 'WGS84'
    units = 'm'
    p_args = ['proj', 'lat_0', 'lon_0', 'x_0', 'y_0', 'ellps', 'datum', 'units']
    grid_names = _Grid._constant[name].keys()

class S2(P2):
    lat_0 = -90
    hemi = 'S'
    
class N2(P2):
    lat_0 = 90
    hemi = 'N'

class P1(_Grid):
    name = 'P1'
    proj = 'laea'
    lon_0 = 0 
    x_0 = 0
    y_0 = 0
    a = 6371228.0
    units = 'm'
    p_args = ['proj', 'lat_0', 'lon_0', 'x_0', 'y_0', 'a', 'units']
    grid_names = _Grid._constant[name].keys()

class S1(P1):
    lat_0 = -90
    hemi = 'S'

class N1(P1):
    lat_0 = 90
    hemi = 'N'
