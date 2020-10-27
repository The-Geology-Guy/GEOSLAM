# SLAM Calculator
# Project and method developed by Vincent Cronin (vince_cronin@baylor.edu)
# Python file created and maintained by Luke Pajer (luke.pajer@gmail.com)
# Last Updated: October, 2020

# Data Manipulation and Storage Packages
import pandas as pd
import geopandas as gpd
import utm
import math
import numpy as np
import itertools
import datetime as dt
import xmltodict
import time
from affine import Affine

# File Management Packages
import os
import requests, io
from io import BytesIO
from tqdm import tqdm
import tempfile
from zipfile import ZipFile
from IPython.display import display, Markdown

# File Rasterization Packages
import rasterio as rio
from rasterio.plot import show, show_hist
from IPython.display import GeoJSON, clear_output

# Mapping Packages
import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
import cartopy.io.shapereader as sreader
import cartopy.feature as cfeature

# Plotting Packages
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource
from matplotlib import patheffects
import matplotlib.ticker as ticker
from obspy.imaging.beachball import beach

# Set Max Number of Columns to show
pd.set_option('display.max_columns', 500)

# Set font for all plots
plt.rcParams['font.family'] = ["Helvetica"]

# Get DEM File

class DEM:
    
    def __init__(self, longitude, latitude, **kwargs):
        self.longitude = longitude
        self.latitude = latitude
        self.lon_del = kwargs.get('lon_del', 0.2)
        self.lat_del = kwargs.get('lat_del', 0.2)
    
    class computation:
        pass
    
    def get_bounds(self, **kwargs):
        if self.latitude >= 0:
            lat_min = self.latitude - self.lat_del
            lat_max = self.latitude + self.lat_del
        else:
            lat_min = self.latitude + self.lat_del
            lat_max = self.latitude - self.lat_del

        if self.longitude >= 0:
            lon_min = self.longitude - self.lon_del
            lon_max = self.longitude + self.lon_del
        else:
            lon_min = self.longitude + self.lon_del
            lon_max = self.longitude - self.lon_del

        return lat_min, lat_max, lon_min, lon_max
    
    def get_dem(self, lat_min, lat_max, lon_min, lon_max, **kwargs):
        computation = DEM.computation()
        display(Markdown('##### Begin Downloading DEM File...'))
        tic = time.perf_counter()
        
        dem_type = kwargs.get('dem_type', 'SRTMGL1_E')
        output_type = kwargs.get('output_type', 'GTiff') 
        
        url = 'https://portal.opentopography.org/API/globaldem'

        query = {
            'demtype': dem_type,
            'south': str(lat_min),
            'north': str(lat_max),
            'west': str(lon_max),
            'east': str(lon_min),
            'outputFormat': 'GTiff'
        }

        fd, path = tempfile.mkstemp()
        
        try:
            with os.fdopen(fd, 'wb') as f:
                ret = requests.get(url, stream=True, params=query)
                for data in ret.iter_content(1024):
                    f.write(data)

            src = rio.open(path)

        finally:
            os.remove(path)
            
        toc = time.perf_counter()
        display(Markdown(f"##### Downloaded DEM file in {toc - tic:0.1f} seconds"))
        display(Markdown('##### Processing...'))
        
        zone = utm.from_latlon(self.latitude, self.longitude)[2]
        src1 = rio.vrt.WarpedVRT(src, crs=f'EPSG:326{str(zone)}')
        cellsize = (src.transform[0] * (111132.0894 * math.cos(np.deg2rad(self.latitude))))
        bounds = src.bounds
        
        computation.src_meters = src1
        computation.bounds = bounds
        computation.cellsize = cellsize
        
        elevation = src.read(1)
        computation.elevation = elevation
        display(Markdown('##### Complete.'))
        
        return computation

# Get Fault Lines from the Quaternary Faults and Folds Database

class qFaults:
    
    def get_qfaults(**kwargs):
        # https://earthquake.usgs.gov/static/lfs/nshm/qfaults/
        display(Markdown("##### _Warning: Be courteous to the USGS API and do not download this file more than once, as this is a very large file._"))
        
        qfaults_file = kwargs.get('qfaults_file', 'Qfaults_GIS')

        tic = time.perf_counter()

        resp = requests.get(f'https://earthquake.usgs.gov/static/lfs/nshm/qfaults/{qfaults_file}.zip')

        path = tempfile.mkdtemp()

        path0 = os.path.join(path, f'{qfaults_file}.shx')
        path1 = os.path.join(path, f'{qfaults_file}.shp')

        bytes_ = BytesIO(resp.content)

        try:
            with ZipFile(bytes_, 'r') as zipObj:

                listOfFileNames = zipObj.namelist()

                for fileName in listOfFileNames:

                    if fileName.endswith('.shx'):

                        with zipObj.open(fileName) as myfile:
                            with open(path0, 'wb') as f:
                                f.write(zipObj.open(fileName).read())

                    if fileName.endswith('.shp'):

                        with zipObj.open(fileName) as myfile1:
                            with open(path1, 'wb') as f1:
                                f1.write(zipObj.open(fileName).read())

                shp = gpd.read_file(path1, driver='ESRI Shapefile')

        finally:
            os.remove(path0)
            os.remove(path1)

        toc = time.perf_counter()
        display(Markdown(f"##### Downloaded Quaternary Faults and Folds Shape File in {toc - tic:0.1f} seconds"))

        return shp

# Query the USGS Database for Earthquake Data

class get_data:
    
    def __init__(self, **kwargs):
        self.lon_min = kwargs.get('lon_min', '0.0')
        self.lon_max = kwargs.get('lon_max', '0.0')
        self.lat_min = kwargs.get('lat_min', '0.0')
        self.lat_max = kwargs.get('lat_max', '0.0')
    
    class computation:
        pass
    
    def event_query(self, **kwargs):
        # Data Pulled from https://earthquake.usgs.gov/fdsnws/event/1/
        tic = time.perf_counter()
        
        lat_min = kwargs.get('lat_min', self.lat_min)
        lat_max = kwargs.get('lat_max', self.lat_max)
        lon_min = kwargs.get('lon_min', self.lon_min)
        lon_max = kwargs.get('lon_max', self.lon_max)

        start_ = kwargs.get('start_time', None)
        end_ = kwargs.get('end_time', None)
        eventid_ = kwargs.get('event_id', None)
        mindepth_ = kwargs.get('min_depth', None)
        maxdepth_ = kwargs.get('max_depth', None)
        minmagnitude_ = kwargs.get('min_mag', None)
        maxmagnitude_ = kwargs.get('max_mag', None)
        eventtype_ = kwargs.get('event_type', None)

        url_geo = 'https://earthquake.usgs.gov/fdsnws/event/1/query'
        query_geo = {
            'format': 'geojson',
            'starttime': start_,
            'endtime': end_,
            'minlatitude': str(lat_min),
            'maxlatitude': str(lat_max),
            'minlongitude': str(lon_max),
            'maxlongitude': str(lon_min),
            'eventid': eventid_,
            'mindepth': mindepth_,
            'maxdepth': maxdepth_,
            'minmagnitude': minmagnitude_,
            'maxmagnitude': maxmagnitude_,
            'orderby': 'time',
            'eventtype': eventtype_,
        }

        geo = requests.get(url_geo, stream=True, params=query_geo)

        data = geo.json()

        df_geo = gpd.GeoDataFrame.from_features(data['features'])

        earthquakes = df_geo[['title', 'mag', 'place', 'time', 'updated', 
                              'tz', 'url', 'detail', 'cdi', 'mmi', 'sig', 
                              'types', 'nst', 'dmin', 'rms', 'gap', 'geometry']].dropna(how='all', axis='columns')

        earthquakes['time'] = pd.to_datetime(earthquakes['time'], unit='ms')

        earthquakes['updated'] = pd.to_datetime(earthquakes['updated'], unit='ms')

        toc = time.perf_counter()
        display(Markdown(f"##### Downloaded Regional Earthquake Data in {toc - tic:0.4f} seconds"))

        return earthquakes
    
    def focal_data(self, earthquakes, **kwargs):
        computation = get_data.computation()
        
        tic = time.perf_counter()
        
        def get_fm_detail(eqs, **kwargs):
            eqs = eqs.reset_index(drop=False)
            fil = kwargs.get('filter_rows', 'no_filter')
            update_ = kwargs.get('update_rows', 'no_update')
            data_ = []
            for i in range(len(eqs)):
                det = requests.get(eqs.detail[i])
                det_check = det.json()["properties"]["products"].get("focal-mechanism", 'no-fm')
                code = det.json()['properties']['code']
                it = itertools.product('ABCDEFGHIJKLMNOPQRSTUVWXYZ', repeat=2)
                li = []
                for y in it:
                    li.append(''.join(y))

                if det_check is not 'no-fm':
                    fms = det.json()["properties"]["products"]["focal-mechanism"]
                    k = eqs['index'][i]
                    if k <= 25:
                        l = chr(97 + k)
                    else:
                        l = li[i - 26]
                    for j in range(len(fms)):
                        if len(fms) > 1:
                            index = str(l).upper() + str(j+1)
                        else:
                            index = str(l).upper()
                        mag, time, geom = eqs['mag'][i], eqs['time'][i], eqs['geometry'][i]
                        updated = pd.to_datetime(fms[j]["updateTime"], unit='ms')
                        if (fms[j]["properties"].get("longitude", False) or fms[j]["properties"].get("latitude", False)):
                            lon = fms[j]["properties"]["longitude"]
                            lat = fms[j]["properties"]["latitude"]
                        else:
                            lon = 0
                            lat = 0
                        
                        try:
                            strike = fms[j]["properties"]["nodal-plane-1-strike"]
                            dip = fms[j]["properties"]["nodal-plane-1-dip"]
                            rake = fms[j]["properties"]["nodal-plane-1-rake"]
                        except:
                            strike = 0
                            dip = 0
                            rake = 0
                            depth = 0
                            
                        try:
                            depth = fms[j]["properties"]["depth"]
                        except:
                            depth = 0
                            
                        try:
                            strike2 = fms[j]["properties"]["nodal-plane-2-strike"]
                            dip2 = fms[j]["properties"]["nodal-plane-2-dip"]
                            rake2 = fms[j]["properties"]["nodal-plane-2-rake"]
                        except:
                            strike2 = 0
                            dip2 = 0
                            rake2 = 0

                        quake = fms[j]['contents']['quakeml.xml']['url']
                        response = requests.get(quake)
                        quake_data = response.content

                        fm_data = xmltodict.parse(quake_data)["q:quakeml"]["eventParameters"]["event"]["focalMechanism"]
                        try:
                            strike_u = fm_data["nodalPlanes"]["nodalPlane1"]["strike"]["uncertainty"]
                            dip_u = fm_data["nodalPlanes"]["nodalPlane1"]["dip"]["uncertainty"]
                            rake_u = fm_data["nodalPlanes"]["nodalPlane1"]["rake"]["uncertainty"]

                            strike2_u = fm_data["nodalPlanes"]["nodalPlane2"]["strike"]["uncertainty"]
                            dip2_u = fm_data["nodalPlanes"]["nodalPlane2"]["dip"]["uncertainty"]
                            rake2_u = fm_data["nodalPlanes"]["nodalPlane2"]["rake"]["uncertainty"]

                        except:
                            strike_u, dip_u, rake_u = 0, 0, 0
                            strike2_u, dip2_u, rake2_u = 0, 0, 0

                        column_names = ['ID', 'mag', 'code', 'time', 'updated', 'longitude', 
                                        'latitude', 'depth', 'np1_strike', 'np1_dip', 'np1_rake', 
                                        'np1S_uncert', 'np1D_uncert', 'np1R_uncert',
                                        'np2_strike', 'np2_dip', 'np2_rake', 'np2S_uncert', 
                                        'np2D_uncert', 'np2R_uncert','geometry']
                        column_value = [index, mag, code, time, updated, lon, lat, depth, 
                                        strike, dip, rake, strike_u, dip_u, rake_u, 
                                        strike2, dip2, rake2, strike2_u, dip2_u, rake2_u, geom]
                        data_.append(dict(zip(column_names, column_value)))
                else:
                    k = eqs['index'][i]
                    if k <= 25:
                        l = chr(97 + k)
                    else:
                        l = li[i - 26]
                    index = str(l).upper()
                    title, mag, time, geom = eqs['title'][i], eqs['mag'][i], eqs['time'][i], eqs['geometry'][i]
                    updated, lon, lat, depth, strike, dip, rake = 0, 0, 0, 0, 0, 0, 0
                    strike2, dip2, rake2, strike_u, dip_u, rake_u = 0, 0, 0, 0, 0, 0
                    strike2_u, dip2_u, rake2_u = 0, 0, 0
                    column_names = ['ID', 'mag', 'code', 'time', 'updated', 'longitude', 
                                    'latitude', 'depth', 'np1_strike', 'np1_dip', 'np1_rake', 
                                    'np1S_uncert', 'np1D_uncert', 'np1R_uncert',
                                    'np2_strike', 'np2_dip', 'np2_rake', 'np2S_uncert', 
                                    'np2D_uncert', 'np2R_uncert','geometry']
                    column_value = [index, mag, code, time, updated, lon, lat, depth, 
                                    strike, dip, rake, strike_u, dip_u, rake_u, 
                                    strike2, dip2, rake2, strike2_u, dip2_u, rake2_u, geom]
                    data_.append(dict(zip(column_names, column_value)))

            data_df = pd.DataFrame(data_)
            if update_ is not 'no_update':
                data_T = data_df.set_index('ID').T
                for key, value in test_json.items():
                    vals = [data_T[key]['title'], data_T[key]['mag'], data_T[key]['code'], data_T[key]['time']]
                    ends = [data_T[key]['geometry']]
                    value_ = vals + value + ends
                    new_df = pd.DataFrame(value_, columns=[key], index=data_T.index)
                    data_T.update(new_df)

                data_df = data_T.T.reset_index(drop=False)

            if fil is not 'no_filter':
                data_df = data_df[data_df['ID'].isin(fil)].reset_index(drop=True)

            return data_df
        
        eqs_fm = get_fm_detail(eqs=earthquakes)
        computation.focal_data = eqs_fm

        def get_err_detail(eqs, **kwargs):
            eqs = eqs.reset_index(drop=False)
            fil = kwargs.get('filter_rows', 'no_filter')
            update_ = kwargs.get('update_rows', 'no_update')
            data_ = []
            for i in range(len(eqs)):
                det = requests.get(eqs.detail[i])
                det_check = det.json()["properties"]["products"].get("origin", 'no-fm')
                it = itertools.product('ABCDEFGHIJKLMNOPQRSTUVWXYZ', repeat=2)
                li = []
                for y in it:
                    li.append(''.join(y))

                if det_check is not 'no-fm':
                    fms = det.json()["properties"]["products"]["origin"]
                    k = eqs['index'][i]
                    if k <= 25:
                        l = chr(97 + k)
                    else:
                        l = li[i - 26]
                    for j in range(len(fms)):
                        if len(fms) > 1:
                            index = str(l).upper() + str(j+1)
                        else:
                            index = str(l).upper()
                        title, mag, time, geom = eqs['title'][i].replace(",", ""), eqs['mag'][i], eqs['time'][i], eqs['geometry'][i]
                        updated = pd.to_datetime(fms[j]["updateTime"], unit='ms')
                        if (fms[j]["properties"].get("longitude", False) or fms[j]["properties"].get("latitude", False)):
                            lon = fms[j]["properties"]["longitude"]
                            lat = fms[j]["properties"]["latitude"]
                        else:
                            lon = 0
                            lat = 0
                        err_check = fms[j]["properties"].get("error-ellipse-plunge", 'no-or')
                        if err_check is not 'no-or':
                            code = fms[j]["properties"]["eventsourcecode"]
                            eplunge = fms[j]["properties"]["error-ellipse-plunge"]
                            eh1 = fms[j]["properties"]["error-ellipse-major"]
                            eh2 = fms[j]["properties"]["error-ellipse-minor"]
                            ehaz = fms[j]["properties"]["error-ellipse-azimuth"]
                            ez = fms[j]["properties"]["vertical-error"]
                        else:
                            code = fms[j]["properties"]["eventsourcecode"]
                            eplunge = 0
                            eh1 = 0
                            eh2 = 0
                            ehaz = 0
                            ez = 0
                        data_.append(dict(zip(['ID', 'title', 'mag', 'code', 'time', 'updated', 'longitude', 'latitude', 'eplunge', 'eh1', 'eh2', 'ehaz', 'ez', 'geometry'],
                                                  [index, title, mag, code, time, updated, lon, lat, eplunge, eh1, eh2, ehaz, ez, geom])))

            data_df = pd.DataFrame(data_)
            if update_ is not 'no_update':
                data_T = data_df.set_index('ID').T
                for key, value in test_json.items():
                    vals = [data_T[key]['title'], data_T[key]['mag'], data_T[key]['code'], data_T[key]['time']]
                    ends = [data_T[key]['geometry']]
                    value_ = vals + value + ends
                    new_df = pd.DataFrame(value_, columns=[key], index=data_T.index)
                    data_T.update(new_df)

                data_df = data_T.T.reset_index(drop=False)

            if fil is not 'no_filter':
                data_df = data_df[data_df['ID'].isin(fil)].reset_index(drop=True)

            return data_df
        
        eqs_err = get_err_detail(eqs=earthquakes)
        computation.error_data = eqs_err

        def full_data(details, errors, **kwargs):
            try:
                err = errors.drop(['ID', 'title', 'mag', 'time', 'updated', 'longitude', 'latitude', 'geometry'], axis=1)
                order_ = ['ID', 'mag', 'code', 'time', 'updated', 'longitude', 'latitude', 'depth', 
                          'np1_strike', 'np1_dip', 'np1_rake', 'np1S_uncert', 'np1D_uncert', 'np1R_uncert', 
                          'np2_strike', 'np2_dip', 'np2_rake', 'np2S_uncert', 'np2D_uncert', 'np2R_uncert', 
                          'eplunge', 'eh1', 'eh2', 'ehaz', 'ez', 'geometry']
                x = pd.merge(details, err, on='code', how='left')
                x = x[order_]
            except:
                display(Markdown('#### Insufficient Data From API.'))
                x = details
            return x
        
        full = full_data(details=eqs_fm, errors=eqs_err)
        
        computation.data = full
        
        toc = time.perf_counter()
        display(Markdown(f"##### Successfully Queried Events in {toc - tic:0.1f} seconds"))
        
        return computation

# Generate Swath Data

class swaths:

    def __init__(self, elevation):
        self.elevation = elevation
    
    class computation:
        pass
    
    def grid_adjustment(self, longitude, latitude, zoneMeridian):
        gridNorthAdjustment = (zoneMeridian - (longitude)) * (np.sin(np.deg2rad(latitude)))
        display(Markdown(f'- **Grid North Adjustment**: {gridNorthAdjustment:0.2f}$^\circ$'))
        return gridNorthAdjustment
    
    def makeVector(plunge, trend):
        a_ = np.cos(np.deg2rad(plunge)) * np.sin(np.deg2rad(trend))
        b_ = np.cos(np.deg2rad(plunge)) * np.cos(np.deg2rad(trend))
        c_ = -np.sin(np.deg2rad(plunge))
        vec = np.array([a_,b_,c_])
        return vec

    def vectorNorm(x):
        return math.sqrt(np.dot(x, x))

    def unitVector(x):
        a_ = (x[0] / swaths.vectorNorm(x))
        b_ = (x[1] / swaths.vectorNorm(x))
        c_ = (x[2] / swaths.vectorNorm(x))
        uVe = np.array([a_,b_,c_])
        return uVe

    def vectorAngle(a, b):
        a_b = np.dot(a, b)
        vec = swaths.vectorNorm(a) * swaths.vectorNorm(b)
        return np.arccos(np.deg2rad(a_b/vec)) 
    
    def pointEvaluator(xCoord, yCoord, zCoord, width, fltUNrml):
        locVect = np.array([xCoord, yCoord, zCoord])
        distToFlt = abs(np.dot(fltUNrml, locVect))
        result = 10 if (distToFlt <= width) else 0
        return result
    
    def findHalfWidth(dipAz, dipAng, hEr1, hEr1Az, hEr2, vEr):
        var1 = hEr1Az - 90
        var2 = np.array([np.cos(np.deg2rad(dipAng)) * np.sin(np.deg2rad(dipAz)), np.cos(np.deg2rad(dipAng)) * np.cos(np.deg2rad(dipAz)), -np.sin(np.deg2rad(dipAng))])
        var3 = np.array([np.sin(np.deg2rad(dipAz - 90)), np.cos(np.deg2rad(dipAz - 90)), 0])
        var4 = np.array([[np.cos(np.deg2rad(var1)), -np.sin(np.deg2rad(var1)), 0], [np.sin(np.deg2rad(var1)), np.cos(np.deg2rad(var1)), 0], [0, 0, 1]])
        var5 = np.cross(var2, var3)
        var6 = np.dot(var4, var5)
        var7 = (1 / (hEr1**2)) * (var6[0] * ((hEr1**2) / 2))**2
        var8 = (1 / (hEr2**2)) * (var6[1] * ((hEr2**2) / 2))**2
        var9 = (1 / (vEr**2)) * (var6[2] * ((vEr**2) / 2))**2
        var10 = math.sqrt(1 / (var7 + var8 + var9))
        var11 = np.array([(var10 * var6[0] * ((hEr1**2) / 2)), (var10 * var6[1] * ((hEr2**2) / 2)), (var10 * var6[2] * ((vEr**2) / 2))])
        var12 = np.dot(var6, var11)
        return var12
    
    def nodal_calc(data, width, fltUNrml, lon, lat, cellsize, cellsize_res, focalDepth, bounds):
        utmCen = utm.from_latlon(lat, lon)
        if (bounds[0] < 360) and (bounds[1] < 360):
            bounds1 = utm.from_latlon(bounds[1], bounds[0])
            bounds2 = utm.from_latlon(bounds[3], bounds[2])
            bounds = [bounds1[0], bounds1[1], bounds2[0], bounds2[1]]
        elev = data.copy()
        for i in range(data.shape[0]):
            yCoord = ((cellsize_res) * (data.shape[0] - i)) - (utmCen[1] - bounds1[1])
            for j in range(data.shape[1]):
                xCoord = ((cellsize) * (j)) - (utmCen[0] - bounds1[0])
                zCoord = elev[i][j] - focalDepth
                elev[i][j] = swaths.pointEvaluator(xCoord, yCoord, zCoord, width, fltUNrml)
                j += 1
            i += 1
        return elev
    
    def focal_metrics(self, lon, lat, depth):
        # Focal Latitude
        focalLat = float(lat)

        # Focal Longitude
        focalLong = float(lon)

        # Focal Depth in Meters
        focalDepthKm = float(depth)
        focalDepth = focalDepthKm * (-1000)
        
        display(Markdown(f'- **Latitude**: {focalLat:0.3f}$^\circ$\n- **Longitude**: {focalLong:0.3f}$^\circ$\n- **Depth**: {focalDepth:0.1f} meters'))
        
        return [focalLat, focalLong, focalDepth]
    
    def err_computations(self, eh1, eh1Az, eh2, ez, grid_adjust):
        # Error calculations
        eh1 = (float(eh1))
        eh1Azimuth = float(eh1Az)
        if (eh1Azimuth + grid_adjust) < 0:
            eh1Az = 360 + (eh1Azimuth + grid_adjust)
        elif (eh1Azimuth + grid_adjust) > 360:
            eh1Az = (eh1Azimuth + grid_adjust) - 360
        else:
            eh1Az = eh1Azimuth + grid_adjust
        eh2 = (float(eh2))
        ez = float(ez) * 1000
        
        display(Markdown(f'- **eh1**: {eh1:0.1f} meters\n- **eh1 Azimuth**: {eh1Az:0.1f}$^\circ$\n- **eh2**: {eh2:0.1f} meters\n- **ez**: {ez:0.1f} meters'))
        
        return eh1, eh1Az, eh2, ez
    
    def light_direction(self, npDipTr):
        npDipTrRad = npDipTr * (math.pi / 180)

        if npDipTrRad > ((3 * math.pi) / 2):
            lightDirection = ((2 * math.pi) - npDipTrRad) * (180 / math.pi)
        else:
            lightDirection = ((math.pi / 2) - npDipTrRad) * (180 / math.pi)
            
        display(Markdown(f'- **Light Direction**: {lightDirection:0.2f}$^\circ$'))
            
        return lightDirection
    
    def nodal_comps(self, strike, grid_adjust, strike_uncert, plunge, dip_ang):
        npDipTrend_ = float(strike) + 90
        if npDipTrend_ < 0:
            npDipTrend = npDipTrend_ + 360
        elif npDipTrend_ > 360:
            npDipTrend = npDipTrend_ - 360
        else:
            npDipTrend = npDipTrend_
            
        npTest = (npDipTrend + grid_adjust)
        if npTest < 0:
            npDipTr = npTest + 360
        elif npTest > 360:
            npDipTr = npTest - 360
        else:
            npDipTr = npTest

        npDipTrendUncert = float(strike_uncert)
        npDipPlunge = float(plunge)
        npDipAngUncert = float(dip_ang)
        
        display(Markdown(f'- **Dip Trend**: {npDipTrend:0.1f}$^\circ$\n- **Dip Trend + North Grid Adjustment**: {npDipTr:0.1f}$^\circ$\n- **Dip Trend Uncertainty**: {npDipTrendUncert:0.1f}$^\circ$\n- **Dip Plunge**: {npDipPlunge:0.1f}$^\circ$\n- **Dip Angle Uncertainty**: {npDipAngUncert:0.1f}$^\circ$'))

        return [npDipTrend, npDipTr, npDipTrendUncert, npDipPlunge, npDipAngUncert]
    
    def dipAngle(n, npDipPlunge, npDipAngUncert):
    
        if n == 0:
            # Swath 1
            angle = npDipPlunge
        elif n < 3:
            # Swaths 2 & 3
            if (npDipPlunge - npDipAngUncert) < 0:
                angle = abs(npDipPlunge - npDipAngUncert)
            else:
                angle = npDipPlunge - npDipAngUncert
        elif n < 6:
            # Swaths 4, 5, & 6
            if (npDipPlunge + npDipAngUncert) > 90:
                angle = 180 - (npDipPlunge + npDipAngUncert)
            else:
                angle = npDipPlunge + npDipAngUncert
        else:
            # Swath 7
            angle = npDipPlunge

        return angle
    
    def dipTrend(n, npDipPlunge, npDipTr, npDipTrendUncert, npDipAngUncert, **kwargs):
        
        if n == 0:
            azimuth = npDipTr
        elif n == 1:
            # Swath 2
            if (npDipPlunge - npDipAngUncert) < 0:
                if (npDipTr - npDipTrendUncert) > 180:
                    azimuth = (npDipTr - npDipTrendUncert) - 180
                else:
                    azimuth = (npDipTr - npDipTrendUncert) + 180
            else:
                if (npDipTr - npDipTrendUncert) < 0:
                    azimuth = 360 + (npDipTr - npDipTrendUncert)
                else:
                    azimuth = (npDipTr - npDipTrendUncert)
        elif n == 2:
            # Swath 3
            if (npDipPlunge - npDipAngUncert) < 0:
                if (npDipTr + npDipTrendUncert) > 180:
                    azimuth = (npDipTr + npDipTrendUncert) - 180
                else:
                    azimuth = (npDipTr + npDipTrendUncert) + 180
            else:
                if (npDipTr + npDipTrendUncert) > 360:
                    azimuth = (npDipTr + npDipTrendUncert) - 360
                else:
                    azimuth = (npDipTr + npDipTrendUncert)
        elif n == 3:
            # Swath 4
            if (npDipPlunge + npDipAngUncert) > 90:
                if (npDipTr - npDipTrendUncert) > 180:
                    azimuth = (npDipTr - npDipTrendUncert) - 180
                else:
                    azimuth = (npDipTr - npDipTrendUncert) + 180
            else:
                if (npDipTr - npDipTrendUncert) < 0:
                    azimuth = 360 + (npDipTr - npDipTrendUncert)
                else:
                    azimuth = (npDipTr - npDipTrendUncert)
        elif n == 4:
            # Swath 5
            if (npDipPlunge + npDipAngUncert) > 90:
                if (npDipTr + npDipTrendUncert) > 180:
                    azimuth = (npDipTr + npDipTrendUncert) - 180
                else:
                    azimuth = (npDipTr + npDipTrendUncert) + 180
            else:
                if (npDipTr + npDipTrendUncert) > 360:
                    azimuth = (npDipTr + npDipTrendUncert) - 360
                else:
                    azimuth = (npDipTr + npDipTrendUncert)
        elif n == 5:
            # Swath 6
            if (npDipPlunge + npDipAngUncert) > 90:
                if npDipTr > 180:
                    azimuth = npDipTr - 180
                else:
                    azimuth = npDipTr + 180
            else:
                azimuth = npDipTr
        else:
            # Swath 7
            azimuth = npDipTr

        return azimuth
    
    def zoneWidth(n, widthFactor, cellsize, angle, azimuth, eh1, eh1Az, eh2, ez, **kwargs):
        multiplier = kwargs.get('multiplier', 1)
        if n==0:
            # Swath 1
            width = widthFactor * (cellsize) * (np.sin(np.deg2rad(angle)))
        else:
            # Swaths 2 through 7
            width = swaths.findHalfWidth(azimuth, angle, eh1, eh1Az, eh2, ez)

        return width * multiplier
    
    def swath_calc(self, npDipPlunge, npDipAngUncert, npDipTr, npDipTrendUncert,
                   widthFactor, cellsize, cellsize_res, eh1, eh1Az, eh2, ez, lat, lon, 
                   focalDepth, bounds, **kwargs):
        
        computation = swaths.computation()
        
        data = self.elevation
        cellsize_res = cellsize_res
        
        tic = time.perf_counter()
        multiplier = kwargs.get('multiplier', 1)
        k = kwargs.get('k', 7)
        thin = kwargs.get('thin', False)
        i_ = kwargs.get('i', 0)
        display(Markdown("#### Evaluate Ground Surface Traces"))
        i = i_
        while i < k:
            tic1 = time.perf_counter()

            angle = swaths.dipAngle(n=i, npDipPlunge=npDipPlunge, npDipAngUncert=npDipAngUncert)

            azimuth = swaths.dipTrend(n=i, npDipPlunge=npDipPlunge, npDipAngUncert=npDipAngUncert,
                                      npDipTr=npDipTr, npDipTrendUncert=npDipTrendUncert)

            width = swaths.zoneWidth(n=i, widthFactor=widthFactor, cellsize=cellsize, angle=angle, 
                              azimuth=azimuth, eh1=eh1, eh1Az=eh1Az, eh2=eh2, ez=ez, multiplier=multiplier)

            dipVect = swaths.unitVector(swaths.makeVector(angle, azimuth))

            if azimuth < 90:
                npStrike = azimuth + 270
            else:
                npStrike = azimuth - 90

            strikeVect = np.array([np.sin(np.deg2rad(npStrike)), np.cos(np.deg2rad(npStrike)), 0])
            fltUNrml = swaths.unitVector(np.cross(dipVect, strikeVect))

            if thin == False:
                data0 = data.copy()
            else:
                data0 = data.copy() / data.min()
            
            if i == 0:
                data0 = data.copy() / data.min()

            elev0 = swaths.nodal_calc(data=data0, width=width, fltUNrml=fltUNrml, lat=lat, lon=lon,
                                      cellsize=cellsize, cellsize_res=cellsize_res, 
                                      focalDepth=focalDepth, bounds=bounds)
            if i == 0:
                computation.middle_road = np.flip(elev0, 0)

            if i == i_:
                elev = elev0
            else:
                elev = elev + elev0

            toc0 = time.perf_counter()
            display(Markdown("- Swath " + str(i + 1) + " of " + str(k) + f" Completed in {toc0 - tic1:0.1f} seconds"))
            i += 1
        toc = time.perf_counter()
        seconds = (toc - tic) % (24 * 3600) 
        hour = seconds // 3600
        seconds %= 3600
        minutes = seconds // 60
        seconds %= 60
        display(Markdown("#### Computed Ground Surface Traces in %01d minutes and %02d seconds." % (minutes, seconds)))
        computation.swaths = np.flip(elev, 0)

        return computation
    
    def fill_swaths(self, swaths, **kwargs):
        
        test = swaths.copy()
            
        li_max = []
        li_min = []
        k = 0
        for i in range(test.shape[0]):
            li_temp = []
            for j in range(test.shape[1]):
                li = test[i][j]
                if li > 0:
                    li_temp.append(j)
                j += 1
            try:
                li_max.append(max(li_temp))
                li_min.append(min(li_temp))
            except:
                li_max.append(-1)
                li_min.append(-1)
            i += 1

        r = np.arange(test.shape[0])[:,None]
        swaths = (((li_min <= r) & (r <= li_max)).astype(int)) * 10
        return swaths.T

    def fill_corners(self, swaths, corners_fill, **kwargs):
        if (corners_fill == 'NE'):
            n_ = range(0, 1)
            row_ = [(swaths.shape[1] - 1), 0] 
        elif (corners_fill == 'SW'):
            n_ = range(1, 2)
            row_ = [(swaths.shape[1] - 1), 0] 
        elif (corners_fill == 'NW'):
            n_ = range(2, 3)
            row_ = [0, 0, (swaths.shape[1] - 1), 0]
            swaths = np.flip(swaths, 1)
        elif (corners_fill == 'SE'):
            n_ = range(3, 4)
            row_ = [0, 0, (swaths.shape[1] - 1), 0]
            swaths = np.flip(swaths, 1)
        elif (corners_fill == 'N') or (corners_fill == 'S'):
            n_ = range(0, 4)
            row_ = [(swaths.shape[1] - 1), 0, (swaths.shape[1] - 1), 0]

        for n in n_:
            for i in range(2):
                if i == 0:
                    test = swaths.copy()
                    li_max = []
                    li_temp = []
                    for i in range(test.shape[0]):
                        li = test[i][row_[n]]
                        if li > 0:
                            li_temp.append(i)
                        i += 1
                    if (n == 0) or (n == 2):
                        li_max.append(max(li_temp))
                    elif (n == 1) or (n == 3):
                        li_max.append(min(li_temp))
                elif i == 1:
                    test = swaths.copy().T
                    li_max2 = []
                    li_temp2 = []
                    for i in range(test.shape[1]):
                        li = test[i][row_[n]]
                        if li > 0:
                            li_temp2.append(i)
                        i += 1
                    if (n == 0) or (n == 2):
                        li_max2.append(max(li_temp2))
                    elif (n == 1) or (n == 3):
                        li_max2.append(min(li_temp2))
                else:
                    break

            swath = swaths.copy()
            if n == 0:
                for i in range(li_max[0], swaths.shape[0]):
                    for j in range(li_max2[0], swaths.shape[1]):
                        swath[i][j] = 10
                        j += 1
                    i += 1
            elif n == 1:
                for i in range(li_max[0], -1, -1):
                    for j in range(li_max2[0], -1, -1):
                        swath[i][j] = 10
                        j += 1
                    i += 1
            elif n == 2:
                for i in range(li_max[0], swaths.shape[0]):
                    for j in range(li_max2[0], swaths.shape[1]):
                        swath[i][j] = 10
                        j += 1
                    i += 1
            elif (n == 3):
                for i in range(li_max[0], -1, -1):
                    for j in range(li_max2[0], -1, -1):
                        swath[i][j] = 10
                        j += 1
                    i += 1

            try:
                swath_ = swath_ + swath
            except:
                swath_ = swaths + swath
        
        if (corners_fill == 'NW') or (corners_fill == 'SE'):
            swath_ = np.flip(swath_, 1)
        
        return swath_
    
    def fill_middle(swaths):
        for n in range(0, 4):
            li_max = []
            for i in range(4):       
                if i == 0:
                    test = swaths.copy()
                    li_temp = []
                    mid_ = int(test.shape[0] / 2)
                    end_ = test.shape[0]
                    for i in range(mid_, end_):
                        li = test[end_ - 1][i]
                        if li > 0:
                            li_temp.append(i)
                        i += 1
                    try:
                        li_max.append(min(li_temp))
                        li_max.append(max(li_temp))
                    except:
                        li_max.append(-1)
                        li_max.append(-1)

                    li_temp2 = []
                    for i in range(0, mid_):
                        li = test[end_ - 1][i]
                        if li > 0:
                            li_temp2.append(i)
                        i += 1
                    try:
                        li_max.append(max(li_temp2))
                        li_max.append(min(li_temp2))
                    except:
                        li_max.append(-1)
                        li_max.append(-1)

                elif i == 1:
                    test = np.flip(test, 1)
                    li_temp = []
                    mid_ = int(test.shape[0] / 2)
                    end_ = test.shape[0]
                    for i in range(mid_, end_):
                        li = test[0][i]
                        if li > 0:
                            li_temp.append(i)
                        i += 1
                    try:
                        li_max.append(min(li_temp))
                        li_max.append(max(li_temp))
                    except:
                        li_max.append(-1)
                        li_max.append(-1)

                    li_temp2 = []
                    for i in range(0, mid_):
                        li = test[0][i]
                        if li > 0:
                            li_temp2.append(i)
                        i += 1
                    try:
                        li_max.append(max(li_temp2))
                        li_max.append(min(li_temp2))
                    except:
                        li_max.append(-1)
                        li_max.append(-1)

                elif i == 2:
                    test = np.flip(test, 0)
                    mid_ = int(test.shape[0] / 2)
                    end_ = test.shape[0]
                    li_temp = []
                    mid_ = int(test.shape[0] / 2)
                    end_ = test.shape[0]
                    for i in range(mid_, end_):
                        li = test[i][end_ - 1]
                        if li > 0:
                            li_temp.append(i)
                        i += 1
                    try:
                        li_max.append(min(li_temp))
                        li_max.append(max(li_temp))
                    except:
                        li_max.append(-1)
                        li_max.append(-1)

                    li_temp2 = []
                    for i in range(0, mid_):
                        li = test[i][end_ - 1]
                        if li > 0:
                            li_temp2.append(i)
                        i += 1
                    try:
                        li_max.append(max(li_temp2))
                        li_max.append(min(li_temp2))
                    except:
                        li_max.append(-1)
                        li_max.append(-1)

                elif i == 3:
                    #test = np.flip(test, 1)
                    mid_ = int(test.shape[0] / 2)
                    end_ = test.shape[0]
                    li_temp = []
                    mid_ = int(test.shape[0] / 2)
                    end_ = test.shape[0]
                    for i in range(mid_, end_):
                        li = test[i][0]
                        if li > 0:
                            li_temp.append(i)
                        i += 1
                    try:
                        li_max.append(min(li_temp))
                        li_max.append(max(li_temp))
                    except:
                        li_max.append(-1)
                        li_max.append(-1)

                    li_temp2 = []
                    for i in range(0, mid_):
                        li = test[i][0]
                        if li > 0:
                            li_temp2.append(i)
                        i += 1
                    try:
                        li_max.append(max(li_temp2))
                        li_max.append(min(li_temp2))
                    except:
                        li_max.append(-1)
                        li_max.append(-1)

                else:
                    break

            swath = np.empty([swaths.shape[0], swaths.shape[1]])
            if n == 0:
                # top left/right from center
                for i in range(li_max[0], li_max[1]):
                    swath[swaths.shape[1] - 1][i] = 10
                    i += 1
                # x 2
                for i in range(li_max[3], li_max[2]):
                    swath[swaths.shape[1] - 1][i] = 10
                    i += 1

            elif n == 1:
                # bottom left/right from center
                for i in range((swaths.shape[1] - li_max[5]), (swaths.shape[1] - li_max[4])):
                    swath[0][i] = 10
                    i += 1
                for i in range((swaths.shape[1] - li_max[6]), (swaths.shape[1] - li_max[7])):
                    swath[0][i] = 10
                    i += 1

            elif n == 2:
                # left up/down from center
                for i in range((swaths.shape[1] - li_max[9]), (swaths.shape[1] - li_max[8])):
                    swath[i][0] = 10
                    i += 1
                # x 2
                for i in range((swaths.shape[1] - li_max[10]), (swaths.shape[1] - li_max[11])):
                    swath[i][0] = 10
                    i += 1

            elif (n == 3):
                # right up/down from center
                for i in range((swaths.shape[1] - li_max[13]), (swaths.shape[1] - li_max[12])):
                    swath[i][swaths.shape[1] - 1] = 10
                    i += 1
                # x 2
                for i in range((swaths.shape[1] - li_max[14]), (swaths.shape[1] - li_max[15])):
                    swath[i][swaths.shape[1] - 1] = 10
                    i += 1

            try:
                swath_ = swath_ + swath
            except:
                swath_ = swaths + swath

        return swath_
    
    def get_shade(self, swaths_, **kwargs):
        tic = time.perf_counter()
        
        display(Markdown("##### Begin computing Seismo-Lineament boundary area"))
        
        corners_fill = kwargs.get('corners_fill', None)
        trend = kwargs.get('trend', 'NS')
        
        if trend is ('EW' or 'WE'):
            swaths_ = swaths_.T
        
        try:
            swaths_1 = swaths.fill_swaths(self, swaths_)
            swaths_2 = swaths.fill_swaths(self, swaths_1)
            swaths_3 = swaths.fill_middle(swaths_2)
            
            if corners_fill is not None:
                for i in corners_fill:
                    try:
                        swaths_c = swaths.fill_corners(self, swaths_c, corners_fill=i)
                    except:
                        swaths_c = swaths.fill_corners(self, swaths_3, corners_fill=i)
            else:
                swaths_c = swaths_3
            
            swaths_z = swaths.fill_swaths(self, swaths_c)

        except:
            swaths_z = swaths_
            
        if trend is ('EW' or 'WE'):
            swaths_z = swaths_z.T
        
        toc = time.perf_counter()
        display(Markdown(f"##### Finished assessment of Seismo-Lineament boundary area in {toc - tic:0.2f} seconds"))
        
        return swaths_z

# Visualize the Data

class SLAM_viz:
    
    def __init__(self, elevation):
        self.elevation = elevation
    
    def elevation_map(self, lightDirection, altDeg, 
                      lon, lat, lon_max, lon_min, lat_min, lat_max, 
                      strike, dip, rake, title, **kwargs):
        
        elevation = self.elevation
        
        faults = kwargs.get('faults', [])
        output_size = kwargs.get('output_size', 1)
        swaths_shade = kwargs.get('swaths_shade', None)
        swaths = kwargs.get('swaths', None)

        # Define the projection
        crs = ccrs.PlateCarree()

        # Set tile type for plotting
        tiler = cimgt.Stamen('terrain')
        mercator = tiler.crs

        fig = plt.figure(figsize=(15, 15))
        #fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1, projection=mercator)
        ax.set_extent([lon_max, lon_min, lat_min, lat_max], crs=crs)

        ax.gridlines(draw_labels=True, linewidth=0)
        ax.set_aspect('auto')

        thin = kwargs.get('thin', True)

        if thin == True:
            image_ = elevation / elevation.min()
        else:
            image_ = elevation

        ls = LightSource(azdeg=lightDirection, altdeg=altDeg)
        im = plt.imshow(ls.hillshade(elevation=image_), cmap='Greys',
                        extent=(lon_max, lon_min, lat_min, lat_max), 
                        transform=crs)

        if len(faults) > 0:
            ax.add_geometries(faults.geometry, crs, alpha=0.5,
                              facecolor='none', edgecolor='yellow', linestyle='solid', linewidth=1.75)
        
        middle_road = kwargs.get('middle_road', None)
        
        if swaths_shade is not None:
            plt.contourf(swaths_shade, extent=(lon_max, lon_min, lat_min, lat_max), 
                         levels=[-1, 1], colors=['k'], alpha = 0.35, transform=crs)
            
            if swaths is not None:
                plt.contour(swaths, extent=(lon_max, lon_min, lat_min, lat_max), levels=[-1, 1],
                            colors=['k'], alpha = 0.25, transform=crs)
            
            if middle_road is not None:
                plt.contour(middle_road, extent=(lon_max, lon_min, lat_min, lat_max), levels=[-1, 1],
                            colors=['k'], linewidths=[2.5], alpha=0.65, transform=crs, zorder=10)

            plt.contour(swaths_shade, extent=(lon_max, lon_min, lat_min, lat_max), levels=[-1, 1], 
                        colors=['lightblue'], linewidths=[3], linestyles=['dashed'], transform=crs)
            
        elif swaths_shade is None:
            if swaths is not None:
                plt.contourf(swaths, extent=(lon_max, lon_min, lat_min, lat_max), 
                             levels=[-1, 1], colors=['k'], alpha = 0.35, transform=crs)
                
                plt.contour(swaths, extent=(lon_max, lon_min, lat_min, lat_max), levels=[-1, 1],
                            colors=['lightblue'], linewidths=[3], linestyles=['dashed'], transform=crs)

        output_size = kwargs.get('output_size', None)
        
        if output_size is not None:
            fig = plt.gcf()
            DPI = fig.get_dpi()
            fig.set_size_inches((swaths_shade.shape[1]/float(DPI)) * output_size, (swaths_shade.shape[0]/float(DPI)) * output_size)
            
        aftershocks_mag = kwargs.get('aftershocks_mag', None)
        aftershocks_lat = kwargs.get('aftershocks_lat', None)
        aftershocks_lon = kwargs.get('aftershocks_lon', None)
        size_multilier = kwargs.get('size_multiplier', 5)
        aftershock_sizeDiff = kwargs.get('sizeDiff', True)
        aftershock_size = kwargs.get('asize', 30)
        
        if aftershocks_mag is not None:
            plt.draw()
            if aftershock_sizeDiff is True:
                for i in range(len(aftershocks_mag)):
                    ax.plot(float(aftershocks_lon[i]), float(aftershocks_lat[i]), marker='o', 
                            markersize=float(aftershocks_mag[i]) * size_multilier, color='red', alpha=0.35, 
                            transform=ccrs.PlateCarree())
            elif aftershock_sizeDiff is False:
                for i in range(len(aftershocks_mag)):
                    ax.plot(float(aftershocks_lon[i]), float(aftershocks_lat[i]), marker='o', 
                            markersize=aftershock_size, color='red', alpha=0.35, 
                            transform=ccrs.PlateCarree())
        
        plt.draw()
        bbox = fig.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        width, height = bbox.width, bbox.height
        my_dpi = fig.dpi

        trans = crs._as_mpl_transform(ax)
        x, y = trans.transform_point((float(lon), float(lat)))
        x_ = ((x/my_dpi))/width
        y_ = ((y/my_dpi))/height
        focal_size = kwargs.get('focal_size', 0.75)
        axi = fig.add_axes([(x_ - (focal_size/width)*0.5), (y_ - (focal_size/height)*0.5), (focal_size/width), (focal_size/height)])
        axi.set_xlim((-7.5, 7.5))
        axi.set_ylim((-7.5, 7.5))

        mt1 = [float(strike), float(dip), float(rake)]

        facecolor = kwargs.get('facecolor', 'k')
        bgcolor = kwargs.get('bgcolor', 'w')

        beach1 = beach(mt1, facecolor=facecolor, bgcolor=bgcolor, width=0.005)

        axi.add_collection(beach1)
        buffer = [patheffects.withStroke(linewidth=4, foreground="whitesmoke")]
        axi.set_title(title, fontdict={'fontsize': 24, 'fontweight': 'bold', 'verticalalignment': 'baseline', 'color': 'k'}, 
                      y=0.95, loc='left', path_effects=buffer)

        axi.axis('equal')
        axi.axis('off')
        
        save_fig = kwargs.get('save_fig', None)
        
        if save_fig is not None:
            plt.savefig(str(save_fig), edgecolor='k', bbox_inches='tight')
            
        plt.show()
        
    def elevation_mult(self, lightDirection, altDeg, 
                      lon, lat, lon_max, lon_min, lat_min, lat_max, 
                      strike, dip, rake, title, **kwargs):
        
        elevation = self.elevation
        
        faults = kwargs.get('faults', [])
        output_size = kwargs.get('output_size', 1)
        swaths_shade = kwargs.get('swaths_shade', None)
        swaths = kwargs.get('swaths', None)

        # Define the projection
        crs = ccrs.PlateCarree()

        # Set tile type for plotting
        tiler = cimgt.Stamen('terrain')
        mercator = tiler.crs

        fig = plt.figure(figsize=(15, 15))
        ax = fig.add_subplot(1, 1, 1, projection=mercator)
        ax.set_extent([lon_max, lon_min, lat_min, lat_max], crs=crs)

        ax.gridlines(draw_labels=True, linewidth=0)
        ax.set_aspect('auto')

        thin = kwargs.get('thin', True)

        if thin == True:
            image_ = elevation / elevation.min()
        else:
            image_ = elevation

        ls = LightSource(azdeg=lightDirection, altdeg=altDeg)
        im = plt.imshow(ls.hillshade(elevation=image_), cmap='Greys',
                        extent=(lon_max, lon_min, lat_min, lat_max), 
                        transform=crs)

        if len(faults) > 0:
            ax.add_geometries(faults.geometry, crs, alpha=0.5,
                              facecolor='none', edgecolor='yellow', linestyle='solid', linewidth=1.75)
        
        middle_road = kwargs.get('middle_road', None)
        
        try:
            if (swaths_shade is not None) and (swaths is not None):
                ran_ = max(len(swaths_shade), len(swaths))
            elif (swaths_shade is None) and (swaths is not None):
                ran_ = len(swaths)
            elif (swaths_shade is not None) and (swaths is None):
                ran_ = len(swaths_shade)
            else:
                ran_ = 0
        except:
            ran_ = 0
        
        for i in range(ran_):
            colors_ = kwargs.get('colors_', ['lightblue'] * len(swaths_shade))
            if swaths_shade is not None:
                plt.contourf(swaths_shade[i], extent=(lon_max, lon_min, lat_min, lat_max), 
                             levels=[-1, 1], colors=['k'], alpha = 0.25, transform=crs)

                if swaths is not None:
                    plt.contour(swaths[i], extent=(lon_max, lon_min, lat_min, lat_max), levels=[-1, 1],
                                colors=['k'], alpha = 0.25, transform=crs)

                if middle_road is not None:
                    plt.contour(middle_road[i], extent=(lon_max, lon_min, lat_min, lat_max), levels=[-1, 1],
                                colors=['k'], linewidths=[2.5], alpha=0.65, transform=crs, zorder=10)

                plt.contour(swaths_shade[i], extent=(lon_max, lon_min, lat_min, lat_max), levels=[-1, 1], 
                            colors=colors_[i], linewidths=[3], linestyles=['dashed'], transform=crs)

            elif swaths_shade is None:
                if swaths is not None:
                    plt.contourf(swaths[i], extent=(lon_max, lon_min, lat_min, lat_max), 
                                 levels=[-1, 1], colors=['k'], alpha = 0.25, transform=crs)

                    plt.contour(swaths[i], extent=(lon_max, lon_min, lat_min, lat_max), levels=[-1, 1],
                                colors=colors_[i], linewidths=[3], linestyles=['dashed'], transform=crs)

        output_size = kwargs.get('output_size', None)
        
        if output_size is not None:
            fig = plt.gcf()
            DPI = fig.get_dpi()
            fig.set_size_inches((swaths_shade[0].shape[1]/float(DPI)) * output_size, (swaths_shade[0].shape[0]/float(DPI)) * output_size)
            
        aftershocks_mag = kwargs.get('aftershocks_mag', None)
        aftershocks_lat = kwargs.get('aftershocks_lat', None)
        aftershocks_lon = kwargs.get('aftershocks_lon', None)
        size_multilier = kwargs.get('size_multiplier', 5)
        aftershock_sizeDiff = kwargs.get('sizeDiff', True)
        aftershock_size = kwargs.get('asize', 30)
        
        if aftershocks_mag is not None:
            plt.draw()
            if aftershock_sizeDiff is True:
                for i in range(len(aftershocks_mag)):
                    ax.plot(float(aftershocks_lon[i]), float(aftershocks_lat[i]), marker='o', 
                            markersize=float(aftershocks_mag[i]) * size_multilier, color='red', alpha=0.35, 
                            transform=ccrs.PlateCarree())
            elif aftershock_sizeDiff is False:
                for i in range(len(aftershocks_mag)):
                    ax.plot(float(aftershocks_lon[i]), float(aftershocks_lat[i]), marker='o', 
                            markersize=aftershock_size, color='red', alpha=0.35, 
                            transform=ccrs.PlateCarree())
        
        plt.draw()
        bbox = fig.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        width, height = bbox.width, bbox.height
        my_dpi = fig.dpi
        
        if (len(lat) > 0) and (len(lon) > 0):
            for i in range(len(lat)):
                colors_ = kwargs.get('colors_', ['w'] * len(lat))
                trans = crs._as_mpl_transform(ax)
                x, y = trans.transform_point((float(lon[i]), float(lat[i])))
                x_ = ((x/my_dpi))/width
                y_ = ((y/my_dpi))/height
                focal_size = kwargs.get('focal_size', 0.75)
                axi = fig.add_axes([(x_ - (focal_size/width)*0.5), (y_ - (focal_size/height)*0.5), (focal_size/width), (focal_size/height)])
                axi.set_xlim((-7.5, 7.5))
                axi.set_ylim((-7.5, 7.5))

                mt1 = [float(strike[i]), float(dip[i]), float(rake[i])]

                facecolor = kwargs.get('facecolor', 'k')
                bgcolor = kwargs.get('bgcolor', 'w')

                beach1 = beach(mt1, facecolor=facecolor, bgcolor=bgcolor, width=0.005)

                axi.add_collection(beach1)
                buffer = [patheffects.withStroke(linewidth=4, foreground="k")]
                axi.set_title(title[i], fontdict={'fontsize': 24, 'fontweight': 'bold', 'verticalalignment': 'baseline', 'color': colors_[i]}, 
                              y=0.95, loc='left', path_effects=buffer)

                axi.axis('equal')
                axi.axis('off')
        
        save_fig = kwargs.get('save_fig', None)
        
        if save_fig is not None:
            plt.savefig(str(save_fig), edgecolor='k', bbox_inches='tight')
        
        plt.show()
        
    def physical_map(self, lon, lat, lon_max, lon_min, lat_min, lat_max, 
                     strike, dip, rake, title, **kwargs):
        
        elevation = self.elevation
        
        faults = kwargs.get('faults', [])
        output_size = kwargs.get('output_size', 1)
        tiler_size = kwargs.get('tiler_size', 1)
        
        swaths_shade = kwargs.get('swaths_shade', None)
        swaths = kwargs.get('swaths', None)

        # Define the projection
        crs=ccrs.PlateCarree()

        # Set tile type for plotting
        tiler = cimgt.Stamen('terrain')
        mercator = tiler.crs

        fig = plt.figure(figsize=(15, 15))
        ax = fig.add_subplot(1, 1, 1, projection=mercator)
        ax.set_extent([lon_max, lon_min, lat_min, lat_max], crs=crs)
        
        try:
            ax.add_image(tiler, tiler_size)
        except:
            ax.add_image(tiler, 1)

        ax.gridlines(draw_labels=True, linewidth=0)
        ax.set_aspect('auto')

        if len(faults) > 0:
            ax.add_geometries(faults.geometry, crs, alpha=0.5,
                              facecolor='none', edgecolor='yellow', 
                              linestyle='solid', linewidth=2)
            
        middle_road = kwargs.get('middle_road', None)

        if swaths_shade is not None:
            plt.contourf(swaths_shade, extent=(lon_max, lon_min, lat_min, lat_max), 
                         levels=[-1, 1], colors=['k'], alpha = 0.35, transform=crs)
            
            if swaths is not None:
                plt.contour(swaths, extent=(lon_max, lon_min, lat_min, lat_max), levels=[-1, 1],
                            colors=['k'], alpha = 0.25, transform=crs)
                
            if middle_road is not None:
                plt.contour(middle_road, extent=(lon_max, lon_min, lat_min, lat_max), levels=[-1, 1],
                            colors=['k'], linewidths=[2.5], alpha=0.65, transform=crs, zorder=10)

            plt.contour(swaths_shade, extent=(lon_max, lon_min, lat_min, lat_max), levels=[-1, 1], 
                        colors=['lightblue'], linewidths=[3], linestyles=['dashed'], transform=crs)
            
        elif swaths_shade is None:
            if swaths is not None:
                plt.contourf(swaths, extent=(lon_max, lon_min, lat_min, lat_max), 
                             levels=[-1, 1], colors=['k'], alpha = 0.35, transform=crs)
                
                plt.contour(swaths, extent=(lon_max, lon_min, lat_min, lat_max), levels=[-1, 1],
                            colors=['lightblue'], linewidths=[3], linestyles=['dashed'], transform=crs)

        output_size = kwargs.get('output_size', None)
        
        if output_size is not None:
            fig = plt.gcf()
            DPI = fig.get_dpi()
            fig.set_size_inches((swaths_shade.shape[1]/float(DPI)) * output_size, (swaths_shade.shape[0]/float(DPI)) * output_size)
            
        aftershocks_mag = kwargs.get('aftershocks_mag', None)
        aftershocks_lat = kwargs.get('aftershocks_lat', None)
        aftershocks_lon = kwargs.get('aftershocks_lon', None)
        size_multilier = kwargs.get('size_multiplier', 5)
        aftershock_sizeDiff = kwargs.get('sizeDiff', True)
        aftershock_size = kwargs.get('asize', 30)
        
        if aftershocks_mag is not None:
            plt.draw()
            if aftershock_sizeDiff is True:
                for i in range(len(aftershocks_mag)):
                    ax.plot(float(aftershocks_lon[i]), float(aftershocks_lat[i]), marker='o', 
                            markersize=float(aftershocks_mag[i]) * size_multilier, color='red', alpha=0.35, 
                            transform=ccrs.PlateCarree())
            elif aftershock_sizeDiff is False:
                for i in range(len(aftershocks_mag)):
                    ax.plot(float(aftershocks_lon[i]), float(aftershocks_lat[i]), marker='o', 
                            markersize=aftershock_size, color='red', alpha=0.35, 
                            transform=ccrs.PlateCarree())

        plt.draw()
        bbox = fig.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        width, height = bbox.width, bbox.height
        my_dpi = fig.dpi

        trans = crs._as_mpl_transform(ax)
        x, y = trans.transform_point((float(lon), float(lat)))
        x_ = ((x/my_dpi))/width
        y_ = ((y/my_dpi))/height
        
        focal_size = kwargs.get('focal_size', 0.75)
        axi = fig.add_axes([(x_ - (focal_size/width)*0.5), (y_ - (focal_size/height)*0.5), (focal_size/width), (focal_size/height)])  

        axi.set_xlim((-7.5, 7.5))
        axi.set_ylim((-7.5, 7.5))

        mt1 = [float(strike), float(dip), float(rake)]
        
        facecolor = kwargs.get('facecolor', 'k')
        bgcolor = kwargs.get('bgcolor', 'w')

        beach1 = beach(mt1, facecolor=facecolor, bgcolor=bgcolor, width=0.05)

        axi.add_collection(beach1)

        buffer = [patheffects.withStroke(linewidth=4, foreground="whitesmoke")]
        axi.set_title(title, fontdict={'fontsize': 24, 'fontweight': 'bold', 'verticalalignment': 'baseline', 'color': 'k'}, 
                      y=0.95, loc='left', path_effects=buffer)

        axi.axis('equal')
        axi.axis('off')
        
        save_fig = kwargs.get('save_fig', None)
        
        if save_fig is not None:
            plt.savefig(str(save_fig), edgecolor='k', bbox_inches='tight')

        plt.show()
        display(Markdown("Map tiles by <a href='http://stamen.com'>Stamen Design</a>, under <a href='http://creativecommons.org/licenses/by/3.0'>CC BY 3.0</a>. Data by <a href='http://openstreetmap.org'>OpenStreetMap</a>, under <a href='http://www.openstreetmap.org/copyright'>ODbL</a>."))
        
    def physical_mult(self, lon, lat, lon_max, lon_min, lat_min, lat_max, 
                     strike, dip, rake, title, **kwargs):
        
        elevation = self.elevation
        
        faults = kwargs.get('faults', [])
        output_size = kwargs.get('output_size', 1)
        tiler_size = kwargs.get('tiler_size', 1)
        
        swaths_shade = kwargs.get('swaths_shade', None)
        swaths = kwargs.get('swaths', None)

        # Define the projection
        crs=ccrs.PlateCarree()

        # Set tile type for plotting
        tiler = cimgt.Stamen('terrain')
        mercator = tiler.crs

        fig = plt.figure(figsize=(15, 15))
        ax = fig.add_subplot(1, 1, 1, projection=mercator)
        ax.set_extent([lon_max, lon_min, lat_min, lat_max], crs=crs)
        
        try:
            ax.add_image(tiler, tiler_size)
        except:
            ax.add_image(tiler, 1)

        ax.gridlines(draw_labels=True, linewidth=0)
        ax.set_aspect('auto')

        if len(faults) > 0:
            ax.add_geometries(faults.geometry, crs, alpha=0.5,
                              facecolor='none', edgecolor='yellow', 
                              linestyle='solid', linewidth=2)
            
        middle_road = kwargs.get('middle_road', None)
        
        try:
            if (swaths_shade is not None) and (swaths is not None):
                ran_ = max(len(swaths_shade), len(swaths))
            elif (swaths_shade is None) and (swaths is not None):
                ran_ = len(swaths)
            elif (swaths_shade is not None) and (swaths is None):
                ran_ = len(swaths_shade)
            else:
                ran_ = 0
        except:
            ran_ = 0
        
        for i in range(ran_):
            colors_ = kwargs.get('colors_', ['lightblue'] * len(swaths_shade))
            if swaths_shade is not None:
                plt.contourf(swaths_shade[i], extent=(lon_max, lon_min, lat_min, lat_max), 
                             levels=[-1, 1], colors=['k'], alpha = 0.25, transform=crs)

                if swaths is not None:
                    plt.contour(swaths[i], extent=(lon_max, lon_min, lat_min, lat_max), levels=[-1, 1],
                                colors=['k'], alpha = 0.25, transform=crs)

                if middle_road is not None:
                    plt.contour(middle_road[i], extent=(lon_max, lon_min, lat_min, lat_max), levels=[-1, 1],
                                colors=['k'], linewidths=[2.5], alpha=0.65, transform=crs, zorder=10)

                plt.contour(swaths_shade[i], extent=(lon_max, lon_min, lat_min, lat_max), levels=[-1, 1], 
                            colors=colors_[i], linewidths=[3], linestyles=['dashed'], transform=crs)

            elif swaths_shade is None:
                if swaths is not None:
                    plt.contourf(swaths[i], extent=(lon_max, lon_min, lat_min, lat_max), 
                                 levels=[-1, 1], colors=['k'], alpha = 0.25, transform=crs)

                    plt.contour(swaths[i], extent=(lon_max, lon_min, lat_min, lat_max), levels=[-1, 1],
                                colors=colors_[i], linewidths=[3], linestyles=['dashed'], transform=crs)

        output_size = kwargs.get('output_size', None)
        
        if output_size is not None:
            fig = plt.gcf()
            DPI = fig.get_dpi()
            fig.set_size_inches((swaths_shade[0].shape[1]/float(DPI)) * output_size, (swaths_shade[0].shape[0]/float(DPI)) * output_size)
        
        aftershocks_mag = kwargs.get('aftershocks_mag', None)
        aftershocks_lat = kwargs.get('aftershocks_lat', None)
        aftershocks_lon = kwargs.get('aftershocks_lon', None)
        size_multilier = kwargs.get('size_multiplier', 5)
        aftershock_sizeDiff = kwargs.get('sizeDiff', True)
        aftershock_size = kwargs.get('asize', 30)
        
        if aftershocks_mag is not None:
            plt.draw()
            if aftershock_sizeDiff is True:
                for i in range(len(aftershocks_mag)):
                    ax.plot(float(aftershocks_lon[i]), float(aftershocks_lat[i]), marker='o', 
                            markersize=float(aftershocks_mag[i]) * size_multilier, color='red', alpha=0.35, 
                            transform=ccrs.PlateCarree())
            elif aftershock_sizeDiff is False:
                for i in range(len(aftershocks_mag)):
                    ax.plot(float(aftershocks_lon[i]), float(aftershocks_lat[i]), marker='o', 
                            markersize=aftershock_size, color='red', alpha=0.35, 
                            transform=ccrs.PlateCarree())
        
        plt.draw()
        bbox = fig.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        width, height = bbox.width, bbox.height
        my_dpi = fig.dpi
        
        if (len(lat) > 0) and (len(lon) > 0):
            for i in range(len(lat)):
                trans = crs._as_mpl_transform(ax)
                x, y = trans.transform_point((float(lon[i]), float(lat[i])))
                x_ = ((x/my_dpi))/width
                y_ = ((y/my_dpi))/height

                focal_size = kwargs.get('focal_size', 0.75)
                axi = fig.add_axes([(x_ - (focal_size/width)*0.5), (y_ - (focal_size/height)*0.5), (focal_size/width), (focal_size/height)])  

                axi.set_xlim((-7.5, 7.5))
                axi.set_ylim((-7.5, 7.5))

                mt1 = [float(strike[i]), float(dip[i]), float(rake[i])]

                facecolor = kwargs.get('facecolor', 'k')
                bgcolor = kwargs.get('bgcolor', 'w')

                beach1 = beach(mt1, facecolor=facecolor, bgcolor=bgcolor, width=0.05)

                axi.add_collection(beach1)

                buffer = [patheffects.withStroke(linewidth=4, foreground="whitesmoke")]
                axi.set_title(title[i], fontdict={'fontsize': 24, 'fontweight': 'bold', 'verticalalignment': 'baseline', 'color': 'k'}, 
                              y=0.95, loc='left', path_effects=buffer)

                axi.axis('equal')
                axi.axis('off')
        
        save_fig = kwargs.get('save_fig', None)
        
        if save_fig is not None:
            plt.savefig(str(save_fig), edgecolor='k', bbox_inches='tight')

        plt.show()
        display(Markdown("Map tiles by <a href='http://stamen.com'>Stamen Design</a>, under <a href='http://creativecommons.org/licenses/by/3.0'>CC BY 3.0</a>. Data by <a href='http://openstreetmap.org'>OpenStreetMap</a>, under <a href='http://www.openstreetmap.org/copyright'>ODbL</a>."))