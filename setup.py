from setuptools import setup
import os
import sys

_here = os.path.abspath(os.path.dirname(__file__))

#if sys.version_info[0] < 3:
#    with open(os.path.join(_here, 'description.rst')) as f:
#        long_description = f.read()
#else:
#    with open(os.path.join(_here, 'description.rst'), encoding='utf-8') as f:
#        long_description = f.read()

version = {}
with open(os.path.join(_here, 'GEOSLAM', 'version.py')) as f:
    exec(f.read(), version)

setup(
    name='GEOSLAM',
    version=version['__version__'],
    description=('Seismo-Lineament Analysis Method'),
    long_description=('The seismo-lineament analysis method is a tool to spatially correlate a shallow-focus earthquake to the surface trace of the fault that generated it. SLAM is the intellectual property and work product of Vince Cronin. The GEOSLAM Python code was a simple translation from Vince Cronin\'s Mathematica files, where this translation was performed by Luke Pajer.'),
    author='Luke Pajer',
    author_email='luke.pajer@gmail.com',
    url='https://github.com/The-Geology-Guy/GEOSLAM',
    download_url='https://github.com/The-Geology-Guy/GEOSLAM/archive/0.0.2.tar.gz',
    license='MIT',
    packages=['GEOSLAM'],
    install_requires=[
        'pandas',
        'geopandas',
        'utm',
        'math',
        'numpy',
        'itertools',
        'datetime',
        'xmltodict',
        'time',
        'affine',
        'os',
        'requests',
        'io',
        'tqdm',
        'tempfile',
        'zipfile',
        'IPython.display',
        'rasterio',
        'cartopy',
        'matplotlib',
        'obspy',
    ],
    include_package_data=True,
    classifiers=[
        'Intended Audience :: Science/Research',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',],
    )
