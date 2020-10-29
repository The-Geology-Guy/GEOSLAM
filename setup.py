from setuptools import setup
import os
import sys

_here = os.path.abspath(os.path.dirname(__file__))

if sys.version_info[0] < 3:
    with open(os.path.join(_here, 'README.md')) as f:
        long_description = f.read()
else:
    with open(os.path.join(_here, 'README.md'), encoding='utf-8') as f:
        long_description = f.read()

version = {}
with open(os.path.join(_here, 'SLAM', 'version.py')) as f:
    exec(f.read(), version)

setup(
    name='GEOSLAM',
    version=version['0.0.1'],
    description=('Seismo-Lineament Analysis Method'),
    long_description=long_description,
    author='Luke Pajer',
    author_email='luke.pajer@gmail.com',
    url='https://github.com/The-Geology-Guy/GEOSLAM',
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
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',],
    )
