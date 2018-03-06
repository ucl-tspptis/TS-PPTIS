from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(

    name='TS-PPTIS',

    version='1.3.2',

    description='A new python implementation of Transition State - Partial Path Transition Interface Sampling',
    long_description=long_description,

    url='https://github.com/ucl-tspptis/TS-PPTIS',
	download_url = 'https://github.com/ucl-tspptis/TS-PPTIS/archive/2.0.0.tar.gz', 
    author='Federico Comitani, Giulio Mattedi',
    author_email='f.comitani@ucl.ac.uk, g.mattedi.16@ucl.ac.uk',

    license='GPL-3.0',

    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'License :: OSI Approved :: GNU General Public License v3.0',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.3'
		],

    keywords='tspptis, md, kinetics',

    packages=find_packages(exclude=['contrib', 'docs', 'tests','testfiles']),

    install_requires=['numpy','mdtraj'],

    extras_require={},

    entry_points={
        'console_scripts': [
            'tspptis=tspptis:main',
        ],
    },
)

