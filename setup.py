"""
PTM-X - a random forest classifier to predict PTM crosstalk intra and inter
proteins.
Web server: http://bioinfo.bjmu.edu.cn/ptm-x/
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path
import PTMXtalk

here = path.abspath(path.dirname(__file__))

# Get the long description from the relevant file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()
    
reqs = ['numpy>=1.12.0', 'scipy>=0.18.1', 'scikit-learn>=0.17', 'joblib>=0.11']

setup(
    name='PTMXtalk',

    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version=PTMXtalk.__version__,

    description='PTMXtalk: PTM cross-talk predictor',
    long_description=long_description,

    # The project's main homepage.
    url='https://github.com/huangyh09/PTM-X',

    # Author details
    author='Yuanhua Huang',
    author_email='yuanhua@ebi.ac.uk',

    # Choose your license
    license='Apache License 2.0',

    # What does your project relate to?
    keywords=['Post-translational modification', 'cross-talk', 
              'co-evolution', 'random forest'],

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages=find_packages(),

    entry_points={
          'console_scripts': [
              'PTMX-predict = PTMXtalk.predict:main',
              'PTMX-feature = PTMXtalk.feature_extractor:main'
              ],
          }, 

    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    
    install_requires=reqs,

    py_modules = ['PTMXtalk']

    # buid the distribution: python setup.py sdist
    # upload to pypi: twine upload dist/...

)