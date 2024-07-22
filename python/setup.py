from setuptools import setup, find_packages

setup(
    name='skyweaver',
    version='0.1.0',
    author='Ewan/Vivek/Colin/Vishnu',
    author_email='ebarr@mpifr-bonn.mpg.de',
    description='A package for offline beamforming for MeerKAT Baseband data',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/erc-compact/skyweaver/tree/python_wrapper_dev',
    packages=find_packages(include=['skyweaver', 'skyweaver.*', 'skyweaver_test_utils', 'skyweaver_test_utils.*']),
    install_requires=[
        'astropy==6.0.1',
        'astropy-iers-data==0.2024.7.8.0.31.19',
        'Cartopy==0.23.0',
        'certifi==2024.7.4',
        'contourpy==1.2.1',
        'cycler==0.12.1',
        'ephem==4.1.5',
        'fonttools==4.53.1',
        'future==1.0.0',
        'geographiclib==2.0',
        'h5py==3.11.0',
        'importlib_resources==6.4.0',
        'katpoint==0.10.1',
        'kiwisolver==1.4.5',
        'markdown-it-py==3.0.0',
        'matplotlib==3.9.1',
        'mdurl==0.1.2',
        'mosaic @ https://github.com/wchenastro/Mosaic/archive/refs/heads/master.zip#sha256=5d2a701fd151b26af26ef4f1acea3fae22dc9ff42a090e48838708a6af059242',
        'numpy==1.26.4',
        'nvector==0.7.7',
        'packaging==24.1',
        'pillow==10.4.0',
        'pyerfa==2.0.1.4',
        'Pygments==2.18.0',
        'pyparsing==3.1.2',
        'pyproj==3.6.1',
        'pyshp==2.3.1',
        'python-dateutil==2.9.0.post0',
        'PyYAML==6.0.1',
        'rich==13.7.1',
        'scipy==1.13.1',
        'shapely==2.0.4',
        'six==1.16.0',
        'typing_extensions==4.12.2'
    ],
    entry_points={
        'console_scripts': [
            'skyweaver=skyweaver.cli:cli',  
            'sw=skyweaver.cli:cli',         # Adding the alias 'sw' 
        ],
    },
    classifiers=[
        'Programming Language :: Python :: 3.9',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='==3.9.*',
)
