import sys
from Cython.Build import cythonize
from Cython.Compiler import Options
from setuptools import Extension, find_packages, setup

# Enable aggressive optimizations
Options.fast_fail = True

ext_modules=cythonize(
    Extension(
        name="bin.kmer_counting._kmer_counter",
        sources=["bin/kmer_counting/_kmer_counter.pyx"]
    ),
    compiler_directives={
                         'language_level': 3,
                         'boundscheck': False,
                         'wraparound': False,
                         'cdivision': True,
                         'nonecheck': False,
                         'overflowcheck':False
                        }
)


pythonVersion = (sys.version_info.major, sys.version_info.minor, sys.version_info.micro)

setup(
    name='primerforge',
    version='2.1.0',
    author='Joseph S. Wirth',
    packages=find_packages(),
    description='software to identify primers that can be used to distinguish genomes',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    python_requires='>=3.9',
    install_requires=[
        'cython'
        'biopython==1.81',
        'numpy',
        'primer3-py>=2.0',
        'pyahocorasick>=2.1',
        'scipy>=1.10'
    ],
    entry_points={
        'console_scripts': [
            'primerForge=bin.main:main',
        ],
    },
    classifiers=[
        'Programming Language :: Python :: 3.11',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    ext_modules=ext_modules,
    zip_safe=False
)
