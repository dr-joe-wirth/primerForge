import sys
from setuptools import setup, find_packages

pythonVersion = (sys.version_info.major, sys.version_info.minor, sys.version_info.micro)

if pythonVersion < (3, 9) or pythonVersion >= (3, 12):
    sys.stderr.write(f"\Incompatible Python version ({'.'.join(map(str, pythonVersion))}). primerForge requires Python >=3.9 and <3.12\n\n")
    sys.exit(1)

setup(
    name='primerforge',
    version='1.3.2',
    author='Joseph S. Wirth',
    packages=find_packages(),
    description='software to identify primers that can be used to distinguish genomes',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    python_requires='>=3.9,<3.12',
    install_requires=[
        'biopython==1.81',
        'khmer>=2.1',
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
)
