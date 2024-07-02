from setuptools import setup, find_packages

setup(
    name='primerforge',
    version='1.2.5',
    author='Joseph S. Wirth',
    packages=find_packages(),
    description='software to identify primers that can be used to distinguish genomes',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    install_requires=[
        'biopython==1.81',
        'khmer>=2.1',
        'numpy',
        'primer3-py>=2.0',
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
