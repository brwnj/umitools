#!/usr/bin/env python
# encoding: utf-8

from setuptools import setup

_locals = {}
with open('umitools/_version.py') as fp:
    exec(fp.read(), None, _locals)
version = _locals['__version__']

setup(name='umitools',
    version=version,
    description='Handle reads with an incorporated UMI',
    author='Joe Brown, Jay Hesselberth',
    author_email='brwnjm@gmail.com',
    license='MIT',
    keywords = "bioinformatics",
    url='https://github.com/brwnjm/umitools',
    entry_points={
        'console_scripts': [
            'umitools = umitools.umitools:main',
        ]
    },
    packages=['umitools'],
    install_requires=['pysam', 'toolshed'],
    zip_safe=False,
    long_description=open('README.md').read(),
    classifiers=[
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Utilities'
    ],
)
