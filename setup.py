#!/usr/bin/env python
# encoding: utf-8

from setuptools import setup

_locals = {}
with open('umitools/_version.py') as fp:
    exec(fp.read(), None, _locals)
version = _locals['__version__']

entry_points = """
[console_scripts]
umitools = umitools.umitools:main
"""

setup(name='umitools',
    version=version,
    description='Handle reads with an incorporated UMI',
    author='Joe Brown, Jay Hesselberth',
    author_email='brwnjm@gmail.com',
    license='MIT',
    keywords = "bioinformatics",
    url='https://github.com/brwnjm/umitools',
    entry_points=entry_points,
    packages=['umitools'],
    install_requires=['editdistance', 'pysam'],
    zip_safe=False,
    long_description=open('README.md').read(),
    classifiers=[
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Utilities'
    ],
)
