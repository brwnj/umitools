#!/usr/bin/env python
# encoding: utf-8
import ez_setup
ez_setup.use_setuptools()

from setuptools import setup

setup(name='umitools',
      version='0.1.1',
      description='handle reads with an incorporated UMI',
      author='Joe Brown',
      author_email='brwnjm@gmail.com',
      license='MIT',
      url='https://github.com/brwnjm/umitools',
      install_requires=['pysam', 'toolshed'],
      scripts=['umitools'],
      long_description=open('README.md').read(),
      classifiers=['License :: OSI Approved :: MIT License'],
)