#!/usr/bin/env python3
# encoding: utf-8

from setuptools import setup

with open('README.rst') as file:
    readme = file.read()

setup(
    name='dbp',
    version='0.0.0',
    author='Kale Kundert',
    long_description=readme,
    packages=[
        'dbp',
    ],
    install_requires=[
    ],
    entry_points={
        'console_scripts': [
            'dbp_relax_b = bind_dna.scripts.dbp_relax_b:main',
        ],
    },
    include_package_data=True,
)
