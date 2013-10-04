from ez_setup import use_setuptools
use_setuptools()

from setuptools import setup, find_packages

setup(name = "prism",
      version = "0.1.0",
      packages = find_packages(),
      entry_points = {
        'console_scripts': [
          'prism = libprism.prism:main',
          'prism-interleaving = libprism.prism_interleaving:main',
        ]
      },

      install_requires=[
        'numpy>=1.5.1',
        'scipy>=0.8.0',
      ],

      author='Volodymyr Kuleshov',
      description='Prism statistical phaser.'
    )
