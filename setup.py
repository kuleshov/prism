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

      # data_files=[('', [
      #   'loops'
      # ])],

      data_files=[('libprism/loops', [
        'libprism/loops/viterbi_factorial.c',
        'libprism/loops/forwards_factorial.c',
        'libprism/loops/backwards_factorial.c',
        'libprism/loops/joint_probability.c',
        'libprism/loops/indices.h',
        'libprism/loops/logaddexp.h',
      ])],

      author='Volodymyr Kuleshov',
      description='Prism statistical phaser.'
    )
