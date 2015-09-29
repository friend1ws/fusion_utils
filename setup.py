#!/usr/bin/env python

from distutils.core import setup

setup(name='fusion_comp',
      version='0.1.0',
      description='Python tools for comparing gene fusion and structural variation results/',
      author='Yuichi Shiraishi',
      author_email='friend1ws@gamil.com',
      url='https://github.com/friend1ws/fusion_comp',
      package_dir = {'': 'lib'},
      packages=['fusion_comp'],
      scripts=['fusion_comp'],
      license='GPL-3'
     )

