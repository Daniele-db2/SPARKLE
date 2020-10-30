from setuptools import setup,find_packages


setup(name='SPARKLE',      
  version='1.0.0',
  description='Proof of concept: alignment of FASTQ files using SPARK',
  requires=['python (>= 3.7)'],
  author='Daniele Bellini',
  author_email='daniele.bellini@stud.unifi.it',
  license='LICENSE.txt',
  install_requires=['pyspark >= 3.0.1', 'numpy >= 1.18.5', 'tabulate >= 0.8.7'],
  zip_safe=False,
  packages=find_packages(),
  include_package_data=True,
  entry_points={'console_scripts': ['SPARKLE=SPARKLE.SPARKLE:main']}          
)

