from setuptools import setup, find_packages
 
classifiers = [
   'Development Status :: 5 - Production/Stable',
  'Operating System :: Microsoft :: Windows :: Windows 10',
  'Operating System :: POSIX :: Linux',
  'License :: OSI Approved :: MIT License',
  'Programming Language :: Python :: 3'
]
 
setup(
  name='ImputeFISH',
  version='0.0.4',
  description='ImputeFISH: A Spatially Aware and Scalable Framework for Gene Imputation in Imaging-Based Spatial Transcriptomics',
  long_description=open('README.md').read() ,
  url='',  
  author='Asish Kumar Swain',
  author_email='swainkasish@gmail.com',
  license='MIT',
  classifiers=classifiers,
  keywords='willupdate', 
  packages=find_packages(),
  install_requires=["scikit-learn>=1.7.2",
        "scipy>=1.15.3",
        "seaborn>=0.13.2"]
)
