from setuptools import setup


def readme():
    with open('README.rst') as f:
        return f.read()

setup(name='egrm',
      version='0.1',
      description='Coalescent Relationship Matrix',
      long_description=readme(),
      classifiers=[
        'Development Status :: 3 - Alpha',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics', 
      ],
      keywords='genetics genome SNP coalescence',
      url='https://github.com/Ephraim-usc/egrm.git',
      author='Caoqi Fan',
      author_email='caoqifan@usc.edu',
      license='USC',
      packages=['egrm'],
      install_requires=[
          'tskit', 'tqdm', 'msprime'
      ],
      scripts=['bin/trees2egrm', 'bin/workflow'],
      zip_safe=False)
