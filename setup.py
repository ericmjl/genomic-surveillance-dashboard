from setuptools import setup  # Always prefer setuptools over distutils
from pip.req import parse_requirements

install_reqs = parse_requirements('requirements.txt', session=False)
reqs = [str(ir.req) for ir in install_reqs]

print(reqs)

packages = ['gsdash']

setup(name='gsdash',
      version='2016-06',
      author='Eric J. Ma',
      author_email='ericmajinglong@gmail.com',
      description=("Genomic surveillance dashboard components."),
      license="MIT",
      keywords="infectious disease, genomics, machine learning",
      url='https://github.com/ericmjl/genomic-surveillance-dashboard',
      packages=packages,
      maintainer='Eric J. Ma',
      maintainer_email='ericmajinglong@gmail.com',
      install_requires=reqs,
      long_description='The components powering a genomic surveillance ' +
                       'dashboard backend.',
      classifiers=["Topic :: Scientific/Engineering :: Visualization"]
      )
