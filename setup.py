# Author: Anand Patil
# Date: 6 Feb 2009
# License: Creative Commons BY-NC-SA
####################################

from setuptools import setup
from numpy.distutils.misc_util import Configuration
config = Configuration('nsst_cov_fun', top_path='.')
config.add_extension(name='fnsst_cov_fun',sources=['fnsst_cov_fun.f'])
config.make_config_py()


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(  version="1.0",
            description="Nonstationary spatiotemporal covariance function",
            license="Creative Commons License",
            requires=['NumPy','PyMC','PyTables','SciPy'],
            long_description="""
            blablabla
            """,
            **(config.todict()))


