from setuptools import setup
from numpy.distutils.misc_util import Configuration
config = Configuration('st_cov_fun', top_path='.')
config.add_extension(name='fst_cov_fun',sources=['fst_cov_fun.f'])
config.make_config_py()


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(  version="1.0",
            description="Spatiotemporal covariance function",
            license="Creative Commons License",
            requires=['NumPy','PyMC','PyTables','SciPy'],
            long_description="""
            blablabla
            """,
            **(config.todict()))


