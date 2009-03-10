from setuptools import setup
def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('st_cov_fun',parent_package,top_path)
    config.add_extension(name='fst_cov_fun',sources=['st_cov_fun/fst_cov_fun.f'])
    config.make_config_py()
    return config

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


