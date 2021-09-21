from setuptools import setup, find_packages


def get_long_description():
    with open("README.md", "r", encoding="utf8") as file:
        long_description = file.read()
    return long_description


setup(
    name='PySDM-examples',
    description='PySDM usage examples reproducing results from literature and depicting how to use PySDM from Python Jupyter notebooks',
    use_scm_version=True,
    setup_requires=['setuptools_scm'],
    install_requires=['PySDM',
                      'PyMPDATA',
                      'atmos-cloud-sim-uj-utils',
                      'pystrict',
                      'joblib',
                      'matplotlib',
                      'ipywidgets',
                      'ghapi'],  # note: test-time-requirement?
    author='https://github.com/orgs/atmos-cloud-sim-uj/people',
    author_email='sylwester.arabas@uj.edu.pl',
    long_description=get_long_description(),
    long_description_content_type="text/markdown",
    url="https://github.com/atmos-cloud-sim-uj/PySDM-examples",
    license="GPL-3.0",
    packages=find_packages(include=['PySDM_examples', 'PySDM_examples.*'])
)
