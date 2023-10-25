import setuptools

with open('README.md','r') as fh:
    long_description=fh.read()

setuptools.setup(
    name='obsplanning',
    version='1.1.0',
    url='https://github.com/pjcigan/obsplanning', #'http://obsplanning.readthedocs.io',
    license='MIT',
    author='Phil Cigan',
    author_email='',
    description='Various utilities to assist with astronomical observation planning',
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(exclude=("docs","images")),
    py_modules=['obsplanning',],
    #python_requires='>2.7',
    install_requires=['numpy', 'matplotlib', 'astropy', 'datetime', 'pytz', 'timezonefinder', 'ephem', 'astroquery', 'tqdm'], 
    extras_require={'scipy':  ["scipy",], 'multicolorfits': ["multicolorfits",],},
)
