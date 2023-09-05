from setuptools import setup, find_packages

setup(
    name='climate_library',
    version='0.1.1',
    packages=find_packages(),
    install_requires=[
        'xarray',
        'xclim'
    ],
    author='Shiv Shankar Singh, Akshit Nanda',
    author_email='shivshankarsingh.py@gmail.com',
    description='A library for calculating climate indices',
    long_description=open('README.md', encoding='utf-8').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/shiv3679/climate-indices',
    license='MIT',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
)
