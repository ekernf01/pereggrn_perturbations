from setuptools import setup
from setuptools import find_packages

with open('README.md', 'r', encoding='utf-8') as fh:
	long_description = fh.read()

setup(
    name='pereggrn_perturbations',	
	py_modules=['pereggrn_perturbations'],
    version='0.0.1',
    description='Efficiently load and validate transcriptome profiles of genetic perturbations',
    long_description=long_description,
	long_description_content_type='text/markdown',
    #url
    author='Eric Kernfeld',
    author_email='eric.kern13@gmail.com',
    install_requires=[
        "pandas",
        "scanpy",
        "anndata",
    ],
    python_requires=">=3.7", 
    url='https://github.com/ekernf01/pereggrn_perturbations',
)