from setuptools import setup, find_packages

setup(
    name='PyGA',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'numpy', 'matplotlib'
    ],
    description='A simple Genetic Algorithm library.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url="https://github.com/Palundy/PyGA",
    author='Palau van Helden',
    author_email='palau@vhelden.nl'
)