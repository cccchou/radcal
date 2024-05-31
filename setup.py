from setuptools import setup, find_packages
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()
setup(
    name='resages',
    version='0.21',
    py_modules=['resages'],
    packages=find_packages(include=['resages', 'resages.*']),
    
    install_requires=[
        
         'numpy',
         'pandas',
         'scipy',
         'tqdm',
    ],
    python_requires='>=3.10',
    package_data={
        '': ['Inputdat/*'],  
    },
    include_package_data=True, 
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    long_description=long_description,
    long_description_content_type="text/markdown",
    author='Alex_chou',
    author_email='alex_chou@sjtu.edu.com',
    description='Methods and codes for reservoir atmosphere 14C age offset calculations',
    url='https://github.com/cccchou/radcal',
)

