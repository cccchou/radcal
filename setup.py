from setuptools import setup, find_packages

setup(
    name='resages',
    version='0.2',
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
    
    author='Alex_chou',
    author_email='alex_chou@sjtu.edu.com',
    description='Methods and codes for reservoir atmosphere 14C age offset calculations',
    url='https://github.com/cccchou/radcal',
)

