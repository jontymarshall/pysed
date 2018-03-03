from setuptools import setup,find_packages

setup(name='pysed',
	version='1.0',
	description='SED fitting or stellar photometry',
	author="Zachary Draper",
	author_email='zackdraper@gmail.com',
	url='',
	packages=find_packages(),
    	classifiers=[
		'Intended Audience :: Science/Research',
		'Topic :: Scientific/Engineering :: Astronomy',

		'License :: OSI Approved :: BSD License',

		'Programming Language :: Python :: 2.7',
		'Programming Language :: Python :: 3.5',
		],
	install_requires=['numpy', 'scipy', 'astropy', 'astroquery', 'matplotlib']
)


