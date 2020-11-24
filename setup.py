import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="osTIR-barricklab", # Replace with your own username
    version="0.0.1",
    author="Barrick Lab",
    author_email="croots@utexas.edu",
    description="Open Source Transcription Initiation",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/barricklab/rbs-calculator",
    packages=setuptools.find_packages(exclude="tests"),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    include_package_data=True,
    entry_points='''
        [console_scripts]
        ostir=ostir.ostir:main
        '''
)