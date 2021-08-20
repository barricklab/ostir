import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="OSTIR",
    version="1.0.6",
    author="Cameron Roots, Jeffrey Barrick, Alexandra Lukasiewicz",
    author_email="croots@utexas.edu",
    description="Open Source Transcription Initiation Rates",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/barricklab/ostir",
    packages=setuptools.find_packages(exclude=['tests', 'calibration']),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.8',
    include_package_data=True,
    test_suite="tests",
    entry_points={
        'console_scripts' : [
          'ostir = ostir.ostir:main',
        ],
    }

)
