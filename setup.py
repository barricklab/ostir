import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="osTIR-barricklab", # Replace with your own username
    version="0.0.1",
    author="Example Author",
    author_email="author@example.com",
    description="Open Source Transcription Initiation",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/barricklab/rbs-calculator",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)