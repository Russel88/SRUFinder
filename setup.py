import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="srufinder", 
    version="0.0.1",
    author="Jakob Russel",
    author_email="russel2620@gmail.com",
    description="SRUFinder: Find and subtype SRUs, mini-arrays, and CRISPR arrays by repeat matching",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Russel88/SRUFinder",
    download_url="https://github.com/Russel88/SRUFinder/archive/v0.0.1.tar.gz",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Development Status :: 4 - Beta"],
    python_requires='>=3.8',
    install_requires=[
        "numpy >= 1.17.5",
        "pandas >= 1.0.1",
        "biopython >= 1.76",
        "multiprocess >= 0.70.9",
        "setuptools"],
    scripts=['bin/srufinder']
)
