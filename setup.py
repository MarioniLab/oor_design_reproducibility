import setuptools

setuptools.setup(
    name="diff2atlas",
    version="0.0.99",
    author="Emma Dann",
    author_email="ed6@sanger.ac.uk",
    description="Differential analysis on atlases",
    url="https://github.com/emdann/diff2atlas",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    license='MIT',
    packages=['diff2atlas'],
    install_requires=[
        "pandas",
        "anndata",
        "scanpy>=1.8.0",
        "scipy",
        "numpy",
        "matplotlib",
        "sklearn"
    ],
    zip_safe=False,
    python_requires=">=3.6"
)
