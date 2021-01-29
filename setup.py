import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="SBSC",
    version="0.0.12",
    license='MIT',
    author="Liam McIntyre",
    author_email="shimbalama@gmail.com",
    description="A simple somatic SNV and MNV caller",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/shimbalama/SBSC",
    packages=setuptools.find_packages(),
    scripts = ['bin/SBSCall.py'], # Add executable to PATH
    include_package_data=True, # for manifest (ie use .txt files)
    zip_safe=False,
    keywords='NA',#describe package
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)

#dependency_links
#dependency_links=['http://github.com/user/repo/tarball/master#egg=package-1.0']

# List PyPI dependencies here
#install_requires=[]
