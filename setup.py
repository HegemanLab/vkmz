import setuptools

setuptools.setup(
    name="vkmz",
    version="1.4.6",
    python_requires=">=3.8",
    description="metabolomics formula prediction and van Krevelen diagram generation",
    author="Mark Esler",
    author_email="eslerm@umn.edu",
    url="https://github.com/HegemanLab/vkmz",
    packages=setuptools.find_packages(),
    entry_points={"console_scripts": ["vkmz = vkmz.__main__:main"]},
    package_data={"vkmz": ["d3.html", "overlay.png", "databases/*"]},
)
