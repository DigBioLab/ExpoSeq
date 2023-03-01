from setuptools import setup, find_packages
setup(
    name = "ExpoSeq",
    version = "0.1.0",
    description = "A pacakge which provides various ways to analyze NGS data from phage display campaigns",
    url = "https://github.com/nilshof01/ExpoSeq",
    author = "Nils Hofmann",
    author_email = "s220672@dtu.dk",
    install_requires = [
                        "numpy",
                        "pandas>=1.5.3",
                        "matplotlib>=3.6.3",
                        "scipy>=1.10.0",
                        "seaborn>=0.12.2",
                        "logomaker",
                        "editdistance",
                        "networkx>=2.6.3",
                        "PyQt5",
                        "sgt>=2.0.3",
                        "scikit-learn>=1.2.1",
                        "openai"
                        ],
    packages = find_packages()

)