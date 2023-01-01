from setuptools import setup

setup(
    name = "ExpoSeq",
    version = "0.1.0",
    description = "A pacakge which provides various ways to analyze phage display data from NGS",
    url = "https://github.com/nilshof01/ExpoSeq",
    author = "Nils Hofmann",
    author_email = "s220672@dtu.dk",
    install_requires = ["tk",
                        "numpy",
                        "pandas",
                        "bio",
                        "matplotlib",
                        "scipy",
                        "seaborn",
                        "logomaker",
                        "editdistance",
                        "networkx",
                        "PyQt5"
                        ],

)