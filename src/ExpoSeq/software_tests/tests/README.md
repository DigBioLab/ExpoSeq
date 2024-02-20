# How to test in ExpoSeq

I use pdm for organizing the package. If you fork or clone the repository all necessary files will be included. Otherwise you need to install it with pip

```bash
pip install pdm
```

If you have setup pdm and build the infrastructure with pdm.lock and pyproject.toml files you can easiliy run the tests using pdm.
I included a command for using **pytest** in the .toml file, so you just need to type:

```bash
pdm test
```

It will automatically run all tests in all sub-directories which start with test. If you want to create a new test file I would suggest to use the folder test_new where you can call the test files with: 

```bash 
pdm test_new
```

## Test structure in ExpoSeq

The files follow the same structure which the package has, so plot files are under plots and classes for uploading data or change data are most likely tested under settings.
