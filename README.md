# Local-Volatility-Model

## Objectives
This project aims to strip the Equity Local Volatility using SABR fitting method and test the calibration quality with Monte-Carlo.

# Setting up the project
Run the file following commands, it will install dependencies and create a virtual environment for the project :
In Windows :
```bat
python -m venv .venv
@REM Could also be python3 -m venv .venv 
.\.venv\Scripts\pip.exe install -r requirements.txt
```
In Linux/MacOS
```bash
python -m venv .venv 
# Could also be python3 -m venv .venv 
.\.venv\bin\pip install -r requirements.txt -U
```

## Project Structure

Here's an overview of the project's structure:

### Description

- **src/**: Contains all the source code for the project.
  - `modele/`: Contains all the modele of the project
- **install_for_windows.bat**: Batch file to set up the project on Windows.
- **README.md**: The readme file.
- **notebook.ipynb**: Contain all the application of the project.
- **requirements.txt**: Contains the list of dependencies for the project.

## License

[MIT](https://choosealicense.com/licenses/mit/)

## Authors

Authors: **Naïm Lehbiben**
