# Local-Volatility-Model

## Objectives
This project aims to strip the Equity Local Volatility using SABR fitting method and test the calibration quality with Monte-Carlo. Arbitrages cleaning on market prices is required since the SABR does not handle them.

## Setting up the project

Run the file `install_for_windows.bat`, it will install dependencies and create a virtual environment for the project.

## Project Structure

Here's an overview of the project's structure:

### Description

- **data/**: Contains datasets and intermediate data files.
- **src/**: Contains all the source code for the project.
  - `data/`: .
  - `utils/`: Additional utility scripts and GUI components.
- **static/**: Contains static files such as the research paper and images.
- **install_for_windows.bat**: Batch file to set up the project on Windows.
- **README.md**: The readme file.
- **requirements.txt**: Contains the list of dependencies for the project.

## License

[MIT](https://choosealicense.com/licenses/mit/)

## Authors

Authors: **Naïm Lehbiben**
