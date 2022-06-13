# Liquidpairs
is a tool that, based on input parameters representing the paths to the files being analyzed, generates a report with the results in which it is possible to observe the differences between the measured characteristics from graphs that are supported by various statistics. The results report also contains the final preprocessed data from which the visualization was created, as well as a table with the results of the statistical test for each measured characteristic. As the final report includes tables with different functions, such as sorting values in columns or searching for certain fragments, it is possible to analyse the fragments separately and there is no need to search for them in complex tables or graphs.

## Program requirements
The following applications are required to run the program successfully: 
-  Python 3 - The program has been tested with version 3.10.4. It is not guaranteed that lower versions will work correctly. The latest version can be downloaded from : https://www.python.org/downloads/
- The final version of the program was tested in Visual Studio Code version 1.67.1, however the program is meant to be executed in other editors as well. Visual Studio Code can be obtained from this link https://code.visualstudio.com/download


## Installation
[Pandas](https://pandas.pydata.org/docs/getting_started/index.html#getting-started)  - library that provides convenient work with data by the means of its data structures. It can be installed using pip : 
```
$ pip install pandas
```
[Plotly](https://plotly.com/python/getting-started/) - library for plotting interactive graphs. It can be installed using pip : 

```
$ pip install plotly==5.8.0
```
[SciPy](https://scipy.github.io/devdocs/getting_started.html#getting-started-ref) - library with a big collection of mathematical functions. It can be installed using pip :

```
$ pip install scipy
```

[NumPy](https://numpy.org/install/) - library used to perform a wide variety of mathematical operations

```
$ pip install numpy
```

## Inputs
The only input to program is a triple of entries that represent paths to given files. Two of these are mandatory files, the third file is optional as is not relevant for certain types of analysis. The first input is the path to the JSON file, which contains the `report_general_stats_data` structure, containing the measured characteristics of the DNA fragments as well as information about their pair reads. The second input is the path to the file, i.e. a table containing the divisions of the control and experimental groups (excel table). An optional path is the path to the files where the information about mitochondrial data is stored in the txt format. 

It is possible to find out what files the program requires by using the `python liquidpairs.py -h ` command, which prints the following information that is described below.

```
positional arguments:
  path1       path to the multiqc_data.json
  path2       path to the table (pacient vs. control)

options:
  -h, --help  show this help message and exit
  -p3 path3   path to the mitochondrial information data
```

## Structure of mandatory and optional files

- First mandatory JSON file : It is important that the input JSON has a structure that contains `report_general_stats_data`. Usually these files are named `multiqc_data.json`, and come as one of the results of the analysis from the MultiQC tool. The structure is described below.

```
{
    {
    ...
    },
    "report_general_stats_data": [
        {
            "DSS3_33_N": {
                "total_reads": 23855690.0,
                "mapped_reads": 23578267.0,
                ....
            },
            ...
}
```

- Second mandatory file, Excel table : here it is important that the table has a precisely defined structure, as described in the example below. It is crucial that it contains a name, which represents the name of the dataset (e.g. where the samples come from), and then two columns, where there is a control and experiments and the samples corresponding to them. *It is very important that at least one cell representing groups is called Control!*
```
    Plasma	                <-------  name
Control	    Tumor           <------- two columns, control and experiment group 
S724	       S531         <------- values 
S726	       S537
S737	       S550
...            ...
```
- Optional input, file containing multiple txt files : here it is important that the folder contains multiple files that contain `genome_results.txt` files. Since the input is a set of files, the path to the input must be specified as `....folder_name/**/genome_results.txt`, since the program works by searching for and working with the given .txt file in each folder. Exaple of structure of folders and the input:

```
input :
.../samples/**/genome_results.txt 
 
where structure : 
 samples
    |_DSS_32_N
        ...
        |_genome_results.txt
        ....
    |_CTRL_12_N
        |_ ...
```

## Output
The output of the program is the result report in HTML format. 

