#!/usr/bin/env python
# coding: utf-8


import pandas as pd
import json 
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.stats import mstats
from scipy.stats import mannwhitneyu
import glob
import re
import sys
from numpy import nan
import os
import argparse
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
from datetime import date
from datetime import datetime
from plotly.subplots import make_subplots
import plotly.graph_objs as go
import os
import tempfile
os.environ['MPLCONFIGDIR'] = tempfile.mkdtemp()



# defining what arguments the program requires

parser = argparse.ArgumentParser(description='Generate report from stats_data')


parser.add_argument('p1', metavar='path1', type=str, # path to the file where the JSON file with report_general_stats_data is located
                    help='path to the multiqc_data.json')

parser.add_argument('p2', metavar='path2', type=str, # path to the file where the table where the samples are divided into groups (control and experiment) is located
                    help='path to the table (pacient vs. control)')

parser.add_argument('-p3', metavar='path3', type=str, # path to the files containing the mtDNA data, example input (samples/**/genome_results.txt, where samples is the name of the file)
                    help='path to the mitochondrial information data')

args = parser.parse_args()


# setting paths to variables to be able to work with them in the program
PATH1_MULTIQC = args.p1
PATH2_TABLE = args.p2

if args.p3 :
    PATH3_MT = args.p3
else :
  
    PATH3_MT = ''


# working with JSON file, loading and data processing

with open(PATH1_MULTIQC ) as file:  
    raw_json = json.load(file)
    
read_one_file = True  # auxiliary conversion to always open max. one file

group_table = pd.read_excel (PATH2_TABLE) # reading data from an excel spreadsheet

sample_name = group_table.columns[0]

control_group = ''
pacient_group = ''


# dividing the groups according to the data from the table into control and experimental

if (group_table.iloc[0][0] == 'Control' or group_table.iloc[0][0] == 'kontrola') : #
    kontrola_group = group_table.iloc[0][0]
    pacient_group = group_table.iloc[0][1]
elif (group_table.iloc[0][1] == 'Control' or group_table.iloc[0][1] == 'kontrola'):
    control_group = group_table.iloc[0][1]
    pacient_group = group_table.iloc[0][0]
else :
    sys.exit('Modify the table to make it clear what is the control group and what is the non-control group.')
    
 

# saving read data to a list and deleting null values from the list if there are any
group_table = group_table.shift(-1)
control = group_table[group_table.columns[0]].tolist()
pacient = group_table[group_table.columns[1]].tolist()

control = [str(x) for x in control if pd.isnull(x) == False and x != 'nan']
pacient = [str(x) for x in pacient if pd.isnull(x) == False and x != 'nan']



# in the case of a JSON file, in some datasets the report_general_stats_data is located in a different place in the structure, so we load the correct one
i = 2
if len(raw_json['report_general_stats_data']) == 2:
    i = 1


general_stats  = {**raw_json['report_general_stats_data'][0], **raw_json['report_general_stats_data'][i]}


general_stats =  pd.DataFrame(general_stats).T
general_stats.index.name = 'control group'
general_stats.reset_index(level=0, inplace=True)


# dividing into CTRL group and EXPR group
for i in general_stats['control group']:
    res = list(filter(i.startswith, control)) != []
    if res == True:
        general_stats.loc[general_stats['control group'] == i, 'control group'] = 'CTRL_' + str(i)

for i in general_stats['control group']:
    res = list(filter(i.startswith, pacient)) != []
    if res == True:
        general_stats.loc[general_stats['control group'] == i, 'control group'] = 'EXPR_' + str(i)  

# erging with paired reads
for control_group in general_stats['control group'] :
    if control_group.upper().endswith('_R1') or control_group.upper().endswith('_R2'):
        general_stats['control group'] = general_stats['control group'].replace(control_group,control_group[:-3])
        


# deleting the `coli` and `ctrl` metrics from table is they exists 

general_stats = general_stats[~general_stats['control group'].astype(str).str.startswith(('coli', 'ctrl'))]


# functions to merge with paired reads
aggregation_functions = {
    'total_reads': 'first', 
    'mapped_reads': 'first', 
    'percentage_aligned': 'first',
    'general_error_rate': 'first',
    'median_coverage': 'first',
    'median_insert_size': 'first', 
    'avg_gc': 'first', 
    '1_x_pc': 'first',
    '5_x_pc': 'first',
    '10_x_pc': 'first',
    '30_x_pc': 'first', 
    '50_x_pc': 'first', 
    'percent_gc': 'mean', 
    'avg_sequence_length': 'mean', 
    'total_sequences': 'max',
    'percent_duplicates': 'max',
    'percent_fails': 'sum',
}

general_stats = general_stats.groupby(general_stats['control group'], as_index=False).aggregate(aggregation_functions).reindex(columns=general_stats.columns)


# adding new column `percetntage_noaliged` and deleting 'mean_coverage'

general_stats['percetntage_noaliged'] = general_stats.apply(lambda x: (((x['total_reads'] - x['mapped_reads'])/x['total_reads'])*100) , axis=1)


if 'mean_coverage' in general_stats.columns :
    general_stats.drop('mean_coverage', inplace=True, axis=1)


# loading and processing of the third file

find_first_pattern = re.compile('(?<=>>>>>>> Coverage per contig\n\n).*', flags = re.S) # searching in a file for Coverage per contig
find_second_pattern = re.compile('(?P<column>N[C|T|W]_\d+.\d)\s\d+\s\d+\s(?P<dna>\d+.\d+)', flags = re.S) #column == name of the sample, dna is the mitochondrial dna

column_names = []
column_values = []
one_column_value = []

third_path = PATH3_MT


# reading data from the third file, if they are
if (third_path) :
    for one_row, file_name in zip(range(0, len(general_stats)),glob.iglob(third_path, recursive=True)):
      
        sample_file = open(file_name, "r").read()
        first_pattern = find_first_pattern.search(sample_file)
        second_pattern = find_second_pattern.findall(first_pattern.group(0))

        for one_pattern in second_pattern :
            column_names.append(one_pattern[0])
            column_values.append(float(one_pattern[1]))

        while read_one_file == True :
            general_stats = general_stats.reindex(general_stats.columns.tolist() + column_names, axis=1) 
            read_one_file = False

        sum_of_all = sum(column_values)

        for one_value in column_values :
            one_column_value.append(one_value/sum_of_all)


        for one_value in range(0, len(column_names)) :
            general_stats.loc[one_row,[ column_names[one_value]]]= [one_column_value[one_value]]

        column_names = []
        column_values = []
        one_column_value = []
# else :
#     # information message that there was no third file in the input
#     print(' The file containing mitochondrial data was not entered')


# WORKING WITH GRAPHS
# functions to add text to graphs, first the EXPR to CTRL group ratio and then the result of Kruskal's stat test
def get_text(control, colitis) :
    
    ratio = get_ratio(control.median(), colitis.median())
    pvalue = get_pValue(control, colitis)
    
    return ('\nratio = {:0.5f}x, p_value = {:0.10f}'.format(ratio, pvalue))

def get_ratio(A,B) :
    if(A == 0 or B == 0) :
        return 0
    else :
        return (B / A)

def get_pValue(A, B) :
    if (A == A.iloc[0]).all() or  (B == B.iloc[0]).all() :
        return 0
    else :
        H, pval = mstats.kruskalwallis(B.tolist(), A.tolist())
    
    return pval

columns = list(general_stats.columns)[1:]
columns_length  = len(columns)
#print(columns_length)

number_of_rows = columns_length // 3

groups = ['CTRL', 'EXPR']

# setting to have max. 3 charts in one row
fig = make_subplots(rows =columns_length // 3, cols= 3 , 
                    shared_yaxes=False,
                    subplot_titles=(columns))

# function to generate individual graphs
def add_to_figure(column_name, row, col) :
    for group in groups :
        fig.add_trace(go.Box(y= general_stats[general_stats['control group'].str.startswith(group)][column_name], 
                             name = group), row = row, col = col)


columns_gen = (c for c in columns)

# generation of individual graphs
for one_row in range(1,( columns_length // 3) + 1) :
    for one_col in range(1, 4) :
        column = next(columns_gen)
        add_to_figure(column, one_row, one_col)


# setting height and width for charts (all graphs, not for individual ones)
fig.update_layout(height= number_of_rows * 350, width=1920)
                
# deleting the legend, since all the information is already added in the graph                 
fig.update_layout(showlegend=False,
                 margin=dict( b=2))

names_dict = {}


# adding text to graphs
def add_value(column_name) :
    names_dict[column_name] = get_text(general_stats[general_stats['control group'].str.startswith('CTRL')][column_name], 
                             general_stats[general_stats['control group'].str.startswith('EXPR')][column_name])

for column_name in columns :
    add_value(column_name)

fig.for_each_annotation(lambda a: a.update(text = f'<b>{a.text}</b>' + '<br>' + names_dict[a.text]))
    
# saving graphs to HTML, which is then  used in HTML part
graphs_html = fig.to_html("file.html")


# adding a table with statistical test results for individual measured characteristics

statistika = pd.DataFrame(columns=['metric','statistic','pvalue'])


for group in general_stats.columns[1:]:
    if(general_stats[group].sum() != 0):
        stat, p= mannwhitneyu(general_stats[general_stats['control group'].str.startswith('CTRL')][group],
                              general_stats[general_stats['control group'].str.startswith('EXPR')][group])
        tmp = {'metric': group ,'statistic': stat, 'pvalue': p, }

        statistika = statistika.append(tmp, ignore_index = True)
      
    
general_stats[general_stats['control group'].str.startswith('CTRL')][group]

statistika = statistika.sort_values(
     by="pvalue",
     ascending=False
)


# deleting the prefixes added above, which helped us to identify individual groups
for control_group in general_stats['control group'] :
    if control_group.upper().startswith('EXPR_'):
        general_stats['control group'] = general_stats['control group'].replace(control_group,control_group[5:])
        
for control_group in general_stats['control group'] :
    if control_group.upper().startswith('CTRL_'):
        general_stats['control group'] = general_stats['control group'].replace(control_group,control_group[5:])
    

# auxiliary variable for storing the date and time for the result report
today_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")



# HTML, JS and CSS parts for the result report
html_string = '''
<!doctype html>
<html>
    <head>
        <meta charset="UTF-8">
        <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.11.5/css/jquery.dataTables.css">
        
        <! -- the libraries we work with --> 

        <script src="https://code.jquery.com/jquery-3.6.0.min.js" integrity="sha256-/xUj+3OJU5yExlq6GSYGSHk7tPXikynS7ogEvDej/m4=" crossorigin="anonymous"></script>
        <script type="text/javascript" charset="utf8" src="https://cdn.datatables.net/1.11.5/js/jquery.dataTables.js"></script>
        <title>Report</title>
        

        <! -- adding styles for individual components --> 
        <style>
            body{
                margin:0 100; 
                background:#ffffff;
                font-family:"Open Sans";
            }
        
            select{
                //position: relative;
                display: block;
                //width: 140%;
                margin: 0 auto;
                font-family: "Open Sans";
                font-size: 18px;
                color: #404040;
                border: 1px solid #404040;
                scroll-behavior: smooth;
                overflow: hidden;
                border-radius: 3px;
                border-color: #989898;
                height : 35px;
                margin-bottom: 25px;
            }
        
        </style>
        
    </head>
        <body>
            <! -- importing DataTables styles--> 
            <style> 
                @import url("https://cdn.datatables.net/1.11.5/css/jquery.dataTables.css")
            </style>        
            

            <! -- adding title to the report and obtating the name and today's date and time --> 
            <h1 style="color:#404040;
                    text-align:left;
                    font-family: "Open Sans";">
                Report from ''' + sample_name +''' data - time of generation ''' + today_date + ''' 
            </h1>
        

            <! -- adding the graphs section --> 
            <hr style="height:1.5px;border-width:0;color:#404040;background-color:#404040">
            <h2 style="color:#404040;
                text-align:left;
                font-family: "Open Sans";">
            Graphs generated using the obtained data
            </h2>
            
            <! -- adding select --> 
            <div class="custom-select" style="width:200px;"><select id="graphControls" > </select></div>
            
            <! -- adding graphs to the HTML, the ones we saved eariler graphs_html --> 
            <div id="graphs">''' + graphs_html + ''' </div>
            
            
            <hr style="height:1.5px;border-width:0;color:#404040;background-color:#404040">

            <! -- adding table with all measured characteristics --> 
            <h2 style="color:#404040;
                text-align:left;
                font-family: "Open Sans";">
            Table of all the fragments and their measured characteristics
            </h2>

            ''' + general_stats.to_html(table_id='general_stats') + '''
         
        
            <hr style="height:1.5px;border-width:0;color:#404040;background-color:#404040">
            <! -- adding table with general statistics --> 
            <h2 style="color:#404040;
                text-align:left;
                font-family: "Open Sans";">
            Table providing information on general statistics
            </h2>
         
            ''' + statistika.to_html(table_id='statistics') + '''
         
            <! -- script for beign able to select graphs--> 
            <script> 

                <!-- adding the options for the select//-->
                $(document).ready( function () {
                    $('table').DataTable();
                
                    $("<option/>", {
                            text: "show all", 
                            value: "all"
                    }).appendTo("#graphControls")
                        
                <!-- obtaining the graphs names for the select options//-->     
                $(".annotation-text").each(function (index){
                    const annotation = $(this).children("tspan:first-of-type"); 
                    const text = annotation.text();
                    const id = `graph-annotation-${index}`;
                    annotation.attr("id", id)
                    $("<option/>", {
                        text,
                        value: index
                    }).appendTo("#graphControls")
                });
                    
                    
                <!-- getting the selected option and displaying it //-->  
                $("#graphControls").on("change",function (){
                    const value = $(this).val();
                    const oldStyle =  $("#graphs svg").attr("style")
                        
                    if(value === "all"){
                        $("#graphs").attr("style", "")
                        //$("#graphs svg").attr("style", oldStyle)
                        $("#graphs svg").attr("style", `transform: translate(0px, 0px);`)
                    }
                    else{
                        const index = Number(value);
                        const graphClass = index > 0 ? `.x${index + 1}y${index + 1}` : ".xy";
                        const graph = $(`g${graphClass}`);
                        const rect =  graph.children("rect").first();
                        const x = Number(rect.attr("x")) - 40;
                        const y = Number(rect.attr("y")) - 40;
                        const width = Number(rect.attr("width")) + 125;
                        const height = Number(rect.attr("height")) + 60;
                            
                        $("#graphs").attr("style",  `overflow:hidden; width: ${width}px; height: ${height}px; margin:0 auto;`)
                        $("#graphs svg").attr("style", `transform: translate(-${x}px, -${y}px);`);
                    }
                    })
                } );
            </script>
    
            <!-- adding scrolling to the tables //-->  
            <script>
                $('#general_stats').DataTable( {
                    "scrollY": "650px",
                    "sScrollX": "100%",
                    "scrollCollapse": true,
                    //"paging":         true,
                } );
            </script>

    </body>
</html>'''

# creating the name for the final report
report_name_string = sample_name + str('_report') + str('_'+ date.today().strftime("%b-%d-%Y")+ str('.html'))


# saving 'html' string into HTML format
with open(report_name_string, 'w', encoding = 'utf8') as f:
    f.write(html_string)
f.close()