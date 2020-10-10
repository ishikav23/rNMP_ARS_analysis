# rNMP_ARS_analysis
Scripts to analysis rNMP incorporation in yeast Autonomously Replication Sequences.

## Dependency

Except Python3 standard libraries, the following packages are needed to run all the scripts:

- Matplotlib
- Seaborn
- Numpy
- Scipy
- Pandas
- Scikit-learn

## Usage

### Calculate rNMP incorporation numbers in ARS region

BED files for all ARS's used are stored in __ARS_bed__ folder. The __get_flanks.py__ script is used to generate ARS flanks. Then, [__RibosePrefereneceAnalysis__](https://github.com/xph9876/RibosePreferenceAnalysis) package is used to count rNMPs inside each ARS region or each ARS flank and generate corresponding background frequencies. You may use __get_region.py__ and __get_bg_region.py__ to select the region you want, use __normalize_ars.py__ for normalization, and use __merge.py__ to merge several normalized frequency files.

### rNMP incorporation rate change simulation

The simulation of rNMP incorporation rate change is performed by __rate_simulation.py__. You may change the parameter inside the scripts if the rNMP incorporation rate for each DNA polymerase is different with the wild-type.

### Plotting

The [__RibosePrefereneceAnalysis__](https://github.com/xph9876/RibosePreferenceAnalysis) package is used to generate the heatmaps. The repository also contains several scripts to generate other figures as following, and you can run scripts with "__--help__" for detailed usage:

1. __draw_bar_plot.py__: Bar charts for rNMP incorporation percentage on the leading or lagging strand. You may use __sort.py__ to change the order.
2. __draw_lela.py__: Scatter plots to compare rNMPs on the leading and lagging strand for each library and bar charts for leading/lagging ratio.
3. __check_time.py__: Scatter plots to discover the relation of rNMP incorporation leading/lagging ratio and corresponding ARS firing time.
4. __draw_ars_split.py__: Line charts for rNMP incorporation rate change and leading/lagging ratio during DNA replication. 
5. __generate_box_plot.py__: Box plots to compare dinucleotide frequency on the leading and lagging strand in a particular range.

## License

This software is under GNU GPL v3.0 license

## Contact

If you have any question, please contact me at [pxu64@gatech.edu](mailto:pxu64@gatech.edu).

