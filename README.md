#Final-project-coding

My lab performs many qRT-PCR experiments to determine the effect of different drug treatments on the expression of different genes. Since the lab works on different types of breast cancers, different treatments are used, and the expression of different genes is of interest. I wanted to write a function that can take the raw data from the qRT-PCR machine, containing the well number and the threshold cycle number (which indicates the gene expression), and a template file which stores the name of the gene of interest, the treatments and the timecourses, and generate a graph showing the fold change in gene expression compared to the vehicle treatment. Since treatments and timecourses can vary depending on who sets up the plate, it was important for the template file to be able to be modified by each person for each experiment.

The format for running the function is:

together_timepoints(datafile, templatefile, plottitle, graphname)

###Description of the function inputs

datafile is the .csv file containing the well number and threshold cycle number.

templatefile is the .csv file containing the template of the plate. The gene name, treatment, timepoint and whether it is a housekeepign gene should be separated by a "-".

plottitle is where you can add the name of the graph which will appear on the graph.

graphname is what the filename that the image of the graph will be saved as. It will be saved into your working directory.

