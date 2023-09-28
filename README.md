# DiArt
![FigureS6_GUIv2](https://github.com/guomics-lab/DIArt/assets/51101922/585495ab-cb2d-4b93-9722-e68e768d36aa)


*indicated the required parameters for running
# Input:
  ○ "mzML dir" button*: Click to select the folder containing the mzML files to be analyzed. Only one folder can be selected, so all mzML files for analysis should be in the same folder.  
  ○ "Clear mzML" button*: Clears the selected mzML folder.  
  ○ "FASTA" button*: Select the FASTA library to be used for analysis.  
  ○ "Lib" button*: Select the Lib library to be used for analysis, in tsv format.param  
# Param
  ○ "IRT":  
    ■ Seed: Random seed for decoy generation.  
    ■ Unit of m/z: Unit for m/z range and tolerance settings, please choose from ('Da', 'ppm').  
    ■ Range of m/z: Minimum and maximum m/z values.  
    ■ Tolerance of m/z: m/z tolerance for MS1 and MS2.  
  ○ "PSMs":  
    ■ Number of PSMs for each precursor to extract  
    ■ Number of fragment ions required for each precursor.  
  ○ "RT cycles":  
    ■ PSMs: Number of RT cycles to search for RSMs for each precursor. Must be greater than (topk + model_cycles - 1).  
    ■ XICs: Number of RT cycles for XICs.  
# Output:  
  ○ "Main output" button: Select the file path for the final result output.  
  ○ Log file: Default location for the log file.  
# Run job:  
  ○ Select the steps to be executed. Default selection is usually sufficient to obtain the final results.  
# Run config:  
  ○ Device: Select the device used for model training. If GPU is available, you can choose "cuda," otherwise, the default is "CPU."  
  ○ #Thread: Number of threads used by DIArt when running DIANN. Choose based on your server's actual conditions.  
  ○ #Epochs: Number of epochs for training the DIArt model. In theory, more training epochs yield better results.  
  ○ Valid prop: Proportion of data used as a validation set during DIArt model training (default is 0.2).  
  ○ Train split num: Number of files the DIArt model is split into during training.
  ○ Batch size: Amount of data loaded in a single batch during DIArt model training (depends on available memory).  
  ○ Task count: Number of subtasks created by DIArt. Once the parameters are set, click the "Run" button to start the process. Clicking the "Stop" button will stop all processes and close the page. The running progress will be displayed in the top-right corner.  
# Log info:  
  ○ The log area currently displays only key step information and progress logs of DIArt. For detailed logs, please refer to the "./diart.log" log file.  

# Demo dataset
The demonstration dataset is accessible through a release in this repository at https://github.com/guomics-lab/DIArt/releases/tag/demodata_v1.0.0. It includes a FASTA library file, a spectral library file, and an mzML file for showcasing the software's functionality.
