# NucleoSeqComparator 

It is a program written in Python, that allows the user to compare and align two sequences, plot them on a graph and save alignment data in a text file.

Libraries used: 
 - tkinter: filedialog
 - matplotlib: plotting plots
 - numpy
 - PySimpleGUI: GUI
 - Biopython: obtaining sequences from NCBI GENBANK gb

Available algorythms:

- DotPlot
- Needlemann-Wunsch
- Smith-Watermann

Supported input types:

 - manual text input: user manually inputs sequences into specified textfields
 - input from FASTA file: user selects desired FASTA file via the file explorer and then specifies which sequences to compare
 - GENBANK database input via sequenceID: user inputs sequence IDs from the NCBI GENBANK database, the program downloads sequences and compares them using desired algorythm.

# Screenshots 

![image](https://user-images.githubusercontent.com/61741336/123684117-5f6d1d80-d84d-11eb-9e21-364c537ee424.png)

Comparing two manually inputted sequences using the DotPlot algorythm

![image](https://user-images.githubusercontent.com/61741336/123684258-8deaf880-d84d-11eb-9b5e-56ea831aed64.png)

A comparison of two, rather short NCBI GENBANK sequences, using the Smith-Watermann algorythm. On the side You can see the textfile generated via the save option.

![image](https://user-images.githubusercontent.com/61741336/123684456-c4c10e80-d84d-11eb-9ce1-f18083fbac4d.png)

A comparison of two longer (approx. 500 characters) sequences obtained from NCBI GENBANK db, utilising Smith-Watermann algorythm.
