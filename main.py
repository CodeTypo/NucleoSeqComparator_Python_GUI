#tkinter, especially filedialog gives us the ability to open up a dialog window and select a file we want to read/modify/parse/etc...
import tkinter
from tkinter import filedialog
from tkinter.filedialog import asksaveasfilename
#matplotlib for plotting 
import matplotlib.pyplot as plt
import numpy as np
#argparse for handling --help argument
import argparse
#Gui 
import PySimpleGUI as sg
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
#Biopython
import os
from Bio import SeqIO
from Bio import Entrez


#A note displaying when the user types in 'python main.py --help'
parser = argparse.ArgumentParser(description='A program made for the purpose of comparing two user-selected sequences of DNA / Aminoacids using Needleman - Wunsch algorythm.'+
' This version supports manual, FASTA fie input modes and NCBI online nucleotide database modes. It comes with a well - described GUI that helps the user to use the application properly.'+
'The program is easy to use and guides the user through the whole time of its execution. Thats why I decided to keep the --help response rather short.')
args = parser.parse_args()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ methods ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#                               FASTA to dict Converter
def fasta2dict(fil):
        """
        THIS CODE FRAGMENT WAS INSPIRED BY: https://gist.github.com/jacob-ogre/5318981 !!
        Parse through the file and look for fasta headers,
        when FASTA header gets encountered, the data structure:
        --->  (key) header : (value) sequence <---
        gets stored in the dictionary, which is later being returned
        as the result of this function

        """
        dic = {}
        current_head = ''
        current_sequence = []
        for line in open(fil):                                  # Parsing through lines in file
            if line.startswith(">") and current_head == '':     # '>' is an indicator of a new entry in the FASTA file
                current_head = line.split(' ')[0]               # The first word after the '>' mark gets stored in the dict as the key of the new key-value pair
            elif line.startswith(">") and current_head != '':   # When the program encounters another entry,
                dic[current_head] = ''.join(current_sequence)   # the previous pair-value pair gets stored in the dict 
                current_head = line.split(' ')[0]               # And the new dict key is being read
                current_sequence = []
            else:
                current_sequence.append(line.rstrip())          #If the line doesn't contain a '>' mark, it means the line contains a sequence that needs to be appended to the current_sequence
        dic[current_head] = ''.join(current_sequence)
        return dic

# A method for drawing a figure with toolbar inside a TKinter window
def draw_figure_w_toolbar(canvas, fig, canvas_toolbar):
    # SOURCE : https://github.com/PySimpleGUI/PySimpleGUI/blob/master/DemoPrograms/Demo_Matplotlib_Embedded_Toolbar.py
    if canvas.children:
        for child in canvas.winfo_children():
            child.destroy()
    if canvas_toolbar.children:
        for child in canvas_toolbar.winfo_children():
            child.destroy()
    figure_canvas_agg = FigureCanvasTkAgg(fig, master=canvas)
    figure_canvas_agg.draw()
    toolbar = Toolbar(figure_canvas_agg, canvas_toolbar)
    toolbar.update()
    figure_canvas_agg.get_tk_widget().pack(side='right', fill='both', expand=1)
#                A toolbar for the figure inside a TKinter window
class Toolbar(NavigationToolbar2Tk):
  # SOURCE : https://github.com/PySimpleGUI/PySimpleGUI/blob/master/DemoPrograms/Demo_Matplotlib_Embedded_Toolbar.py
    def __init__(self, *args, **kwargs):
        super(Toolbar, self).__init__(*args, **kwargs)
        
#A method capable of checking if the diagonal of the windows contains at least as many corresponding characters as the user provided threshold value    
def countMatches(hwindow,vwindow):
    a = hwindow [0]
    b = vwindow [0]
    if(a == b):
        matchCounter = 0
        iterator = 0
        for c in hwindow:
            if(hwindow[iterator]==vwindow[iterator]):
                matchCounter +=1
            iterator +=1 
    else:
        return 0
    return matchCounter 
        
def dotPlotWindowed(sHorizontal,sVertical,window,threshold):
    """
    A dotplot algorythm upgraded with the moving-window integration feature
    """
    global outputArray                      #A reference to a global array of size [s1][s2] filled with zeroes at that time
    row = 0                                #An iterator counting the rows processed so far
    vWindow = []                           #An array representing the vertical part of the smoothening window
    for c2 in sVertical:                   #Iterating through chars in the vertical sequence
        vWindow.append(c2)                 #Adding to the vWindow one letter a time
        if(row>=window-1):                 #If the size of window exceeds user-input size,
            if(row>=window):               
                vWindow.pop(0)             #The first item in the window (at index [0]) gets popped
            hWindow=[]                     #An array representing the horizontal part of the smoothening window
            column = 0                     #An iterator counting the columns processed so far
            for c1 in sHorizontal:         #Iterating through every single character in each row
                hWindow.append(c1)         #Adding to the window one letter a time  
                if(column>=window-1):      #If the size exceeds the declared size
                    if(column>=window):
                        hWindow.pop(0)     #The first char in the window get popped
                    if countMatches(hWindow,vWindow) >= threshold:  #Using an external method to count the amount of matches in the smoothening window
                        outputArray[row-1][column-1]=1 #If the count of maches meets the required threshold, a 1 is being set in the array[row-1][column-1]
                column += 1                #The algorythm proceeds to the next column
        row+=1                             #The algorythm proceeds to the next row
    return outputArray

def needlemanWunsch(sHorizontal,sVertical,matchScore,missmatchPenalty,gapPenalty):
    
    #local variables
    lenHor = len(sHorizontal)
    lenVert = len(sVertical)

    align1 =''
    align2 = ''
    
    nwMatchCount = 0
    nwGapCount = 0
    nwMissmatchCount = 0
    nwChainLength = 0


    #Initialize an array filled with zeros
    Scores = np.zeros((lenVert + 1,lenHor + 1))
    #Fill in the first vertical row with numbers 0, -1, -2, -3, (...)
    Scores[:,0] = np.linspace(0,lenVert * gapPenalty, lenVert+1)
    #Fill in the first horizontal row with numbers 0, -1, -2, -3, (...)
    Scores[0,:] = np.linspace(0, lenHor * gapPenalty, lenHor+1)

    #An array for storing local scores (local[0] - diagonal, local [1] - horizontal, local [2] - vertical)
    local = np.zeros(3)

    #An array storing a trace back through the best scored route
    traceback = np.empty((lenVert + 1,lenHor + 1))
    traceback[:] = np.NaN
    traceback[0][0] = 0


    #Populating the scores array with values
    for i in range (lenVert):
        for j in range (lenHor):
            if sVertical[i] == sHorizontal[j]:
                local[0] = Scores[i][j] + matchScore
            else:
                local[0] = Scores[i][j] +missmatchPenalty
            local[1] = Scores[i][j+1] + gapPenalty
            local[2] = Scores[i+1][j] + gapPenalty
            localMax = np.max(local)
            Scores[i+1][j+1] = localMax

    # After the whole score array is populated, a most optimal way to the beginning needs to be found
    k = lenVert
    l = lenHor
    while k > 0 and l > 0:
        score       = Scores[k][l]
        score_up    = Scores[k -1 ][l]
        score_diag  = Scores[k-  1][l - 1]
        traceback[k][l] = 0
        nwChainLength+=1


        # Check which cell was used before we got to this particular cell.
        if score == score_diag + matchScore and sVertical[k-1] == sHorizontal[l-1]:
            align1 += sVertical[k-1]
            align2 += sHorizontal[l-1]
            nwMatchCount += 1
            k -= 1
            l -= 1

        elif score == score_diag + missmatchPenalty:
            align1 += sVertical[k-1]
            align2 += sHorizontal[l-1]
            nwMissmatchCount +=1
            k -= 1
            l -= 1

        elif score == score_up + gapPenalty:
            align1 += sVertical[k-1]
            align2 += '-'
            nwGapCount +=1
            k -= 1

        else:
            align1 += '-'
            align2 += sHorizontal[l-1]
            nwGapCount +=1
            l -= 1

    #Populating the alignment strings
    while k > 0:
        align1 += sVertical[k-1]
        align2 += '-'
        traceback[k][l] = 0
        nwGapCount +=1
        nwChainLength+=1
        k -= 1


    while l > 0:
        align1 += '-'
        align2 += sHorizontal[l-1]
        traceback[k][l] = 0 
        nwGapCount +=1
        nwChainLength+=1       
        l -= 1
    
    traceback[k][l] = 0



    #Reversing the alignment string so they will match the user-input order
    align1 = align1[::-1]
    align2 = align2[::-1]

    nwFinalScore = Scores[lenVert][lenHor]
    return [Scores, align1, align2,nwFinalScore,nwMatchCount,nwMissmatchCount,nwGapCount,nwChainLength,traceback]
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def smithWaterman(sHorizontal,sVertical,matchScore,missmatchPenalty,gapPenalty):
    
    #local variables
    lenHor = len(sHorizontal)
    lenVert = len(sVertical)

    align1 =''
    align2 = ''
    
    nwMatchCount = 0
    nwGapCount = 0
    nwMissmatchCount = 0
    nwChainLength = 0


    #Initialize an array filled with zeros
    Scores = np.zeros((lenVert + 1,lenHor + 1))

    #An array for storing local scores (local[0] - diagonal, local [1] - horizontal, local [2] - vertical)
    local = np.zeros(3)

    #An array storing a trace back through the best scored route
    traceback = np.empty((lenVert + 1,lenHor + 1))
    traceback[:] = np.NaN
    
    #Populating the scores array with values
    highScore = 0
    highScoreCount = 0
    for i in range (lenVert):
        for j in range (lenHor):
            if sVertical[i] == sHorizontal[j]:
                local[0] = Scores[i][j] + matchScore
            else:
                local[0] = Scores[i][j] +missmatchPenalty
            local[1] = Scores[i][j+1] + gapPenalty
            local[2] = Scores[i+1][j] + gapPenalty
            localMax = np.max(local)
            if localMax < 0:
                localMax = 0
            Scores[i+1][j+1] = localMax
            if Scores[i+1][j+1] == highScore:
                highScoreCount += 1
            if Scores[i+1][j+1] > highScore:
                highScore = localMax
                highScoreCount = 1
    
    solutions = np.argwhere(Scores == highScore)
    
    seq1list = []
    seq2list = []
    

    # After the whole score array is populated, a most optimal way to the beginning needs to be found
    for solution in solutions:
        verticalIndex = solution[0]
        horizontalIndex = solution[1]
        k = verticalIndex
        l = horizontalIndex
        traceback[k][l] = 0
        align1= ''
        align2=''
        while  Scores [k][l] > 0 :
            score       = Scores[k][l]
            score_up    = Scores[k -1 ][l]
            score_diag  = Scores[k-  1][l - 1]
            traceback[k][l] = 0
            nwChainLength+=1


            # Check which cell was used before we got to this particular cell.
            if score == score_diag + matchScore and sVertical[k-1] == sHorizontal[l-1]:
                align1 += sVertical[k-1]
                align2 += sHorizontal[l-1]
                nwMatchCount += 1
                k -= 1
                l -= 1

            elif score == score_diag + missmatchPenalty:
                align1 += sVertical[k-1]
                align2 += sHorizontal[l-1]
                nwMissmatchCount +=1
                k -= 1
                l -= 1

            elif score == score_up + gapPenalty:
                align1 += sVertical[k-1]
                align2 += '-'
                nwGapCount +=1
                k -= 1

            else:
                align1 += '-'
                align2 += sHorizontal[l-1]
                nwGapCount +=1
                l -= 1
        align1 = align1[::-1]
        align2 = align2[::-1]
        seq1list.append(align1)
        seq2list.append(align2)

        print(align1)
        print(align2)
    
    traceback[k][l] = 0



    #Reversing the alignment string so they will match the user-input order


    nwFinalScore = highScore
    return [Scores,nwFinalScore,traceback,seq1list,seq2list]



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ globals ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
global letterSet
letterSet = {}

global oneLtrAcids
oneLtrAcids = {"A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V"}
global dnaLtr
dnaLtr      = {"A","C","G","T"}

global fastaDictionary

global manualSelected, genbankSelected, fastaSelected, calculated, fastaFirstSequenceSelected, plotInverted
manualSelected, genbankSelected, fastaSelected, calculated, fastaFirstSequenceSelected, plotInverted = False, False, False, False, False, False

global sHorizontal, sVertical, alignHorizontal, alignVertical,matchString, finalScore, matchCount, missmatchCount, gapCount, chainLength
sHorizontal, sVertical, alignHorizontal, alignVertical,matchString, finalScore, matchCount, missmatchCount, gapCount, chainLength = "", "","","","", 0, 0, 0, 0, 0


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ layout ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

columnLeft = [
    [sg.Text("~~~~~~~~~~~~~~~~~~~~~~~~ INPUT MODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~",size=(63, 1))],
    [sg.B("manual",size=(62,1))],
    [sg.FileBrowse("FASTA",key="-path-",change_submits=True, target='-FASTAinput-',size=(62,1))],
    [sg.Button("genbank",size=(62,1))],
    [sg.I(key='-FASTAinput-', enable_events=True, visible=False)],

   
    [sg.Text("~~~~~~~~~~~~~~~~~~~~~~ SEQUENCE TYPE ~~~~~~~~~~~~~~~~~~~~~~~~~",size=(63, 1))],
    [sg.Radio('Nucleotides', "RADIO1", default=True, key="dictsNucleo"),
    sg.T("                                                                    "), 
    sg.Radio('Aminoacids', "RADIO1", default=False, key="dictsAmino")],

    [sg.Text("~~~~~~~~~~~~~~~~~~~~~~ ALGORYTHM TYPE ~~~~~~~~~~~~~~~~~~~~~~~~~",size=(63, 1))],
    [sg.Radio('Dot plot', "RADIO2", default=True, key="algosDotPlot"),
    sg.T("                  "), 
    sg.Radio('Needlemann-Wunsch', "RADIO2", default=False, key="algosNW"),
    sg.T("          "), 
    sg.Radio('Smith-Waterman', "RADIO2", default=False, key="algosSW")],

    [sg.Text("~~~~~~~~~~~~~~~~~~~~~~~~~~ SCORES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~",size=(63, 1))],
    [sg.T('match score:',      size=(55,1)),sg.I(key='-matchScore-',      size=(6,1),change_submits=True)],
    [sg.T('missmatch penalty:',size=(55,1)),sg.I(key='-missmatchPenalty-',size=(6,1),change_submits=True)],
    [sg.T('gap penalty:',      size=(55,1)),sg.I(key='-gapPenalty-',      size=(6,1),change_submits=True)],
    
    [sg.Text("~~~~~~~~~~~~~~~~~~~~~~~~~ CONTROLS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~",size=(63, 1))],
    [sg.B("Calculate",size=(14,1)),sg.B("Reset",size=(14,1)),sg.B("Save",size=(14,1)), sg.B("Exit",size=(14,1))]
    ]
 
columnRight = [
    [sg.T('~~~~~~~~~~~~~~~~~~~~~~~ DOT PLOT ~~~~~~~~~~~~~~~~~~~~~~~~~~')],
    [sg.T('Window size:',   size=(50,1)),sg.I(key='-winSize-',      size=(6,1),change_submits=True)],
    [sg.T('threshold:',     size=(50,1)),sg.I(key='-threshold-',    size=(6,1),change_submits=True)],
    
    [sg.T('~~~~~~~~~~~~~~~~~~~~~ MANUAL INPUT ~~~~~~~~~~~~~~~~~~~~~~~~')],
    [sg.T('First sequence:',key='firstlbl',size=(20,1)),sg.I(key='-firstSeq-', size=(41, 1))],
    [sg.T('Second sequence:',key='secondlbl',size=(20,1)),sg.I(key='-secondSeq-',size=(41, 1))],
    
    [sg.T('~~~~~~~~~~~~~~~~~~~~~ FASTA INPUT ~~~~~~~~~~~~~~~~~~~~~~~~~')],
    [sg.T('Select the first sequence',size=(50, 1), visible=False, key='listboxDesc')],
    [sg.Listbox(values = [],size=(64, 5), key='list',select_mode='single', change_submits='True')],

    [sg.T('~~~~~~~~~~~~~~~~~~~~~~~ NCBI INPUT ~~~~~~~~~~~~~~~~~~~~~~~~~')],
    [sg.T('database to search: (default - nucleotide)',       size=(34,1), visible = True, key = "-dbasesLabel-"),sg.I(size=(25,1),change_submits=True, visible = True, key='-blastDatabases-'  )],
    [sg.T('First sequence NCBI id:' ,key='ncbifirstlbl',size=(34,1)),sg.I(key='-ncbifirstSeq-', size=(25, 1),change_submits=True)],
    [sg.T('Second sequence NCBI id:',key='ncbisecondlbl',size=(34,1)),sg.I(key='-ncbisecondSeq-',size=(25, 1),change_submits=True)]
]

bottom = [
    [sg.T("first seq",size=(100,1), key="alignedStrHor",font='Courier 10')],
    [sg.T("match markers",size=(100,1), key="matchMarkers",font='Courier 10')],
    [sg.T("second seq",size=(100,1), key="alignedStrVer",font='Courier 10')]
]

plotFragment = [
    [sg.T('Figure:')],
    [sg.Column(
        layout=[
            [sg.Canvas(key='fig_cv',
                        # it's important that you set this size
                        size=(500 * 2, 450)
                        )]
            ],
            background_color='#DAE0E6',
            pad=(0, 0)
    )],
    [sg.Canvas(key='controls_cv')]
]

layout = [
    [sg.Frame(layout=columnLeft, title=''), sg.Frame(layout=columnRight, title='')],
    [sg.Frame(layout=plotFragment, title='')],
    [sg.Frame(layout = bottom,title='')]
    ]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ GUI ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
window = sg.Window("Demo", layout)
# Create an event loop
while True:
    event, values = window.read()
    # End program if user closes window or
    # presses the OK button
    if event == "Exit" or event == sg.WIN_CLOSED:
        break
    if event == "Reset":
        calculated = False
        manualSelected = False
        fastaSelected = False
        genbankSelected = False
        fastaFirstSequenceSelected = False

        sHorizontal = ""
        sVertical = ""

        window.FindElement('-path-').Update(button_color=('white', 'blue'))
        window.FindElement('manual').Update(button_color=('white', 'blue'))
        window.FindElement('genbank').Update(button_color=('white', 'blue'))

        window.FindElement('-firstSeq-').Update("")
        window.FindElement('-secondSeq-').Update("")
        window.FindElement('-blastDatabases-').Update("")
        window.FindElement('-ncbifirstSeq-').Update("")
        window.FindElement('-ncbisecondSeq-').Update("")

        window.FindElement('-matchScore-').Update("")
        window.FindElement('-missmatchPenalty-').Update("")
        window.FindElement('-gapPenalty-').Update("")

    elif event == "-FASTAinput-":
        calculated = False
        manualSelected = False
        fastaSelected = True
        genbankSelected = False
        
        window.FindElement('-path-').Update(button_color=('white', 'green'))
        window.FindElement('manual').Update(button_color=('white', 'blue'))
        window.FindElement('genbank').Update(button_color=('white', 'blue'))

        window.FindElement('-firstSeq-').Update("")
        window.FindElement('-secondSeq-').Update("")
        window.FindElement('-blastDatabases-').Update("")
        window.FindElement('-ncbifirstSeq-').Update("")
        window.FindElement('-ncbisecondSeq-').Update("")
        
        try:
            fastaDictionary = fasta2dict(values["-FASTAinput-"])
            window.FindElement('list'). Update(fastaDictionary)
            window.FindElement('listboxDesc').Update(visible=True)
            window.FindElement('listboxDesc').Update("Select first sequence")

        except:
            print("FASTA input error")


    elif event == 'list' and len(values['list']):     # if a list item is chosen
        sg.Popup('Selected ', values['list'])
        if fastaFirstSequenceSelected == False:
            subkey = values['list']
            key = subkey[0]
            sHorizontal = fastaDictionary.get(key)
            fastaFirstSequenceSelected = True
            window.FindElement('listboxDesc').Update("Select second sequence")

        elif  fastaFirstSequenceSelected == True:
            subkey = values['list']
            key = subkey[0]
            sVertical = fastaDictionary.get(key)
            window.FindElement('listboxDesc').Update("Fill in the params and click Calculate")
    
    elif event =="manual":
        calculated = False
        manualSelected = True
        fastaSelected = False
        genbankSelected = False
        
        window.FindElement('manual').Update(button_color=('white', 'green'))
        window.FindElement('-path-').Update(button_color=('white', 'blue'))
        window.FindElement('genbank').Update(button_color=('white', 'blue'))

        window.FindElement('-firstSeq-').Update("")
        window.FindElement('-secondSeq-').Update("")
        window.FindElement('-blastDatabases-').Update("")
        window.FindElement('-ncbifirstSeq-').Update("")
        window.FindElement('-ncbisecondSeq-').Update("")

    elif event =="genbank":
        calculated = False
        genbankSelected = True
        fastaSelected = False
        manualSelected = False
        
        window.FindElement('genbank').Update(button_color=('white', 'green'))
        window.FindElement('manual').Update(button_color=('white', 'blue'))
        window.FindElement('-path-').Update(button_color=('white', 'blue'))

        window.FindElement('-firstSeq-').Update("")
        window.FindElement('-secondSeq-').Update("")
        window.FindElement('-blastDatabases-').Update("")
        window.FindElement('-ncbifirstSeq-').Update("")
        window.FindElement('-ncbisecondSeq-').Update("")

    elif event == "Calculate":
        
        #If the user didn't select input method he / she is unable to use the calculate button
        if manualSelected == False and genbankSelected == False and fastaSelected == False or bool(values["-matchScore-"])== False or bool(values["-missmatchPenalty-"])== False or bool(values["-gapPenalty-"])== False:
            sg.Popup('Invalid action!')
            continue


        #If the user selected manual input, he / she needs to provide date for all the necessary input fields
        elif manualSelected == True:
            if bool(values["-firstSeq-"]) == False or bool(values["-secondSeq-"])== False:
                sg.Popup('Invalid action!')
                continue
            else:
                sHorizontal      = values["-firstSeq-"]
                sVertical        = values["-secondSeq-"]
        
        #If the user selected ganbank input, he / she needs to provide date for all the necessary input fields
        elif genbankSelected == True:
            if bool(values["-ncbifirstSeq-"])== False or bool(values["-ncbifirstSeq-"])== False:
                sg.Popup('Invalid action!')
                continue
            if bool(values['-blastDatabases-']):
                bDatabases = values["-blastDatabases-"]
            else:
                bDatabases = "nucleotide"

            idNCBI = values["-ncbifirstSeq-"]
            idNCBIsecond = values["-ncbisecondSeq-"]


            Entrez.email = "250065@student.pwr.edu.pl"  # Always tell NCBI who you are
            filename = "temp.fasta"

            # Downloading...
            try:
                net_handle = Entrez.efetch(
                    db=bDatabases, id=idNCBI, rettype="fasta", retmode="text"
                )
                out_handle = open(filename, "w")
                out_handle.write(net_handle.read())
                out_handle.close()
                net_handle.close()
                record = SeqIO.read(filename, "fasta")
                sHorizontal = str(record.seq)
                os.remove(filename)
            except:
                sg.popup("Bad Request!")
                continue

            try:
                net_handle = Entrez.efetch(
                    db="nucleotide", id=idNCBIsecond, rettype="fasta", retmode="text"
                )
                out_handle = open(filename, "w")
                out_handle.write(net_handle.read())
                out_handle.close()
                net_handle.close()
                record = SeqIO.read(filename, "fasta")
                sVertical = str(record.seq)
                os.remove(filename)
            except:
                sg.popup("Bad request!")
                continue


        if values["dictsNucleo"] == True:
            letterSet = dnaLtr
        elif values["dictsAmino"] == True:
            letterSet = oneLtrAcids

        inputCheck = True

        for c1 in sHorizontal:
            if c1.upper() in letterSet:
                continue
            else:
                inputCheck = False

        for c2 in sVertical:
            if c2.upper() in letterSet:
                continue
            else:
                inputCheck = False
        
        if inputCheck == False:
            sg.Popup("Invalid Input!")
            continue
        
        
        matchScore       = float(values["-matchScore-"])
        missmatchPenalty = float(values["-missmatchPenalty-"])
        gapPenalty       = float(values["-gapPenalty-"])

        if values["algosNW"] == True:
        
            result,align1,align2,nwFinalScore,nwMatchCount,nwMissmatchCount,nwGapCount,nwChainLength, traceback = needlemanWunsch(sHorizontal,sVertical,matchScore, missmatchPenalty, gapPenalty)
            finalScore = nwFinalScore
            matchCount = nwMatchCount
            gapCount = nwGapCount
            missmatchCount = nwMissmatchCount
            alignHorizontal = align1
            alignVertical = align2
            chainLength = nwChainLength
            

            matchString = ""
            for i in range (len(alignVertical)):
                if alignVertical[i] == alignHorizontal [i]:
                    matchString += "*"
                elif alignVertical[i] != "-" and alignHorizontal[i] != "-":
                    matchString += "|"
                else:
                    matchString += (" ") 

            window.FindElement('alignedStrHor').Update(align1.upper())
            window.FindElement('matchMarkers').Update(matchString)
            window.FindElement('alignedStrVer').Update(align2.upper())
            
        elif values["algosDotPlot"] == True:
            outputArray = [[0 for x in range(len(sHorizontal))] for y in range(len(sVertical))] 
            winSize     =  int(values["-winSize-"])
            threshold   = int(values["-threshold-"])
            result = dotPlotWindowed(sHorizontal,sVertical,winSize,threshold)

        # ------------------------------- PASTE YOUR MATPLOTLIB CODE HERE
        elif values["algosSW"] == True:
            result,swFinalScore, traceback, align1list,align2list = smithWaterman(sHorizontal,sVertical,matchScore, missmatchPenalty, gapPenalty)
        
        plt.figure(1)
        plt.clf()
        plt.xlabel("Sequence 2")
        plt.ylabel("Sequence 1")
        fig = plt.gcf()
        DPI = fig.get_dpi()
        # ------------------------------- you have to play with this size to reduce the movement error when the mouse hovers over the figure, it's close to canvas size
        fig.set_size_inches(505 * 2 / float(DPI), 404 / float(DPI))
        # -------------------------------
        if len(sVertical) < 100 and len(sHorizontal) < 100:
            c= plt.pcolor(result, edgecolors='k', linewidths=0.1)
            fig.colorbar(c)

        else:
            c = plt.pcolor(result)

        if values["algosNW"] == True:
            plt.pcolor(traceback, color="black")
            
        if values["algosSW"] == True:
            for i in range (len(traceback)):
                for j in range (len(traceback[i])):
                    if(traceback[i][j] == 0):
                        plt.plot(j + 0.5,i + 0.5,"ro")

        plt.gca().invert_yaxis()
        # ------------------------------- Instead of plt.show()
        draw_figure_w_toolbar(window['fig_cv'].TKCanvas, fig, window['controls_cv'].TKCanvas)

        calculated = True
    
    elif event == 'Save':
        if calculated == True and values["algosNW"] == True:
            
            files = [('Text Document', '*.txt'), ('All Files', '*.*')]
            filename = asksaveasfilename(filetypes = files, defaultextension = files)


            fileToSave = open(filename,"w")

            fileToSave.write("seq1:" + sHorizontal + "\n")
            fileToSave.write("seq2:" + sVertical + "\n")
            fileToSave.write("Match:" + str(matchCount) + "\n")
            fileToSave.write("Missmatch:" + str(missmatchCount) + "\n")
            fileToSave.write("Gap:" + str(gapCount) + "\n")
            fileToSave.write("Score:" + str(finalScore) + "\n")
            fileToSave.write("Length:" + str(chainLength) + "\n")
            fileToSave.write("Identity: " + str(matchCount) + "/" + str(chainLength) + " (" + str(int((matchCount / chainLength) * 100)) + "%)\n")
            fileToSave.write("Gaps: " + str(gapCount) + "/" + str(chainLength) + " (" + str(int((gapCount / chainLength) * 100)) + "%)\n")
            fileToSave.write(alignHorizontal + "\n")
        
            fileToSave.write(matchString + "\n")
            fileToSave.write(alignVertical)
            fileToSave.close()
            
            sg.Popup('Saved!')
            
        elif calculated == True and values["algosSW"] == True:
            
            files = [('Text Document', '*.txt'), ('All Files', '*.*')]
            filename = asksaveasfilename(filetypes = files, defaultextension = files)


            fileToSave = open(filename,"w")

            fileToSave.write("seq1:" + sHorizontal + "\n")
            fileToSave.write("seq2:" + sVertical + "\n")
            fileToSave.write("Score:" + str(swFinalScore) + "\n")
            
            
            
            for i in range(len(align1list)):
                
                matchString = ""
                align1 = align1list[i]
                align2 = align2list[i]
                chainLength = len(align1)
                matchCount = 0
                missmatchCount = 0
                gapCount = 0
                
                for j in range (len(align1)):
                    if align1[j] == align2 [j]:
                        matchString += "*"
                        matchCount += 1
                    elif align1[j] != "-" and align2[j] != "-":
                        matchString += "|"
                        missmatchCount += 1
                    else:
                        matchString += (" ") 
                        gapCount += 1
                
                fileToSave.write("\nOPTION " + str(i + 1) + ": \n")
                fileToSave.write("Match:" + str(matchCount) + "\n")
                fileToSave.write("Missmatch:" + str(missmatchCount) + "\n")
                fileToSave.write("Gap:" + str(gapCount) + "\n")
                fileToSave.write("Length:" + str(chainLength) + "\n")
                fileToSave.write("Identity: " + str(matchCount) + "/" + str(chainLength) + " (" + str(int((matchCount / chainLength) * 100)) + "%)\n")
                fileToSave.write("Gaps: " + str(gapCount) + "/" + str(chainLength) + " (" + str(int((gapCount / chainLength) * 100)) + "%)\n")
                fileToSave.write(align1 + "\n")
                fileToSave.write(matchString + "\n")
                fileToSave.write(align2 + "\n")
            fileToSave.close()
            window.FindElement('alignedStrHor').Update(align1.upper())
            window.FindElement('matchMarkers').Update(matchString)
            window.FindElement('alignedStrVer').Update(align2.upper())
            
            sg.Popup('Saved!')
            
            
            
        else:
            sg.Popup('Invalid action!')


            
            

window.close() 