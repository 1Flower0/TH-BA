import itertools
import time
import tkinter as tk
from tkinter import Menu, filedialog, messagebox, scrolledtext

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

from bm_search import PatternSearch as PS
from decodeDict import codonAminoDict as cad
class MyGUI:
    def __init__(self):
        self.root = tk.Tk()
        self.seqName = ""
        self.sequence = ""
        self.plotIndex = 0
        self.vicinityWidth = tk.IntVar()
        self.filterVal = tk.IntVar()
        self.longestCheck = tk.BooleanVar()
        self.cummulativeEval = tk.BooleanVar()
        self.root.geometry("1700x1000")
        self.root.title("Bitte Titel Ändern -- Wenn passender Name gefunden")
        self.root.config(background='lightgray')
        
        self.myMenu = Menu(self.root)
        self.openFile = Menu(self.myMenu,tearoff=False)
        self.openFile.add_command(label='Open File',command=self.openSeqFile)
        self.myMenu.add_cascade(label="File",menu=self.openFile)
        self.root.config(menu=self.myMenu)
       
        self.fig,self.axes = plt.subplots(4,1)
                
        self.fig.set_figheight(8.5)
        self.fig.set_figwidth(10)
        self.fig.subplots_adjust(left=0.05,bottom=0.1,right=0.95,top=0.925,hspace=0.75)

        self.leftFrame = tk.Frame(self.root,width=850,height=1000,background='lightgray')
        self.leftFrame.grid(row=0,column=0,padx=2,pady=5,sticky="news")

        self.rightFrame = tk.Frame(self.root,width=850,height=1000,background='lightgray')
        self.rightFrame.grid(row=0,column=1,padx=2,pady=5)

        self.canvas = FigureCanvasTkAgg(self.fig,self.rightFrame)
        self.canvas.get_tk_widget().grid(column=0,row=0,padx=5)
        self.rightFrame.place(relwidth=0.65,relx=0.4)

        self.toolbar = NavigationToolbar2Tk(self.canvas,self.rightFrame,pack_toolbar=False)
        self.toolbar.update()
        self.toolbar.grid(column=0,row=1,padx=5,pady=5)
        
        self.leftUpperFrame = tk.Frame(self.leftFrame,width=850,height=500,background='lightgray')
        self.leftUpperFrame.grid(column=0,row=0)
        self.leftLowerFrame = tk.Frame(self.leftFrame,width=850,height=500,background='lightgray')
        self.leftLowerFrame.grid(column=0,row=1)

        self.textplace = scrolledtext.ScrolledText(master=self.leftUpperFrame,width=92,height=25,borderwidth=3,font=("Arial",10))
        self.textplace.grid(column=0,row=0,padx=10,pady=10,sticky='N')

        self.motifFrame = tk.LabelFrame(master=self.leftLowerFrame,text="Motif",background='lightgray')
        self.motifFrame.grid(column=0,row=0,sticky="news")        
        self.vicinityFrame = tk.LabelFrame(master=self.leftLowerFrame,text="Umgebung",background='lightgray',width=800)
        self.vicinityFrame.grid(column=0,row=1,sticky="news")
        self.filterungFrame = tk.LabelFrame(master=self.leftLowerFrame,text="Filterung",background='lightgray',width=800)
        self.filterungFrame.grid(column=0,row=2,sticky="news")
        #LF1 SubFrames
        self.lf1_1 = tk.Frame(master=self.motifFrame,background='lightgray')
        self.lf1_1.grid(column=0,row=0,sticky="news")
        self.lf1_2 = tk.Frame(master=self.motifFrame,background='lightgray')
        self.lf1_2.grid(column=0,row=1,sticky="news")
        #LF1
        self.codonLabel = tk.Label(master=self.lf1_1,background="lightgray",text="Aus welchen Codons sollen Motifs erstellt werden?")
        self.codonLabel.grid(column=0,row=0)
        self.codons = tk.Entry(master=self.lf1_1)
        self.codons.grid(column=1,row=0,sticky="N")
        self.codons.bind('<Key>',self.checkCodonInput)
        self.kLabel = tk.Label(self.lf1_2,background='lightgray',text='Tupellänge bis:')
        self.kLabel.grid(column=0,row=0)
        self.kLowerLabel =tk.Label(self.lf1_2,background='lightgray',text='von: ')
        self.kLowerLabel.grid(column=0,row=1,sticky="w")
        self.kSpinboxUpper = tk.Spinbox(self.lf1_2,from_=0,to=5,width=10)
        self.kSpinboxUpper.grid(column=1,row=0)
        self.kSpinboxLower = tk.Spinbox(self.lf1_2,from_=0,to=5,width=10)
        self.kSpinboxLower.grid(column=1,row=1)
        self.findLongestMotif = tk.Checkbutton(self.lf1_2,text="Längste(s) Motif(s) finden",background='lightgray',onvalue=True,offvalue=False,variable=self.longestCheck,command=self.longestClicked)
        self.findLongestMotif.grid(column=2,row=0)
        self.cumEvalCheck = tk.Checkbutton(self.lf1_2,text='Kummulative Auswärtung',background='lightgray',onvalue=True,offvalue=False,variable=self.cummulativeEval)
        self.cumEvalCheck.grid(column=2,row=1)
        #LF2
        self.neighborLabel = tk.Label(self.vicinityFrame,background='lightgray',text="Umgebung:")
        self.neighborLabel.grid(column=0,row=0)
        self.codonNeighbors = tk.Spinbox(self.vicinityFrame,from_=0,to=50,width=10)
        self.codonNeighbors.grid(column=1,row=0)    
        self.codonRadio = tk.Radiobutton(master=self.vicinityFrame,background='lightgray',text='Codons',value=3,variable=self.vicinityWidth)    
        self.codonRadio.grid(column=0,row=1)
        self.dinukleoRadio = tk.Radiobutton(master=self.vicinityFrame,background='lightgray',text='Dinukleotide',value=2,variable=self.vicinityWidth)
        self.dinukleoRadio.grid(column=1,row=1)
        self.baseRadio = tk.Radiobutton(master=self.vicinityFrame,background='lightgray',text='Basen',value=1,variable=self.vicinityWidth)
        self.baseRadio.grid(column=2,row=1)
        #LF3
        self.allCodons = tk.Radiobutton(master=self.filterungFrame,background='lightgray',text='Alle Codons',value=1,variable=self.filterVal,command=self.customEnable)
        self.allCodons.grid(column=0,row=0)
        self.customCodons = tk.Radiobutton(master=self.filterungFrame,background='lightgray',text='Custom',value=2,variable=self.filterVal,command=self.customEnable)
        self.customCodons.grid(column=1,row=0)
        self.customCodonsInput = tk.Entry(master=self.filterungFrame,state="disabled")
        self.customCodonsInput.grid(column=1,row=2)

        self.pltBtn = tk.Button(self.leftLowerFrame,text="Evaluate Data",command=self.evalData,state="disabled")
        self.pltBtn.grid(column=0,row=3,padx=5)
        self.closeBtn = tk.Button(self.leftLowerFrame,command=self.closeWindow,text="Close")
        self.closeBtn.grid(column=0,row=4,ipadx=21)
        
        self.prevBtn = tk.Button(self.rightFrame,text="Zurück",command=self.prevPlot,state="disabled")
        self.prevBtn.grid(column=0,row=2,padx=5,pady=2)
        self.nxtBtn= tk.Button(self.rightFrame,text="Nächstes",command=self.nextPlot,state="disabled")
        self.nxtBtn.grid(column=0,row=3,padx=5,pady=2)

        self.root.mainloop()

    def filterFinds(self,myFinds:dict)->dict:
        myKeys = list(myFinds.keys())
        print(len(myKeys))
        myFilteredFinds = {}
        for index in range(len(myKeys)-1,0,-1):
            key = myKeys[index]
            print(key)
            keyLength = len(key)
            for upperElement in myFinds[key]:
                for lowerIndex in range(0,index):
                    lowerKey = myKeys[lowerIndex]                        
                    for lowerElement in myFinds[lowerKey]:
                        newList = []
                        if lowerElement == upperElement:
                            newList.append(lowerElement)                            
                        elif (lowerElement+len(lowerKey)) == (upperElement+keyLength):
                            newList.append(lowerElement)
                        elif lowerElement > upperElement and lowerElement+len(lowerKey) < upperElement+keyLength:
                            newList.append(lowerElement)
                        filteredList = [x for x in myFinds[lowerKey] if x not in newList]
                        myFilteredFinds[lowerKey] = filteredList
        myFilteredFinds = {k:v for k,v in myFilteredFinds.items() if not v == []}
        return myFilteredFinds

    def longestClicked(self):
        if self.longestCheck.get() == True:
            self.kSpinboxUpper.config(state='disabled')
        else:
            self.kSpinboxUpper.config(state='normal')

    def customEnable(self):
        if self.filterVal.get() == 2:
            self.customCodonsInput.config(state="normal")
        else:
            self.customCodonsInput.config(state="disabled")

    def closeWindow(self):
        self.warning = messagebox.askyesno("Fenster schließen?","Möchten Sie das Programm beenden?")
        if self.warning == 1:
            self.root.quit()
            self.root.destroy()
    
    def nFoundlingsWindow(self):
        msg = "".join([str(int(len(key)/self.vicinityWidth.get()))+"Codons "+str(key)+": "+str(len(self.cleanedIndicies[key]))+"\n" for key in self.cleanedIndicies.__reversed__()])
        messagebox.showinfo(title="Folgende Funde",message=msg)

    def nothingFound(self):
        messagebox.showerror(title="Keine Funde",message="Keine Teilsequenz gefunden")

    def notFoundInVicinityWarning(self,where:str):
        messagebox.showerror(title='Fehler',message='Gesuchte Codons konnten nicht im Bereich '+where+' gefunden werden')

    def checkCodonInput(self,event):
        if len(self.codons.get()) > 1 and not self.sequence == "":
            self.pltBtn.config(state="active")           
        else:
            self.pltBtn.config(state="disabled")

    def openSeqFile(self):
        self.textplace.config(state="normal")
        self.textplace.delete(1.0,tk.END)
        self.sequence=""
        newSeq = ""
        self.files = filedialog.askopenfilenames(filetypes=[('Fasta',['*.fasta','*.fna','*.fa'])])
        self.lastLineLength = []
        if len(self.files) == 1:
            file = self.files[0]
            f = np.genfromtxt(file,dtype=str,delimiter='\n',autostrip=True)
            self.seqName = f[0]
            self.sequence = ''.join([f[x] for x in range(1,len(f))])
            newSeq = ''.join([self.sequence[x:x+70]+'\n' for x in range(0,len(self.sequence),70)])
            self.lastLineLength.append(70-len(newSeq[len(newSeq)%71:-1]))
            newSeq+='+\n'
            self.textplace.insert(tk.INSERT,newSeq)
            self.textplace.config(state='disabled')
        elif len(self.files) > 1:
            for file in self.files:
                with open(file,'r') as f:
                    for line in f:
                        if(line[0]=='>'):
                            self.seqName = line
                            continue
                        if line[0] == '\n':
                            self.sequence+='\n+\n'
                            continue
                        if len(line) < 70:
                            self.lastLineLength.append((70-len(line)))
                        self.sequence+=line.strip()
            self.textplace.insert(tk.INSERT,self.sequence)
            self.textplace.config(state='disabled')        
            newSeq = ""
            for line in self.sequence:
                newSeq += line.strip()
            self.sequence = newSeq
    
    def translateCodons(self,sequence:str):
        if self.filterVal.get() == 1:
            print("Im in translate FilterVal1")
            if len(sequence)%3==0:
                print("Im in 0 mod 3") 
                return " ".join(cad[sequence[i:i+3]] for i in range(0,len(sequence),3))    
            else:
                return " ".join(cad[sequence[i:i+3]] if i+3<=len(sequence) and not sequence.__contains__('+') else '' for i in range(0,len(sequence),3))
        elif self.filterVal.get() == 2:
            print("Im in translate FilterVal2")
            return " ".join(cad[c] for c in self.customCodonArray)
             
    def createSearchPattern(self,txt,perms):
        print("Im in Search Now")
        print(len(perms))
        ps = PS
        incidies = {}
        if type(perms) == 'list':
            for perm in perms:
                part= "".join([x for x in perm])
                print("searching for motif",part)
                incidies[part]=ps.search(ps,txt,part)
        else:
            print("Im in else")
            for perm in perms:
                part= "".join([element for element in perm])
                incidies[part]=ps.search(ps,txt,part)
                
        if self.longestCheck.get() == False:
            return incidies
        else:
            myList = self.dynFinds(txt,incidies,self.codArr)
            longestLength = len(myList[-1][0])
            newDict = {}
            for entry in myList:
                if len(entry[0]) == longestLength:
                    newDict[entry[0]] = entry[1]
            return newDict

    def showNeighbors(self,txt:str,indexDict:dict):
        print("Im showing my Plot and Highlights")
        self.foundIndicies = list(indexDict.keys())
        self.clearPlot()
        self.keyLen = len(self.foundIndicies[self.plotIndex])
        if self.cummulativeEval.get() == False:            
            self.showPlot(indexDict[self.foundIndicies[self.plotIndex]],txt,self.keyLen)
            self.highlightFounds()
        else:
            self.cumulativeEval(txt,self.cleanedIndicies)


    def clearPlot(self):
        [x.clear() for x in self.axes]
        if self.vicinityWidth.get() == 1:
            self.axes[0].set_title("Basen davor")
            self.axes[1].set_title("Basen danach")
        elif self.vicinityWidth.get() == 2:
            self.axes[0].set_title("Dinukleoide davor")
            self.axes[1].set_title("Dinukleoide danach")
        else:               
            self.axes[0].set_title("Codons davor")
            self.axes[1].set_title("Codons danach")
            self.axes[2].set_title("Aminosäuren davor")
            self.axes[3].set_title("Aminosäuren danach")

    def dynFinds(self,text,oldDict,codonArray):
        myList = [(k,v) for k,v in oldDict.items() if not v == []]
        ps = PS
        for key in myList:
            for cod in codonArray:
                newKey = key[0]+cod
                returnValue = ps.search(ps,text,newKey)
                if not returnValue == []:
                    myList.append((newKey,returnValue))
                newKey = newKey[:-3]
        return myList 
    
    def cumulativeEval(self,text:str,motifIndicies:dict):
        codonsBefore = ""
        codonsAfter = ""
        delta = self.offset*self.vicinityWidth.get()
        for key,indexList in motifIndicies.items():
            keyLength = len(key)
            for index in indexList:
                # Filter Codons which occure before found Motif
                start = index-delta
                # Startpoint < 0 ?
                if start < 0:
                    # Get as close as posible to index 0
                    start = index%self.vicinityWidth.get()
                    codonsBefore += text[start:index]
                elif start == 0:
                    codonsBefore += ""
                elif text[index-1] == '+':# Motif is At Start of nth File
                    codonsBefore+=""
                elif '+' in text[start:index]:
                   plusIndex = text[start:index].find('+')
                   if not plusIndex == -1:
                        tmp = start+(start+plusIndex)
                        count = 0
                        while tmp > self.vicinityWidth.get():
                            tmp -= self.vicinityWidth.get()
                            count+=1
                        codonsBefore += text[index-(self.vicinityWidth.get()*count):index]
                else:
                    codonsBefore+=text[start:index]
                
                start = index+keyLength
                end = start+delta

                if start >= len(text):
                    codonsAfter +=""
                elif text[start] == '+' or text[start+1] == '+':
                    codonsAfter +=""
                elif '+' in text[start:end]:
                    plusIndex = text[start:end].find('+')
                    if not plusIndex == -1:
                        pass
                    
                else:
                    codonsAfter+=text[start:end]
                
        codonsBeforeDict = self.sequenceToDict(codonsBefore)
        codonsAfterDict = self.sequenceToDict(codonsAfter)
        aminosBeforeDict = self.aminosToDict(self.translateCodons(codonsBefore))
        aminosAfterDict = self.aminosToDict(self.translateCodons(codonsAfter))

        codonsBeforeRel = [x/sum(codonsBeforeDict.values()) for x in codonsBeforeDict.values()]
        codonsAfterRel = [x/sum(codonsBeforeDict.values()) for x in codonsAfterDict.values()]
        aminosBeforeRel = [x/sum(aminosBeforeDict.values()) for x in aminosBeforeDict.values()]
        aminosAfterRel = [x/sum(aminosAfterDict.values()) for x in aminosAfterDict.values()]
        
        self.axes[0].bar([*codonsBeforeDict.keys()],codonsBeforeRel,width=0.15)
        self.axes[1].bar([*codonsAfterDict.keys()],codonsAfterRel,width=0.15)
        
        self.axes[0].set_xticks(self.axes[0].get_xticks(), self.axes[0].get_xticklabels(), rotation=60, ha='right')
        self.axes[1].set_xticks(self.axes[1].get_xticks(), self.axes[1].get_xticklabels(), rotation=60, ha='right')
        
        self.axes[2].bar([*aminosBeforeDict.keys()],aminosBeforeRel,width=0.15)
        self.axes[3].bar([*aminosAfterDict.keys()],aminosAfterRel,width=0.15)
        
        self.axes[2].set_xticks(self.axes[2].get_xticks(), self.axes[2].get_xticklabels(), rotation=45, ha='right') 
        self.axes[3].set_xticks(self.axes[3].get_xticks(), self.axes[3].get_xticklabels(), rotation=45, ha='right')
        print("MainWindow")
        self.canvas.draw()
                
    def showPlot(self,indexlist,txt:str,keyLen:int):
        print("#####In showPlot Method#####")
        print("Current Motif",self.foundIndicies[self.plotIndex])
        print("Count of Motif found: ",len(indexlist))
        print([x for x in indexlist])
        self.clearPlot()
        codonsBefore = ""
        codonsAfter = ""
        for i,index in enumerate(indexlist):
            #Filter Codons which occure before found Motif
            if (index-(self.offset*self.vicinityWidth.get())) < 0: # Check if Vicinity os not out of Bounce | Negative
                start = index%self.vicinityWidth.get()
                codonsBefore += txt[start:index]
            elif txt[index-1] == '+': # Check if Character right before Motif is '+'
                codonsBefore+=""
            elif "+" in txt[index-(self.offset*self.vicinityWidth.get()):index]: #Check if Subsequence contains '+'
                plusIndex = txt[index-(self.offset*self.vicinityWidth.get()):index].find('+') #Read Sequence until Index right after '+'
                if not plusIndex == -1:
                    tmp = (self.offset*self.vicinityWidth.get())-plusIndex
                    count = 0
                    while tmp > self.vicinityWidth.get():
                        tmp -= self.vicinityWidth.get()
                        count+=1
                    codonsBefore += txt[index-(self.vicinityWidth.get()*count):index]
            else: #Everything is fine and i can just read the Subsequence
                codonsBefore += txt[index-(self.offset*self.vicinityWidth.get()):index]
        for index in indexlist:
            #Filter Codons which occure after found Motif
            if txt[index+keyLen-1] == '\n+' or txt[index+keyLen-1] == '+':
                codonsAfter+=""
            if not "+" in txt[index+keyLen:index+keyLen+self.offset*self.vicinityWidth.get()]: #Check if '+' in Subsequence after found Motif 
                if (index+keyLen+self.offset*self.vicinityWidth.get()) <= len(txt): #Check if Vicinity is not out of Bounce | Positive
                    codonsAfter += txt[index+keyLen:index+keyLen+self.offset*self.vicinityWidth.get()]
                else:   
                    codonsAfter += txt[index+keyLen:len(txt)-(len(txt)%self.vicinityWidth.get())]
            else: #Everything is fine and i can just read the Subsequence
                codonsAfter += txt[index+keyLen:index+keyLen+(self.offset*self.vicinityWidth.get())]


        codonsBeforeDict = self.sequenceToDict(codonsBefore)
        codonsAfterDict = self.sequenceToDict(codonsAfter)
        if int(self.vicinityWidth.get()) == 1:
            match int(self.filterVal.get()):
                case 1:
                    relBefore = [x/sum(codonsBeforeDict.values()) for x in codonsBeforeDict.values()]
                    relAfter = [x/sum(codonsAfterDict.values()) for x in codonsAfterDict.values()]
                    self.axes[0].bar([*codonsBeforeDict.keys()],relBefore,width=0.15)
                    self.axes[1].bar([*codonsAfterDict.keys()],relAfter,width=0.15)
        elif int(self.vicinityWidth.get()) == 2:
            match int(self.filterVal.get()):
                case 1:
                    relBefore = [x/sum(codonsBeforeDict.values()) for x in codonsBeforeDict.values()]
                    relAfter = [x/sum(codonsAfterDict.values()) for x in codonsAfterDict.values()]
                    self.axes[0].bar([*codonsBeforeDict.keys()],relBefore,width=0.15)
                    self.axes[1].bar([*codonsAfterDict.keys()],relAfter,width=0.15)
        elif int(self.vicinityWidth.get()) == 3:    
            match int(self.filterVal.get()):
                case 1:
                    print("Im in FilterVal1")
                    start = time.time()
                    aminosBeforeDict = self.aminosToDict(self.translateCodons(codonsBefore))
                    aminosAfterDict = self.aminosToDict(self.translateCodons(codonsAfter))
                    end = time.time()
                    print("I translated to AminoAcids",(end-start))
                    sumBeforeDict = sum(codonsBeforeDict.values())
                    relBefore = [x/sumBeforeDict for x in codonsBeforeDict.values()]
                    
                    sumAfterDict = sum(codonsAfterDict.values())
                    relAfter = [x/sumAfterDict for x in codonsAfterDict.values()]
                    print("Here comes the Plot")
                    start = time.time()
                    self.axes[0].bar([*codonsBeforeDict.keys()],relBefore,width=0.15)
                    self.axes[1].bar([*codonsAfterDict.keys()],relAfter,width=0.15)
                    #self.axes[1].sharex(self.axes[0])
                    self.axes[0].set_xticks(self.axes[0].get_xticks(), self.axes[0].get_xticklabels(), rotation=60, ha='right')
                    self.axes[1].set_xticks(self.axes[1].get_xticks(), self.axes[1].get_xticklabels(), rotation=60, ha='right')
                    
                    print('start')
                    print(relBefore)
                    print(relAfter)                    
                    print(self.axes[1].get_xticks())                    

                    sumAminoBefore = sum(aminosBeforeDict.values())
                    print("Sum of Aminos Before ",sumAminoBefore)
                    print("Printing Aminoacids with its count")
                    print(aminosBeforeDict)
                    relBefore = [x/sumAminoBefore for x in aminosBeforeDict.values()]
                    sumAminoAfter = sum(aminosAfterDict.values())
                    print("Sum of Aminos After ",sumAminoAfter)
                    print("Printing Aminoacids with its count")
                    print(aminosAfterDict)
                    relAfter = [x/sumAminoAfter for x in aminosAfterDict.values()]
                    self.axes[2].bar([*aminosBeforeDict.keys()],relBefore,width=0.15)
                    self.axes[3].bar([*aminosAfterDict.keys()],relAfter,width=0.15) 
                    self.axes[2].set_xticks(self.axes[2].get_xticks(), self.axes[2].get_xticklabels(), rotation=45, ha='right') 
                    self.axes[3].set_xticks(self.axes[3].get_xticks(), self.axes[3].get_xticklabels(), rotation=45, ha='right')
                    print('start')
                    print(relBefore)
                    print(relAfter)
                    print(self.axes[1].get_xticks())
                    print(self.axes[2].get_xticks())
                    print(self.axes[3].get_xticks())
                    print(self.axes[1].get_xticklabels())
                    end = time.time()
                    print("Plot is finished",(end-start))
                case 2:
                    print("Im in FilterVal2")
                    aminosBeforeDict = self.aminosToDict(self.translateCodons(codonsBefore))
                    aminosAfterDict = self.aminosToDict(self.translateCodons(codonsAfter))
                    sumRelBefore = sum(codonsBeforeDict.values())
                    if not sumRelBefore == 0:
                        relBefore = [x/sumRelBefore for x in codonsBeforeDict.values()]
                    else:
                        relBefore = [0 for x in codonsBeforeDict.values()]
                    relAfter = [x/sum(codonsAfterDict.values()) for x in codonsBeforeDict.values()]
                    self.axes[0].bar([*codonsBeforeDict.keys()],relBefore,width=0.15)
                    self.axes[1].bar([*codonsAfterDict.keys()],relAfter,width=0.15)
                    relBefore = [x/sum(aminosBeforeDict.values()) for x in aminosBeforeDict.values()]
                    relAfter = [x/sum(aminosAfterDict.values()) for x in aminosAfterDict.values()]
                    self.axes[2].bar([*aminosBeforeDict],[*relBefore],width=0.15)
                    self.axes[3].bar([*aminosAfterDict],[*relAfter],width=0.15)
                    print(self.axes)
        self.fig.savefig('test.png')                 
        print("Canvas should be drawn")
        self.canvas.draw()

    def nextPlot(self):
        if self.plotIndex == len(self.foundIndicies)-1:
            self.nxtBtn.config(state="disabled")
        self.plotIndex+=1
        self.showPlot(self.cleanedIndicies[self.foundIndicies[self.plotIndex]],self.sequence,self.keyLen)
        if str(self.prevBtn['state']) == 'disabled':
            self.prevBtn.config(state="active")
        self.highlightFounds()

    def prevPlot(self):
        self.plotIndex-=1
        self.showPlot(self.cleanedIndicies[self.foundIndicies[self.plotIndex]],self.sequence,self.keyLen)
        self.highlightFounds()
        if self.plotIndex == 0:
            self.prevBtn.config(state="disabled")
        if str(self.nxtBtn['state']) == 'disabled':
            self.nxtBtn.config(state="active")
        
    def sequenceToDict(self,subsequence:str)->dict:
        seqDict = dict()
        delta = self.vicinityWidth.get()
        print(delta)
        if delta == 3:
            if self.filterVal.get() == 1:
                print("Im in FilterVal1 in Dict creation")
                codonBase = list(cad.keys())
            elif self.filterVal.get() == 2:
                print("Im in FilterVal2 in Dict creation")
                codonBase = np.array([str.upper(x.strip()) for x in self.customCodonsInput.get().split(",")])
                print('Codon Base')
                print(codonBase)
                self.customCodonArray = codonBase
            for x in codonBase:
                if x not in seqDict:
                    seqDict[x] = 0
            for x in range(0,len(subsequence),delta):
                if subsequence[x:x+delta] in seqDict:
                    seqDict[subsequence[x:x+delta]]+=1
            print("Im returning the sequenceToDict Dict")
            return seqDict
        else:
            for x in range(0,len(subsequence),delta):
                if subsequence[x:x+delta] not in seqDict:
                    seqDict[subsequence[x:x+delta]] = 1
                else:
                    seqDict[subsequence[x:x+delta]]+=1
            return seqDict

    def aminosToDict(self,aminos:str)->dict:
        print("Im creating a Amino Dict")
        aminoDict = dict()
        print(aminos)
        for x in aminos.split():
            if x not in aminoDict:
                aminoDict[x] = 1
            else:
                aminoDict[x]+=1
        return aminoDict

    def combinationOfCodons(self,codonArray):
        if self.longestCheck.get() == False:
            myList = []
            if int(self.kSpinboxLower.get()) == int(self.kSpinboxUpper.get()):
                return np.array([x for x in itertools.product(codonArray,repeat=int(self.kSpinboxLower.get()))])
            else:
                for i in range(int(self.kSpinboxLower.get()),int(self.kSpinboxUpper.get())+1):
                    tmpList = list(itertools.product(codonArray,repeat=i))
                    myList.append(tmpList)
            return myList
        else:
            return np.array([x for x in itertools.product(codonArray,repeat=3)])

    def evalData(self):
        print("Im in Eval now")
        start = time.time()
        self.plotIndex = 0
        self.fig.suptitle(self.seqName)
        self.codArr = np.array([str.upper(x.strip()) for x in self.codons.get().split(",")])
        indicies = self.createSearchPattern(self.sequence,self.combinationOfCodons(self.codArr))
        self.cleanedIndicies = {k:v for k,v in indicies.items() if not v == []}
        #self.cleanedIndicies = self.filterFinds(self.cleanedIndicies)
        print(self.cleanedIndicies)
        print("I filtered my Finds")
        if len(self.cleanedIndicies) == 0:
            self.nothingFound()
            return
        self.nxtBtn.config(state="active")
        self.offset = int(self.codonNeighbors.get())
        self.showNeighbors(self.sequence,self.cleanedIndicies)
        end = time.time()
        print("Completely Finished in: ",end-start)
        self.nFoundlingsWindow()
        self.highlightFounds()
        
    def highlightFounds(self):
        self.textplace.tag_remove('highlight',1.0,tk.END)
        indexlist = self.cleanedIndicies[self.foundIndicies[self.plotIndex]]
        self.xlist = [i for i, ltr in enumerate(self.sequence) if ltr == '+']
        for index in indexlist:
            count=0
            for x in self.xlist:
                if x < index:
                   count+=1
                   index += self.lastLineLength[count-1] 
            tmp = (int(index)/70)
            tmp2 = index - int(tmp)*70
            begin_vk = int(tmp)+1
            end_nk = tmp2+self.keyLen
            tmpStr1 = str(begin_vk+count)+"."+str(tmp2)
            if end_nk > 70:    
                while(end_nk>70):
                    begin_vk+=1
                    end_nk -= 70
                tmpStr2 = str(begin_vk)+"."+str(end_nk)
            else:               
                tmpStr2 = str(int(tmp)+1+count)+"."+str(tmp2+len(self.foundIndicies[self.plotIndex]))
            self.textplace.tag_add('highlight',tmpStr1,tmpStr2)
            self.textplace.tag_config('highlight',background='green',foreground='black')
            self.textplace.see(float(tmpStr1))
MyGUI()