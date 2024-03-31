import streamlit as st
import os
import io
import shutil
from zipfile import ZipFile
import itertools
import numpy as np
import time
from math import *
import matplotlib.pyplot as plt, mpld3
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
#from Bio import SeqIO
from bm_search import PatternSearch as PS
from decodeDict import codonAminoDict as cad
import streamlit.components.v1 as components


#semsettin kullular


class MyGUI:
    def __init__(self):
        st.set_page_config(
            page_title="BA",
            page_icon="🧊",
            layout="wide",
            initial_sidebar_state="expanded")
        
        self.filterSeite, self.diagramSeite = st.columns([3, 8])
        self.vicinityWidth = ""
        self.filterVal = ""
        self.longestCheck :bool
        self.cummulativeEval :bool
        self.text_Content='''Foo, Bar
                              123, 456
                              789, 000
                            '''
        


        if "nr" not in st.session_state:
            st.session_state.nr = 0

        if "txtArea" not in st.session_state:
            st.session_state.txtArea=''

        if "info" not in st.session_state:
            st.session_state.info=''
       
        if "max" not in st.session_state:
            st.session_state.max = 0
       
        if "min" not in st.session_state:
            st.session_state.min = 0
                  
        if "codons_input" not in st.session_state:
            st.session_state.codons_input=''
            
        if "seqName" not in st.session_state:
            st.session_state.seqName=''
            
        if "umgebung" not in st.session_state:
            st.session_state.umgebung=0
       
        if "sequence" not in st.session_state:
            st.session_state.sequence=''

        if "newSeq" not in st.session_state:
            st.session_state.newSeq=''
       
        if "plotIndex" not in st.session_state:
            st.session_state.plotIndex=0
       
        if "eva_Data" not in st.session_state:
            st.session_state.eva_Data=True
       
        if "filter_Input" not in st.session_state:
            st.session_state.filter_Input = True
       
        if "cutoms_Input" not in st.session_state:
            st.session_state.cutoms_Input = ""
            
        if "tb_disabled" not in st.session_state:
            st.session_state.tb_disabled = False

        if "back" not in st.session_state:
            st.session_state.back = True
            
        if "next" not in st.session_state:
            st.session_state.next = True
            
        if "Custom_Codon" not in st.session_state:
            st.session_state.Custom_Codon = True
        if "counter" not in st.session_state:
            st.session_state.counter = 0
        if "zipi" not in st.session_state:
            st.session_state.zipi = True
        if "lo" not in st.session_state:
            st.session_state.lo = b'example content'
            
            
        self.fig,self.axes = plt.subplots(4,1)
        
        
        self.fig.set_figheight(6.3)
        self.fig.set_figwidth(13.9)
        self.fig.subplots_adjust(left=0.05,bottom=0.1,right=0.95,top=0.925,hspace=0.75)
            
        if "canvas" not in st.session_state:
            st.session_state.canvas = mpld3.fig_to_html(self.fig)
        if "canvas1" not in st.session_state:
           st.session_state.canvas1 = self.fig
        

        with self.filterSeite.container():
            st.session_state.up = st.file_uploader("Upload a FASTA file", type=["fasta", "fa","fna"], on_change=self.openSeqFile, key="upload")
            box1 = st.container(border=True)
            box1.write("Motiv")
            box1.text_input("Aus welchen Codons sollen Motifs erstellt werden?", key="codons_input" , on_change=self.checkCodonInput)
            motiv = box1.columns(2)
            motiv[1].write("")
            motiv[1].write("")
            self.longestCheck = motiv[1].checkbox("Längste(s) Motif(s) finden", value=False, key="disable", on_change=self.longestClicked)
            motiv[1].write("")
            self.cummulativeEval = motiv[1].checkbox("Kummulative Auswärtung",value=False) 
            motiv[0].number_input(label="Tupellänge bis:" , min_value=0 ,max_value=5, key="max", disabled=st.session_state.tb_disabled)
            
            motiv[0].number_input("von:", min_value=0 ,max_value=5,key= "min")
            
            box2 = st.container(border=True)
            box2.number_input("Umgebung:" , min_value=0 ,max_value=50,key="umgebung")
            self.vicinityWidth = box2.radio("", options=["Codons", "Dinukleotide", "Basen"], horizontal = True)       
                        
            box3 = st.container(border=True)
            box3.write("Filterung")
            self.filterVal = box3.radio("", options=["All Codons", "Custom"], horizontal = True, on_change=self.customEnable)     
            box3.text_input(label="Customs", key='cutoms_Input', disabled= st.session_state.filter_Input)
            box3.file_uploader("Upload txt", type=["txt", "doc"], key="cCodons", disabled=st.session_state.filter_Input,on_change=self.openCusFile)
            self.folders = box3.columns(2)
            self.folders[0].button(label='Zip erstellen',disabled=st.session_state.Custom_Codon, on_click=self.multipleQuest)
            self.folders[1].download_button(
            label="Download ZIP",
            key="ZIP",
            data=st.session_state.lo,
            file_name="myfile.zip",
            mime="application/zip",
            disabled=st.session_state.zipi
        )
            #self.folders[1].button(label="download", on_click=self.create_download_link_for_folder )# = self.create_download_button(temp_folder)
            but = st.columns(2)
            but[0].button(label="Evaluate Data", disabled = st.session_state.eva_Data, on_click= self.evalData)
            but[1].button(label="Refresh Page", on_click = st.session_state.clear) 
            
             
        with self.diagramSeite.container():
            st.text_area(label="sequence", label_visibility="hidden", value=st.session_state.newSeq,key=st.session_state.txtArea, height=5)

            
            
            #components.html(st.session_state.canvas, height=700)
            st.pyplot(st.session_state.canvas1)
            
            uDiagram = st.columns(2)
            uDiagram[0].button(label="Zurück", disabled=st.session_state.back, on_click=self.prevPlot)
            uDiagram[1].button(label="Nächstes", disabled=st.session_state.next, on_click= self.nextPlot)    
 


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
        if self.longestCheck == True:
            st.session_state.tb_disabled = False
        else:
            st.session_state.tb_disabled = True

    def customEnable(self):
        if self.filterVal =='Custom':
            st.session_state.filter_Input = True
        else:
            st.session_state.filter_Input = False

    def nFoundlingsWindow(self):
        try:
            msg = "Folgende Funde".join([str(int(len(key)/st.session_state.nr))+"Codons "+str(key)+": "+str(len(self.cleanedIndicies[key]))+ "\n" for key in self.cleanedIndicies.__reversed__()])    
            if st.session_state.Custom_Codon:
                st.success( icon="💯", body=msg)
            else:
                msg=msg+"\n"+"Line "+str(st.session_state.counter)
                st.toast( icon="💯", body=msg)
        except:
            if st.session_state.Custom_Codon:
                st.error(icon= "🚨", body="Keine Teilsequenz gefunden")
            else:
                st.toast(icon= "🚨", body="Keine Teilsequenz gefunden")
            
            print('#######################################')
            print("object has no attribute cleanedIndicies")
            print('#######################################')
        
      
    def nothingFound(self):
        st.error(icon= "🚨", body="Keine Teilsequenz gefunden")

    def notFoundInVicinityWarning(self,where:str):
        st.error( icon=  "🚨", body='Gesuchte Codons konnten nicht im Bereich '+where+' gefunden werden')
  
    def checkCodonInput(self):
        if len(st.session_state.codons_input) > 1 and not st.session_state.newSeq == "":
            st.session_state.eva_Data = False          
        else:
            st.session_state.eva_Data = True     

    def openSeqFile(self):   
        st.session_state.sequence=""
        newSeq = ""
        self.lastLineLength = []
        files = st.session_state.upload
        counter = 1

        #if tb == 1:
        if True:
            file = files
            try:
                f = np.genfromtxt(files,dtype=str,delimiter='\n',autostrip=True)
                st.session_state.seqName = f[0]
                st.session_state.sequence = ''.join([f[x] for x in range(1,len(f))])
                st.session_state.info = f
                newSeq = ''.join([st.session_state.sequence[x:x+70]+'\n' for x in range(0,len(st.session_state.sequence),70)])
                self.lastLineLength.append(70-len(newSeq[len(newSeq)%71:-1]))
                newSeq+='+\n'
                print(newSeq)
                st.session_state.newSeq=newSeq
            except:
                st.session_state.newSeq =''
                print('#######################################')
                print("datei entfernt")
                print('#######################################')
        elif counter > 1:
            for file in files:
                with open(file,'r') as f:
                    for line in f:
                        if(line[0]=='>'):
                            self.seqName = line
                            continue
                        if line[0] == '\n':
                            st.session_state.sequence+='\n+\n'
                            continue
                        if len(line) < 70:
                            self.lastLineLength.append((70-len(line)))
                        st.session_state.sequence+=line.strip()     
            newSeq = ""
            for line in st.session_state.sequence:
                newSeq += line.strip()
            st.session_state.sequence = newSeq

    def translateCodons(self,sequence:str):
        if self.filterVal == "All Codons":
            print("Im in translate FilterVal1")
            if len(sequence)%3==0:
                print("Im in 0 mod 3") 
                return " ".join(cad[sequence[i:i+3]] for i in range(0,len(sequence),3))    
            else:
                return " ".join(cad[sequence[i:i+3]] if i+3<=len(sequence) and not sequence.contains('+') else '' for i in range(0,len(sequence),3))
        elif self.filterVal == "Custom":
            print("Im in translate FilterVal2")
            return " ".join(cad[c] for c in self.customCodonArray)
    
    def createSearchPattern(self,txt,perms):
        print("Im in Search Now")
        print(len(perms))
        ps = PS
        incidies = {}
        try:
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
                    
            if st.session_state.tb_disabled == False:
                return incidies
            else:
                myList = self.dynFinds(txt,incidies,self.codonArray)
                longestLength = len(myList[-1][0])
                newDict = {}
                for entry in myList:
                    if len(entry[0]) == longestLength:
                        newDict[entry[0]] = entry[1]
                return newDict
        except:
            print('#######################################')
            print("sequence item 0: expected str instance, tuple found")
            print('#######################################')
    
    def showNeighbors(self,txt:str,indexDict:dict):
        print("Im showing my Plot and Highlights")
        self.foundIndicies = list(indexDict.keys())
        self.clearPlot()
        self.keyLen = len(self.foundIndicies[self.plotIndex])
        if self.cummulativeEval == False:            
            self.showPlot(indexDict[self.foundIndicies[self.plotIndex]],txt,self.keyLen)
            self.highlightFounds()
        else:
            self.cumulativeEval(txt,self.cleanedIndicies)
 
    def clearPlot(self):
        #[x.clear for x in self.axes]
        if self.vicinityWidth == "Basen":
            self.axes[0].set_title("Basen davor")
            self.axes[1].set_title("Basen danach")
        elif self.vicinityWidth == "Dinukleotide":
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

    ### yapilmasi gerek
    def cumulativeEval(self,text:str,motifIndicies:dict):
        if self.vicinityWidth == "Codons": st.session_state.nr = 3
        elif self.vicinityWidth == "Dinukleotide":st.session_state.nr = 2
        else :st.session_state.nr = 1
        codonsBefore = ""
        codonsAfter = ""
        delta = self.offset*st.session_state.nr
        for key,indexList in motifIndicies.items():
            keyLength = len(key)
            for index in indexList:
                # Filter Codons which occure before found Motif
                start = index-delta
                # Startpoint < 0 ?
                if start < 0:
                    # Get as close as posible to index 0
                    start = index%st.session_state.nr
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
                        while tmp > st.session_state.nr:
                            tmp -= st.session_state.nr
                            count+=1
                        codonsBefore += text[index-(st.session_state.nr*count):index]
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
                
        if st.session_state.Custom_Codon==True:
            codonsBeforeDict=""
            codonsAfterDict=""
            codonsBeforeRel =""
            codonsAfterRel = ""
            aminosBeforeRel = ""
            aminosAfterRel = ""

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
        print("Hauptfenster")
        st.session_state.canvas1 =self.fig
        
   
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
            if (index-(self.offset*st.session_state.nr)) < 0: # Check if Vicinity os not out of Bounce | Negative
                start = index%st.session_state.nr
                codonsBefore += txt[start:index]
            elif txt[index-1] == '+': # Check if Character right before Motif is '+'
                codonsBefore+=""
            elif "+" in txt[index-(self.offset*st.session_state.nr):index]: #Check if Subsequence contains '+'
                plusIndex = txt[index-(self.offset*st.session_state.nr):index].find('+') #Read Sequence until Index right after '+'
                if not plusIndex == -1:
                    tmp = (self.offset*st.session_state.nr)-plusIndex
                    count = 0
                    while tmp > st.session_state.nr:
                        tmp -= st.session_state.nr
                        count+=1
                    codonsBefore += txt[index-(st.session_state.nr*count):index]
            else: #Everything is fine and i can just read the Subsequence
                codonsBefore += txt[index-(self.offset*st.session_state.nr):index]
        for index in indexlist:
            #Filter Codons which occure after found Motif
            if txt[index+keyLen-1] == '\n+' or txt[index+keyLen-1] == '+':
                codonsAfter+=""
            if not "+" in txt[index+keyLen:index+keyLen+self.offset*st.session_state.nr]: #Check if '+' in Subsequence after found Motif 
                if (index+keyLen+self.offset*st.session_state.nr) <= len(txt): #Check if Vicinity is not out of Bounce | Positive
                    codonsAfter += txt[index+keyLen:index+keyLen+self.offset*st.session_state.nr]
                else:   
                    codonsAfter += txt[index+keyLen:len(txt)-(len(txt)%st.session_state.nr)]
            else: #Everything is fine and i can just read the Subsequence
                codonsAfter += txt[index+keyLen:index+keyLen+(self.offset*st.session_state.nr)]

        codonsBeforeDict = self.sequenceToDict(codonsBefore)
        codonsAfterDict = self.sequenceToDict(codonsAfter)
        filterVal = 0
        if self.filterVal == "All Codons" : filterVal=1
        if self.filterVal == "Custom" : filterVal=2
        if st.session_state.nr == 1:
            match filterVal:
                case 1:
                    relBefore = [x/sum(codonsBeforeDict.values()) for x in codonsBeforeDict.values()]
                    relAfter = [x/sum(codonsAfterDict.values()) for x in codonsAfterDict.values()]
                    self.axes[0].bar([*codonsBeforeDict.keys()],relBefore,width=0.15)
                    self.axes[1].bar([*codonsAfterDict.keys()],relAfter,width=0.15)
        elif st.session_state.nr == 2:
            match filterVal:
                case 1:
                    relBefore = [x/sum(codonsBeforeDict.values()) for x in codonsBeforeDict.values()]
                    relAfter = [x/sum(codonsAfterDict.values()) for x in codonsAfterDict.values()]
                    self.axes[0].bar([*codonsBeforeDict.keys()],relBefore,width=0.15)
                    self.axes[1].bar([*codonsAfterDict.keys()],relAfter,width=0.15)
        elif st.session_state.nr == 3:    
            match filterVal:
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
        self.fig.savefig('Hatest.png')  
        print("Canvas should be drawn")
        #st.session_state.canvas = mpld3.fig_to_html(self.fig)
        st.session_state.canvas1 = self.fig


    def nextPlot(self):
        if self.plotIndex == len(self.foundIndicies)-1:
            st.session_state.next=False
        self.plotIndex+=1
        self.showPlot(self.cleanedIndicies[self.foundIndicies[self.plotIndex]],st.session_state.sequence,self.keyLen)
        if st.session_state.back == False:
            st.session_state.back == True
        self.highlightFounds()  

    def prevPlot(self):
        self.plotIndex-=1
        self.showPlot(self.cleanedIndicies[self.foundIndicies[self.plotIndex]],self.sequence,self.keyLen)
        self.highlightFounds()
        if self.plotIndex == 0:
           st.session_state.back = True
        if st.session_state.next == False:
            st.session_state.next = True
    
    def sequenceToDict(self,subsequence:str)->dict:
        seqDict = dict()
        delta = st.session_state.nr
        s:str
        if delta == 3:
            if self.filterVal == "All Codons":
                print("Im in FilterVal1 in Dict creation")
                codonBase = list(cad.keys())
            elif self.filterVal == "Custom":
                s = st.session_state.cutoms_Input
                print("Im in FilterVal2 in Dict creation")
                codonBase = np.array([str.upper(x.strip()) for x in s.split(" ")])
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
        if st.session_state.tb_disabled == False:
            myList = []
            if st.session_state.min == st.session_state.max:
                return np.array([x for x in itertools.product(codonArray,repeat=st.session_state.min)])
            else:
                for i in range(st.session_state.min,st.session_state.max+1):
                    tmpList = list(itertools.product(codonArray,repeat=i))
                    myList.append(tmpList)
            return myList
        else:
            return np.array([x for x in itertools.product(codonArray,repeat=3)])    
                         
    def evalData(self):
        if self.vicinityWidth == "Codons": st.session_state.nr = 3
        elif self.vicinityWidth == "Dinukleotide":st.session_state.nr = 2
        else :st.session_state.nr = 1
        print("Im in Eval now")
        start = time.time()
        self.plotIndex = 0
        self.fig.suptitle(st.session_state.seqName)
        s=st.session_state.codons_input
        try:
            self.codonArray = np.array([str.upper(x.strip()) for x in str(s).split("+")])
            indicies = self.createSearchPattern(st.session_state.sequence, self.combinationOfCodons(self.codonArray))
            self.cleanedIndicies = {k:v for k,v in indicies.items() if not v == []}
            #self.cleanedIndicies = self.filterFinds(self.cleanedIndicies)
            print(self.cleanedIndicies)
            print("I filtered my Finds")
            if len(self.cleanedIndicies) == 0:
                self.nothingFound()
                return
            st.session_state.next = True
            self.offset = int(st.session_state.umgebung)
            self.showNeighbors(st.session_state.sequence,self.cleanedIndicies)
            end = time.time()
            print("Completely Finished in: ",end-start)
        except:
            print('#######################################')
            print('object has no attribute items')
            print('#######################################')
        self.nFoundlingsWindow()
        self.highlightFounds()   

        
    def highlightFounds(self):
        try:
            indexlist = self.cleanedIndicies[self.foundIndicies[self.plotIndex]]
            self.xlist = [i for i, ltr in enumerate(st.session_state.sequence) if ltr == '+']
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
                #annotated_text((tmpStr1), background='green',foreground='black' )
        except:
            print('#######################################')
            print('highligth')
            print('#######################################')


    def multipleQuest(self):
        try:
            s=st.session_state.newSeq
            with open("folder/readme.txt", "w") as f:
                f.write("Die genutzte Sequence"+'\n')
                f.write(st.session_state.newSeq)
                f.write("MOTIV:"+'\n')
                f.write(st.session_state.codons_input+'\n')
                f.write("TUPELLÄNGE: VON-BIS"+'\n')
                f.write(str(st.session_state.min)+"-"+str(st.session_state.max)+'\n')
                f.write("Umgebung:"+'\n')
                f.write("")
                f.write("Typ:"+'\n')
                f.write(self.vicinityWidth+'\n')
            for line in st.session_state.cCodons:
                st.session_state.counter = st.session_state.counter + 1
                print(st.session_state.counter)
                print(line.decode("utf-8").strip())
                st.session_state.cutoms_Input=line.decode("utf-8").strip()
                self.evalData()
                s=st.session_state.counter
                self.fig.savefig('folder/'+str(s)+st.session_state.cutoms_Input+'.png')
                st.session_state.cutoms_Input=""
                self.axes[0].cla()
                self.axes[1].cla()
                self.axes[2].cla()
                self.axes[3].cla()

                        
            #schleife einfügen und die deien Speichern 
            st.session_state.counter=0
            st.session_state.zipi=False
            self.create_download_link_for_folder()

            
            
        except:
            st.write("datei")
            print('Datei erstellung')
            
    def openCusFile(self):
        if st.session_state.cCodons== 0:
            st.session_state.Custom_Codon = True          
        else:
            st.session_state.Custom_Codon = False  
            
            
    def create_download_link_for_folder(self):
        folder_path = "folder"
        
        st.write("hallo<")
         # Create in-memory buffer
        zip_buffer = io.BytesIO()

        # Zip the folder into the buffer
        with ZipFile(zip_buffer, "w") as zipf:
            for root, _, files in os.walk(folder_path):
                for file in files:
                    zipf.write(os.path.join(root, file), os.path.relpath(os.path.join(root, file), folder_path))

        # Serve the zip file for download
        zip_buffer.seek(5)

        st.markdown(
            f'<a href="folder_download.zip" download>Click here to download the folder</a>',
            unsafe_allow_html=True
        )
        st.session_state.lo=zip_buffer
        st.session_state.zipi=False
       # st.session_state.ZIP=self.folders[1].download_button(
       #     label="Download ZIP",
       #     key="ZIP",
       #     data=zip_buffer,
       #     file_name="myfile.zip",
       #     mime="application/zip",
       # )
        
        
            
MyGUI()