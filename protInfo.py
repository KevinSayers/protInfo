from Bio.SeqUtils.ProtParam import ProteinAnalysis
import csv

#Kevin Sayers

#function getSequence
#Input: list of each line in a given accession record
#OutPut: return the complete sequence string
def getSequence(recordList):
    startSeq = 0
    endSeq = 0
    seqStr = ""
    
    for i, line in enumerate(recordList):
        lineString = line
        processList = []
        processList = lineString.split('   ')
        if(processList[0] == 'SQ'):
            startSeq = i
        endSeq = i

    startSeq = startSeq + 1
    for j in range(startSeq, endSeq+1):
        tempStr = str(recordList[j])
        tempList = []
        tempList = tempStr.split()
        for g in tempList:
            seqStr += g


    return seqStr

#function getName
#Input: string for individual chains from accession
#Output: string with the name associated with a given chain
def getName(chain):
    nameStr = ""
    tempList = []

    tempList = chain.split()
    for j in range(2, len(tempList)):
        nameStr = nameStr + " " + str(tempList[j])


    return nameStr
    
#function getID
#Input: list of lines for a given accession
#Output: string with the accession ID
def getID(recordList):
    idStr = ""
    tempStr = str(recordList[1])
    idList = tempStr.split()
    tempStr = idList [1]
    idList = tempStr.split(';')
    idStr = idList[0]
    return idStr

#function getChain
#Input: the list of lines for a given accession record
#Output: a list of chains for a accession record including sequence range and name
def getChain(recordList):
    chainList = []
    lineList = []
    for line in recordList:
        finalOut = ""
        finalList = []
        lineString = line
        processLIst = []
        processList = lineString.split('   ')
        if(processList[0] == 'FT' and processList[1] == 'CHAIN'):
            lineAfterChain = True
            lineList = line.split()
            chainStr1 = ""
            chainStr2 = ""
            outString = "" 
            for i in range(2,len(processList)):
                chainStr1 += str(processList[i]) + " "
            tempIndex = (recordList.index(line))+1
            for j in range(0,2):
                tempList = str(recordList[tempIndex+j]).split('   ')
                if(tempList[1] == ''):
                    for l in range(2,len(tempList)):
                        chainStr2 += str(tempList[l]) + " "
            outString += chainStr1 + " " + chainStr2


            chainList.append(outString)

    return chainList
#function getSubSeq
#Input: a chain for an accession, and the total sequence
#Output: the sequence 
def getSubSeq(chain, sequence):
    seqStart = 0
    seqEnd = 0
    tempChainList = []
    tempSeqList = []

    #final output
    subSeq = []
    subStr = ""
    
    try:
        tempChainList = chain.split()
        seqStart = int(tempChainList[0])
        seqEnd = int(tempChainList[1])

        tempSeqList = list(sequence)
        for i in range(seqStart-1, seqEnd):
            subSeq.append(str(tempSeqList[i]))

        for j in subSeq:
            subStr += str(j)
    except ValueError:
        subStr = "ERROR"

    return subStr
    #print len(tempSeqList)
    #print seqStart
    #print seqEnd

#Function getMF
#Input: subsequence string
#Output: string with molecular formula
def getMF(subSeq):
    listofaminoacids = []
    #Dictionary for each amino acid with atoms for each
    A = {'C':3, 'H':7, 'N':1, 'O':2, 'S':0}
    R = {'C':6, 'H':14,'N':4, 'O':2, 'S':0}
    N = {'C':4, 'H':8, 'N':2, 'O':3, 'S':0}
    D = {'C':4, 'H':7, 'N':1, 'O':4, 'S':0}
    C = {'C':3, 'H':7, 'N':1, 'O':2, 'S':1}
    Q = {'C':5, 'H':10,'N':2, 'O':3, 'S':0}
    E = {'C':5, 'H':9, 'N':1, 'O':4, 'S':0}
    G = {'C':2, 'H':5, 'N':1, 'O':2, 'S':0}
    H = {'C':6, 'H':9, 'N':3, 'O':2, 'S':0}
    I = {'C':6, 'H':13,'N':1, 'O':2, 'S':0}
    L = {'C':6, 'H':13,'N':1, 'O':2, 'S':0}
    K = {'C':6, 'H':14,'N':2, 'O':2, 'S':0}
    M = {'C':5, 'H':11,'N':1, 'O':2, 'S':1}
    F = {'C':9, 'H':11,'N':1, 'O':2, 'S':0}
    P = {'C':5, 'H':9, 'N':1, 'O':2, 'S':0}
    S = {'C':3, 'H':7, 'N':1, 'O':3, 'S':0}
    T = {'C':4, 'H':9, 'N':1, 'O':3, 'S':0}
    W = {'C':11,'H':12,'N':2, 'O':2, 'S':0}
    Y = {'C':9, 'H':11,'N':1, 'O':3, 'S':0}
    V = {'C':5, 'H':11,'N':1, 'O':2, 'S':0}
    
    dictOfAmino = {'A':A,'R':R,'N':N,'D':D,'C':C,'Q':Q, 'E':E, 'G':G,'H':H,'I':I,'L':L,'K':K,'M':M,'F':F,'P':P,'S':S,'T':T,'W':W,'Y':Y,'V':V}
    mySeq = subSeq
    analysis = ProteinAnalysis(mySeq)
    listofaminoacids.append(analysis.count_amino_acids())

    for i in listofaminoacids:
        carbonTotal = 0
        hydrogenTotal = 0
        oxygenTotal = 0
        nitrogenTotal = 0
        sulfurTotal = 0
        peptideBonds = 0
        
        for value in i:
                for amino in dictOfAmino:
                        
                        if value == amino:
                                peptideBonds = peptideBonds + i[value]
                                thisAmino = {}
                                thisAmino = dictOfAmino[amino]
                                carbonTotal = carbonTotal + (i[value]*thisAmino['C'])
                                hydrogenTotal = hydrogenTotal + (i[value]*thisAmino['H'])
                                oxygenTotal = oxygenTotal + (i[value]*thisAmino['O'])
                                nitrogenTotal = nitrogenTotal + (i[value]*thisAmino['N'])
                                sulfurTotal = sulfurTotal + (i[value]*thisAmino['S'])
                                                             

        #Correcting totals for peptide bond loss of water
        peptideBonds = peptideBonds - 1
        hydrogenTotal = hydrogenTotal -(peptideBonds*2)
        oxygenTotal = oxygenTotal - (peptideBonds*1)
        outString = "C" + str(carbonTotal) + "H" + str(hydrogenTotal) + "N" + str(nitrogenTotal) + "O" + str(oxygenTotal) + "S" + str(sulfurTotal)
        return outString

    
def getMW_mono(subSeq):
    peptideBonds = 0
    molecularWeight = 0.0
    waterLoss = 18.015
    
    listofaminoacids = []

    #MONOISOTOPIC MW FOR EACH AMINO ACID CURRENTLY
    dictOfAmino = {'A':71.03711,
                   'R':156.10111,
                   'N':114.04293,
                   'D':115.02694,
                   'C':103.00919,
                   'Q':128.05858,
                   'E':129.04259,
                   'G':57.02146,
                   'H':137.05891,
                   'I':113.08406,
                   'L':113.08406,
                   'K':128.09496,
                   'M':131.04049,
                   'F':147.06841,
                   'P':97.05276,
                   'S':87.03203,
                   'T':101.04768,
                   'W':186.07931,
                   'Y':163.06333,
                   'V':99.06841}
    mySeq = subSeq
    analysis = ProteinAnalysis(mySeq)
    listofaminoacids.append(analysis.count_amino_acids())

    for i in listofaminoacids:
        for value in i:
            for amino in dictOfAmino:
                if value == amino:
                    peptideBonds = peptideBonds + i[value]
                    #print dictOfAmino[value]
                    #print i[value]
                    molecularWeight = molecularWeight + (i[value]*dictOfAmino[value])

    #peptideBonds = peptideBonds - 1 
    #molecularWeight = molecularWeight - (peptideBonds*waterLoss)
    molecularWeight =  molecularWeight+waterLoss
    return molecularWeight

def getMW_average(subSeq):
    peptideBonds = 0
    molecularWeight = 0.0
    waterLoss = 18.015
    
    listofaminoacids = []

    #AVERAGE MW FOR EACH AMINO ACID CURRENTLY
    dictOfAmino = {'A':71.0788,
                   'R':156.1875,
                   'N':114.1038,
                   'D':115.0886,
                   'C':103.1388,
                   'Q':128.1307,
                   'E':129.1155,
                   'G':57.0519,
                   'H':137.1411,
                   'I':113.1594,
                   'L':113.1594,
                   'K':128.1741,
                   'M':131.1926,
                   'F':147.1766,
                   'P':97.1167,
                   'S':87.0782,
                   'T':101.1051,
                   'W':186.2132,
                   'Y':163.1760,
                   'V':99.1326}
    mySeq = subSeq
    analysis = ProteinAnalysis(mySeq)
    listofaminoacids.append(analysis.count_amino_acids())

    for i in listofaminoacids:
        for value in i:
            for amino in dictOfAmino:
                if value == amino:
                    peptideBonds = peptideBonds + i[value]
                    #print dictOfAmino[value]
                    #print i[value]
                    molecularWeight = molecularWeight + (i[value]*dictOfAmino[value])

    #peptideBonds = peptideBonds - 1 
    #molecularWeight = molecularWeight - (peptideBonds*waterLoss)
    molecularWeight =  molecularWeight+waterLoss
    return molecularWeight


def getModRes(chain, recordList):
    seqStart = 0
    seqEnd = 0
    tempChainList = []
    modPosition = 0

    modifiedResidues = ""
    
    try:
        tempChainList = chain.split()
        seqStart = int(tempChainList[0])
        seqEnd = int(tempChainList[1])

    except ValueError:
        modifiedResidues = "ERROR"
    for line in recordList:
        lineString = line
        processList = []
        processList = lineString.split('   ')
        
        if(processList[0] == 'FT' and processList[1] == 'MOD_RES'):
            modPosition = int(processList[3])
            if(seqStart <= modPosition <= seqEnd):
                modifiedResidues = modifiedResidues + "Residue:" + str(modPosition) + " "+ "MOD:" + processList[-1] 

    return modifiedResidues


def outputCSV(outCSVName, tempID, tempName, tempMW_mono, tempMW_average, tempMF, tempModRes, tempSubSeq):


    csvFile = open(outCSVName, 'a')
    writer = csv.writer(csvFile)
    writer.writerow([tempID, tempName, tempMW_mono, tempMW_average, tempMF, tempModRes, tempSubSeq])
    csvFile.close()
    

def main():
    totProcessed = 0
    numErrors = 0

    inputName = raw_input("Input file: ")
    outCSVName = raw_input("Output CSV file name (include '.csv'): ")
    
    inputFile = open(inputName, 'r')
    outputFile = open('output.txt', 'a')
    outputErrors = open('errors.txt', 'a')


    fileList = []

    for line in inputFile:
        if(line != '//\n'):
            fileList.append(line)
        else:
            tempSeq = getSequence(fileList)
            tempID = getID(fileList)
            tempChain = getChain(fileList)
            tempSubSeq = ""
            tempMW_mono = 0.0
            tempMW_average = 0.0
            tempMF = ""
            tempModRes = ""

            print "Processing:" + tempID

            #print tempID, ",", tempChain
            for i in tempChain:
                #getSubSeq(i, tempSeq)
                #outString = tempID + "   #     " + i
                
                tempName = getName(i)
                tempModRes = getModRes(i, fileList)
                tempSubSeq = getSubSeq(i, tempSeq)
                tempMF = getMF(tempSubSeq)
                tempMW_mono = getMW_mono(tempSubSeq)
                tempMW_average = getMW_average(tempSubSeq)
                outString = tempID + "    "  + "    " + str(tempMW_mono) + "    " + tempMF
                if(tempSubSeq == "ERROR"):
                    numErrors = numErrors + 1
                    tempName = "Subsequence Error"
                    tempMW_mono = 0
                    tempMW_average = 0
                    tempMF = ""
                    outputErrors.write(tempID)
                    outputErrors.write("\n")
                    
                outputCSV(outCSVName,tempID, tempName, tempMW_mono, tempMW_average, tempMF, tempModRes, tempSubSeq)
                totProcessed = totProcessed + 1
                outputFile.write(outString)
                outputFile.write("\n")

            #print tempID, "", tempName
            #print tempSeq
            
            fileList = []

        
            
    print "Total sequences processed:" + str(totProcessed)
    print "Number of errors:"  + str(numErrors)            

    inputFile.close()
    outputFile.close()
    outputErrors.close()


main()







    
