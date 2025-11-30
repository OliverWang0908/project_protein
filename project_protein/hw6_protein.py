"""
15-110 F25 Hw6 - Protein Sequencing Project
Name: Guanyang Wang
AndrewID: guanyanw
"""

import hw6_protein_tests as test

project = "Protein" # don't edit this

### WEEK 1 ###

'''
readFile(filename)
#1 [Check6-1]
Parameters: str
Returns: str
'''
def readFile(filename):
    f = open(filename, "r")
    text = f.read()
    text = text.replace("\n", "")
    return text



'''
dnaToRna(dna, startIndex)
#2 [Check6-1]
Parameters: str ; int
Returns: list of strs
'''
def dnaToRna(dna, startIndex):
    rna = dna.replace("T", "U")
    rna = rna[startIndex:]
    ls = []
    for i in range(0,len(rna),3):
        if rna[i:i+3] in ["UAA","UAG", "UGA"]:
            ls.append(rna[i:i+3])
            break
        ls.append(rna[i:i+3])
    return ls


'''
makeCodonDictionary(filename)
#3 [Check6-1]
Parameters: str
Returns: dict mapping strs to strs
'''
def makeCodonDictionary(filename):
    import json
    f = open(filename, "r")
    codon_table = json.load(f)
    dic = {}
    for aa in codon_table:
        for codon in codon_table[aa]:
            codon = codon.replace("T", "U")
            dic[codon] = aa
    return dic


'''
generateProtein(codons, codonD)
#4 [Check6-1]
Parameters: list of strs ; dict mapping strs to strs
Returns: list of strs
'''
def generateProtein(codons, codonD):
    ls = []
    for codon in codons:
        if codonD[codon] == "Stop":
            ls.append("Stop")
            break
        ls.append(codonD[codon])
    ls[0] = "Start"
    return ls


'''
synthesizeProteins(dnaFilename, codonFilename)
#5 [Check6-1]
Parameters: str ; str
Returns: 2D list of strs
'''
def synthesizeProteins(dnaFilename, codonFilename):
    dna = readFile(dnaFilename)
    ls = []
    i = 0
    while i <= len(dna)-2:
        if dna[i:i+3] == "ATG":
            codons = dnaToRna(dna,i)
            codonD = makeCodonDictionary(codonFilename)
            protein_lst = generateProtein(codons,codonD)
            ls.append(protein_lst)
            i = i + 3*len(protein_lst)
        else:
            i = i + 1
       
    return ls


def runWeek1():
    print("Human DNA")
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    print("Elephant DNA")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")


### WEEK 2 ###

'''
commonProteins(proteinList1, proteinList2)
#1 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs
Returns: 2D list of strs
'''
def commonProteins(proteinList1, proteinList2):
    ls = []
    for i in proteinList1:
        if i in proteinList2:
            if i not in ls:
                ls.append(i)

    return ls


'''
combineProteins(proteinList)
#2 [Check6-2]
Parameters: 2D list of strs
Returns: list of strs
'''
def combineProteins(proteinList):
    ls = []
    for i in proteinList:
        for j in i:
            ls.append(j)
    
    return ls


'''
aminoAcidDictionary(aaList)
#3 [Check6-2]
Parameters: list of strs
Returns: dict mapping strs to ints
'''
def aminoAcidDictionary(aaList):
    dic = {}
    for i in aaList:
        if i not in dic:
            dic[i] = 1
        else:
            dic[i] = dic[i]+1

    return dic


'''
findAminoAcidDifferences(proteinList1, proteinList2, cutoff)
#4 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs ; float
Returns: 2D list of values
'''

def findAminoAcidDifferences(proteinList1, proteinList2, cutoff):
    ls = []
    ls1 = combineProteins(proteinList1)
    ls2 = combineProteins(proteinList2)
    
    ls1_n = list(set(ls1))
    ls2_n = list(set(ls2))
    ls3 = list(set(ls1_n + ls2_n))
    print(ls3)
    for i in ls3:
        if i not in ["Start", "Stop"]:
            fr1 = ls1.count(i)/len(ls1)
            fr2 = ls2.count(i)/len(ls2)
            if abs(fr1 - fr2) > cutoff:
                ls.append([i, fr1, fr2])
    
    

    return ls
    


'''
displayTextResults(commonalities, differences)
#5 [Check6-2]
Parameters: 2D list of strs ; 2D list of values
Returns: None
'''
def displayTextResults(commonalities, differences):
    common = []
    for i in commonalities:
        if len(i) > 2:
            common.append(i)
    print("The following proteins occured in both DNA sequences:")
    for i in common:
        st = ""
        for j in i:
            if j not in ["Start", "Stop"]:
                st = st + j + "-"
        st = st[:-1]
        print(st)
    print("The following amino acids occurred at very different rates in the two DNA sequences:")
    for i in differences:
        print("{}: {}% in Seq1, {}% in Seq2".format(i[0], round(i[1]*100,2), round(i[2]*100,2)))


def runWeek2():
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")

    commonalities = commonProteins(humanProteins, elephantProteins)
    differences = findAminoAcidDifferences(humanProteins, elephantProteins, 0.005)
    displayTextResults(commonalities, differences)


### WEEK 3 ###

'''
makeAminoAcidLabels(proteinList1, proteinList2)
#2 [Hw6]
Parameters: 2D list of strs ; 2D list of strs
Returns: list of strs
'''
def makeAminoAcidLabels(proteinList1, proteinList2):
    
    ls1 = combineProteins(proteinList1)
    ls2 = combineProteins(proteinList2)
    
    ls1 = list(set(ls1))
    ls2 = list(set(ls2))
    ls = list(set(ls1 + ls2))
    ls.sort()
    return ls


'''
setupChartData(labels, proteinList)
#3 [Hw6]
Parameters: list of strs ; 2D list of strs
Returns: list of floats
'''
def setupChartData(labels, proteinList):
    ls1 = combineProteins(proteinList)
    ls = []
    for i in labels:
        fr = ls1.count(i)/len(ls1)
        ls.append(fr)
    
    return ls


'''
createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None)
#4 [Hw6] & #5 [Hw6]
Parameters: list of strs ; list of floats ; str ; list of floats ; str ; [optional] list of strs
Returns: None
'''
def createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None):
    import matplotlib.pyplot as plt
    import numpy as np
    
    # Create x positions
    x = np.arange(len(xLabels))  # [0,1,2,...]
    width = 0.35                 # width < 0.5 so bars sit side-by-side

    # Create the figure
    fig, ax = plt.subplots()

    # If an edgeList (colors for outlines) is provided, use it for bar edges
    if edgeList is not None:
        rects1 = ax.bar(x - width/2, freqList1, width, label=label1, edgecolor=edgeList, linewidth=2)
        rects2 = ax.bar(x + width/2, freqList2, width, label=label2, edgecolor=edgeList, linewidth=2)
    else:
        # First bar group (left shift by width/2)
        rects1 = ax.bar(x - width/2, freqList1, width, label=label1)

        # Second bar group (right shift by width/2)
        rects2 = ax.bar(x + width/2, freqList2, width, label=label2)

    # X-axis labels
    ax.set_xticks(x)
    ax.set_xticklabels(xLabels)

    # Add legend
    ax.legend()

    # OPTIONAL: print values above bars (nice for grading)
    # ax.bar_label(rects1, padding=3)
    # ax.bar_label(rects2, padding=3)

    plt.show()


'''
makeEdgeList(labels, biggestDiffs)
#5 [Hw6]
Parameters: list of strs ; 2D list of values
Returns: list of strs
'''
def makeEdgeList(labels, biggestDiffs):
    ls = []
    for i in biggestDiffs:
        ls.append(i[0])
    
    color_ls = []
    for i in labels:
        if i in ls:
            color_ls.append("black")
        else:
            color_ls.append("white")
    
    return color_ls


'''
runFullProgram()
#6 [Hw6]
Parameters: no parameters
Returns: None
'''
def runFullProgram():
    # Load and synthesize proteins from the two p53 DNA files
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")

    # Generate textual comparison (common proteins and amino-acid differences)
    commonalities = commonProteins(humanProteins, elephantProteins)
    differences = findAminoAcidDifferences(humanProteins, elephantProteins, 0.005)
    displayTextResults(commonalities, differences)

    # Prepare chart data and highlight sufficiently-different amino acids
    labels = makeAminoAcidLabels(humanProteins, elephantProteins)
    freqHuman = setupChartData(labels, humanProteins)
    freqElephant = setupChartData(labels, elephantProteins)
    edgeList = makeEdgeList(labels, differences)

    # Create bar chart comparing the two genes
    createChart(labels, freqHuman, "Human p53", freqElephant, "Elephant p53", edgeList)

    return


### RUN CODE ###

# This code runs the test cases to check your work
if __name__ == "__main__":
    print("\n" + "#"*15 + " WEEK 1 TESTS " +  "#" * 16 + "\n")
    test.week1Tests()
    print("\n" + "#"*15 + " WEEK 1 OUTPUT " + "#" * 15 + "\n")
    runWeek1()

    ## Uncomment these for Week 2 ##
  
    print("\n" + "#"*15 + " WEEK 2 TESTS " +  "#" * 16 + "\n")
    test.week2Tests()
    print("\n" + "#"*15 + " WEEK 2 OUTPUT " + "#" * 15 + "\n")
    runWeek2()
   

    ## Uncomment these for Week 3 ##
    
    print("\n" + "#"*15 + " WEEK 3 TESTS " +  "#" * 16 + "\n")
    test.week3Tests()
    print("\n" + "#"*15 + " WEEK 3 OUTPUT " + "#" * 15 + "\n")
    runFullProgram()
    
