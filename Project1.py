from Bio.Blast import NCBIXML
from collections import defaultdict

class CovidVariant:
   def __init__(self, nm, genome, s, other):
      self.name = nm
      self.whole = genome
      self.sprotein = s
      self.other = other

def readFile(filename):
   f = open(filename)
   genome = filename.split(".",1)[0]
   curr_seq = []
   for line in f:
      line = line.strip()
      if len(line) != 0 and not line.startswith(">"): curr_seq.append(line)
   return ''.join(curr_seq)

def getMutationRate(part, sbjcr,query,refseg, changes, w):
   w.write(f"\nNucleotide mutation of {part} {sbjcr} genome vs {query} genome\n")
   total = 0
   for sample1Nuc, sample2 in changes.items():
      for sample2Nuc, count in sample2.items():
         if sample1Nuc == "N" or sample2Nuc == "N":
            continue
         else:
            if (sample1Nuc in validNucleotide and sample2Nuc in validNucleotide):
               total += count
            percentage = (count/refseg[sample1Nuc])*100
            x,y = sample1Nuc, sample2Nuc
            if x == "-":
               x = "Addition"
               percentage = count/refseg[sample2Nuc]*100
            elif y == "-":
               y = "Deletion"
            w.write(f"{x} -> {y}: {round(percentage,3)}%\n")
   return total

def getFrequency(sbjcr,query,refseg, changes, w2,totalPair,totalSub):
   chi = 0
   expected = 0
   w2.write(f"\nThe frequency of the transition and tranversion nucleotides of {sbjcr} genome vs {query} genome:\n")
   for sample1Nuc, sample2 in changes.items():
      for sample2Nuc, count in sample2.items():
         if (sample1Nuc in ("A","G") and sample2Nuc in ("A","G")) or (sample1Nuc in ("C","T") and sample2Nuc in ("C","T")):
            percentage = (count/refseg[sample1Nuc])*100
            w2.write(f"{sample1Nuc} -> {sample2Nuc}:\n total {sample2Nuc} in {sbjcr}: {refseg[sample2Nuc]}\n substitution count: {count} \n frequency: {round(percentage,3)}%\n")
   w2.write(f"Total Subsitution vs Total Pairs: {totalSub} and {totalPair}\n")
   # w2.write(f"Expected:{expected}\n")
   # w2.write(f"ChiSquare: {round(chiSquare,3)}\n")
   

def compareGenomes(time, variants, filename, query, sbjcr):
   genomeChanges = defaultdict(lambda: defaultdict(lambda:0))
   SProteinChanges = defaultdict(lambda: defaultdict(lambda:0))
   OtherProteinChanges = defaultdict(lambda: defaultdict(lambda:0))
   sample1 = defaultdict(lambda:0)
   totalPair = 0
   totalSub = 0
   chi = 0
   translatedQuery = ""
   translatedSubject = ""
   silent = 0
   nonsilent = 0    
   
   f = open(filename,"r")
   item = next(NCBIXML.parse(f))
   w4.write(f"{query} and {sbjcr} C-T silent mutation vs nonsilent:\n")
   for alignment in item.alignments:
      for hsp in alignment.hsps:
         if hsp.expect < 0.01:
            for num, nuc in enumerate(hsp.query):
               if (hsp.query[num] in validNucleotide or hsp.sbjct[num] in validNucleotide):
                  if (hsp.query[num] in validNucleotide): translatedQuery += transcript[hsp.query[num]]
                  if (hsp.sbjct[num] in validNucleotide): translatedSubject += transcript[hsp.sbjct[num]]
               else :
                  if (hsp.query[num] not in validNucleotide): translatedQuery += hsp.query[num]
                  if (hsp.sbjct[num] not in validNucleotide): translatedSubject += hsp.sbjct[num]
               totalPair +=1
               sample1[nuc] +=1
               if hsp.query[num] != hsp.sbjct[num]:
                  if num in range(21563, 25385) or num in range(266, 21556) or num in range(25393, 26221) or num in range(26245, 26473) or num in range(26523, 27192) or num in range(27202, 27388) or num in range(27394, 27888) or num in range(27894, 28260) or num in range(28274, 29534) or num in range(29558, 29675):
                     if (hsp.query[num] == "C" and hsp.sbjct[num] == "T"):
                        location = (num-21563)%3
                        if location == 0:
                           queryCodon = hsp.query[num] + hsp.query[num+1] + hsp.query[num+2]
                           sbjctCodon = hsp.sbjct[num] + hsp.sbjct[num+1] + hsp.sbjct[num+2]
                        elif location == 1:
                           queryCodon = hsp.query[num-1] + hsp.query[num] + hsp.query[num+1]
                           sbjctCodon =  hsp.sbjct[num-1] +  hsp.sbjct[num] +  hsp.sbjct[num+1]
                        elif location == 1:
                           queryCodon = hsp.query[num-2] + hsp.query[num-1] + hsp.query[num]
                           sbjctCodon =  hsp.sbjct[num-2] +  hsp.sbjct[num-1] +  hsp.sbjct[num]
                        queryCodon = queryCodon.replace("T","U")   
                        sbjctCodon = sbjctCodon.replace("T","U")
                        if proteinTranslate[queryCodon] == proteinTranslate[sbjctCodon]:
                           w4.write(f"Silent C-U mutation {queryCodon} to {sbjctCodon} ({proteinTranslate[queryCodon]} to {proteinTranslate[sbjctCodon]}) at location {num} on Refseg\n")
                           silent +=1
                        else:
                           w4.write(f"Non-Silent C-U mutation {queryCodon} to {sbjctCodon} ({proteinTranslate[queryCodon]} to {proteinTranslate[sbjctCodon]}) at location {num} on Refseg\n")
                           nonsilent+=1
                  genomeChanges[hsp.query[num]][hsp.sbjct[num]] += 1
                  if num in range(21563, 25385):
                     SProteinChanges[hsp.query[num]][hsp.sbjct[num]] += 1
                  else:
                     OtherProteinChanges[hsp.query[num]][hsp.sbjct[num]] += 1     
   w4.write(f" Total silent C-U: {silent}\n Total nonsilent C-U:{nonsilent}\n\n")
   
   totalSub = getMutationRate("whole",sbjcr,query,sample1,genomeChanges, w)
   getMutationRate("S protein of",sbjcr,query,sample1,SProteinChanges, w)
   getMutationRate("other protein of",sbjcr,query,sample1,OtherProteinChanges, w)
   chi = getFrequency(sbjcr,query, sample1, genomeChanges,w2, totalPair,totalSub)
   
   
   if time == 0: 
      w3.write(f"\n{query} RNA:\n{translatedQuery} \n")  
   w3.write(f"\n{sbjcr} RNA: \n{translatedSubject} \n")
   
   variants[sbjcr].append(CovidVariant(sbjcr,genomeChanges, SProteinChanges, OtherProteinChanges))

#===============================================================================#
genomes = {}
f = open("genome.txt")
w = open("result.txt", "w")
w2 = open("frequency.txt", "w")
w3 = open("transcription.txt", "w")
w4 = open("translation.txt", "w")
variants = defaultdict(list)

validNucleotide = ("A","C","G","T")
transcript = {"A":"A", "T":"U", "G":"G", "C":"C"}

proteinTranslate = {
   "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
   "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R", "AGA":"R", "AGG":"R",
   "AAU":"N", "AAC":"N",
   "GAU":"D", "GAC": "D",
   "AAU":"B", "AAC":"B", "GAU":"B", "GAC":"B",
   "UGU":"C", "UGC":"C",
   "CAA":"Q", "CAG":"Q",
   "GAA":"E", "GAG":"E",
   "CAA":"Z", "CAG":"Z", "GAA":"Z", "GAG":"Z",
   "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G",
   "CAU":"H", "CAC":"H",
   "AUG": " M",
   "AUU":"I", "AUC":"I", "AUA":"I",
   "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L", "UUA":"L", "UUG":"L",
   "AAA":"K", "AAG":"K",
   "UUU":"F", "UUC":"F",
   "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
   "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S", "AGU":"S", "AGC":"S",
   "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
   "UGG":"W",
   "UAU": "Y", "UAC": "Y",
   "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
   "UAA" :"STOP", "UGA":"STOP", "UAG":"STOP"
   }
time = 0

for line in f:
   elements = line.strip()
   if len(elements) != 0:
      elements = elements.split()
      filename, query, sbjcr , g1, g2 = elements[0], elements[1], elements[2], elements[3], elements[4]
      compareGenomes(time, variants, filename, query, sbjcr)
      time +=1
   