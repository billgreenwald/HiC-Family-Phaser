
# coding: utf-8

# This script will:
# 
#     fix switch errors
#     
#     phase unphased variants that are distinguishable from parental haplotypes 
#     
#     unphase phased variants that are switch errors if people are locked
#     
#     remove sites from file containing multiallilic SNVs (can be turned off)

# # Input Parameters

# In[ ]:


import argparse
import os
from __future__ import print_function
import sys

parser=argparse.ArgumentParser()
parser._optionals.title = "Flag Arguments"
parser.add_argument('-i', help='Input vcf.  Required', required=True)
parser.add_argument('-o', help='Output file name.  Required', required=True)
parser.add_argument('-s', help='Sample to phase.  Required', required=True)
parser.add_argument('-p', help='Pedigree File.  Required', required=True)
parser.add_argument('-l', help='Samples where phase errors should not be fixed, provided as a comma separated list.  Optional', required=False, default="")
args = vars(parser.parse_args())


# In[ ]:


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


# In[ ]:


if not os.path.exists(args['i']):
    eprint("Input vcf does not exist.  Exiting")
    exit()
if not os.path.exists(args['p']):
    eprint("Pedigree File does not exist.  Exiting")
    exit()
lc=0
for line in open(args['p']):
    lc+=1
    if len(line.split())!=3:
        eprint("WARNING: Line "+str(lc)+" from pedgiree file does not have exactly 3 fields.")
if len(args['l'])!= and and "," not in args['l']:
    eprint("WARNING: -l argument provided, but no commas were found.  If only one sample is provided, there is no error.  If intended use was to input multiple samples, there is an error.")


# In[45]:


vcfFile=args['i']
outFile=args['o']
person=args['s']
family={line.split()[0]:(line.split()[1],line.split()[2]) for line in open(args['p'])} #c \t p \t m
LOCKED_PEOPLE=set(args['l'].split(","))
DEBUG=False


# # Functions

# In[46]:


def fastFindMins(nonab):
    wind=1000
    wind2=200
    t=[]
    t2=[]
    pre=np.mean(nonab[0:wind])
    post=np.mean(nonab[wind+1:2*wind+1])
    for i in range(wind+1,len(nonab)-wind):
        pre=pre-float(nonab[i-wind-1])/wind+float(nonab[i-1])/wind
        post=post-float(nonab[i])/wind+float(nonab[wind+i])/wind
        if pre>nonab[i] and post>nonab[i]:
            t.append(nonab[i])
            t2.append(i) 

            
    final={}
    for i in range(wind2,len(t)-wind2):
        if min(t[i-wind2:i])>=t[i] and min(t[i+1:i+wind2+1])>=t[i]:
            if t[i] not in final:
                final[t[i]]=[t2[i]]
            else:
                final[t[i]].append(t2[i])

    return[final[key][len(final[key])/2] for key in final]

def fastFindMax(nonab):
    wind=1000
    wind2=200
    t=[]
    t2=[]
    pre=np.mean(nonab[0:wind])
    post=np.mean(nonab[wind+1:2*wind+1])
    for i in range(wind+1,len(nonab)-wind):
        pre=pre-float(nonab[i-wind-1])/wind+float(nonab[i-1])/wind
        post=post-float(nonab[i])/wind+float(nonab[wind+i])/wind
        if pre<nonab[i] and post<nonab[i]:
            t.append(nonab[i])
            t2.append(i) 

            
    final={}
    for i in range(wind2,len(t)-wind2):
        if max(t[i-wind2:i])<=t[i] and max(t[i+1:i+wind2+1])<=t[i]:
            if t[i] not in final:
                final[t[i]]=[t2[i]]
            else:
                final[t[i]].append(t2[i])

    return[final[key][len(final[key])/2] for key in final]


# In[47]:


def getStartingHaps(chrom,blocks,chromInd,startingParentForEachChrom,startingParentHapForEachChrom):
    chromblocks=[x for x in blocks if x[3]==chrom]
    p1hap=2
    p2hap=2
    for x in chromblocks:
        if x[1]==0 and p1hap==2:
            if x[2]=='min':
                p1hap=1
            else:
                p1hap=0
        if x[1]==1 and p2hap==2:
            if x[2]=='min':
                p2hap=1
            else:
                p2hap=0
    if p1hap==2:
        p1hap=startingParentHapForEachChrom[chromInd][0]
    if p2hap==2:
        p2hap=startingParentHapForEachChrom[chromInd][1]
        
    return p1hap,p2hap


# In[48]:


def phaseUnphasedVar(childGenotype,leftMatchingParentHaplotype,rightMatchingParentHaplotype):
    cL=childGenotype[0]
    cR=childGenotype[2]
    
    if cL==leftMatchingParentHaplotype and cR==rightMatchingParentHaplotype:
        return cL+"|"+cR
    elif cR==leftMatchingParentHaplotype and cL==rightMatchingParentHaplotype:
        return cR+"|"+cL
    else:
        return childGenotype


# In[49]:


def phaseChildToParents(unph1,unph2,phase1,phase2):
    if unph1==phase1 and unph2 == phase2:
        return unph1+"|"+unph2
    elif unph1==phase2 and unph2==phase1:
        return unph2+"|"+unph1
    else:
        return unph1+"/"+unph2
    
def phaseParentToChild(unph1,unph2,cL):
    if unph1==cL:
        return unph1+"|"+unph2
    elif unph2==cL:
        return unph2+"|"+unph1
    else:
        return unph1+"/"+unph2
    
def phaseChildToSingleParent(unph1,unph2,pL,pHap):
    if unph1==pL:
        if pHap==0:
            return unph1+"|"+unph2
        else:
            return unph2+"|"+unph1
    elif unph2==pL:
        if pHap==0:
            return unph2+"|"+unph1
        else:
            return unph1+"|"+unph2
    else:
        return unph1+"/"+unph2
    


# In[50]:


def phaseUnphasedVar2(c,p1,p2,p1L,p2L,p1Hap,p2Hap):
    if c[1]=="/":
        if p1[1]=="/":
            if p2[1]=="/":
                return "^  ","^  ","^  " #cant phase any
            else: #phase c via p2, then phase p1 via c
                c=phaseChildToSingleParent(c[0],c[2],p2L,p2Hap)
                return phaseParentToChild(p1[0],p1[2],c[p1Hap*2]),"^  ",c
        else:
            if p2[1]=="/": #phase c via p1, then phase p2 via c
                c=phaseChildToSingleParent(c[0],c[2],p1L,p1Hap)
                return "^  ",phaseParentToChild(p2[0],p2[2],c[p2Hap*2]),c
            else: #phase child via p1 and p2
                return "^  ","^  ",phaseChildToParents(c[0],c[2],p1L,p2L)
    else:
        if p1[1]=="/":
            if p2[1]=="/": #phase p1 via c and then p2 via c
                return phaseParentToChild(p1[0],p1[2],c[p1Hap*2]),phaseParentToChild(p2[0],p2[2],c[p2Hap*2]),"^  "
            else: #phase p1 via c
                return phaseParentToChild(p1[0],p1[2],c[p1Hap*2]),"^  ","^  "
        else: #phase p2 via c
            return "^  ",phaseParentToChild(p2[0],p2[2],c[p2Hap*2]),"^  "


# In[51]:


def fixSwitchErrors2(cL,cR,p1L,p1R,p2L,p2R):
    
    #cL and cR match p1L and p2L
    if cL==p1L and cR==p2L:
        return 0 #no switch
    
    #c is hom
    if cL==cR:

        #cL matches p1L, cR does not match p2L, cR matches p2R
        if cL==p1L and cR != p2L and cR==p2R:
            return 1 #p2 switch
        #cR does not match either p2
        if cR!=p2L and cR!=p2R:
            return -1 #genotype error
        #cL does not match p1L or p1R
        if cL!=p1L and cL!=p1R:
            return -1 #genotype error
        #cL does not match p1L, cL matches p1R, cR matches p2L
        if cL!=p1L and cL==p1R and cR==p2L:
            return 2 #p1 switch
        #cL does not match p1L, cL does match p1R, cR does not match p2L, cR does match p2R
        if cL!=p1L and cL==p1R and cR!=p2L and cR==p2R:
            return 3 # double parental switch
        

    #c is het
    else:
        
        #cL doesnt match either p1, cR matches both p2
        if cL!=p1L and cL!=p1R and cR==p2L and cR==p2R:
            return -1 #genotype error
        #cL doesnt match p1L, cR doesnt match p2L, cL matches p2L, cR matches p1L.  
        if cL==p2L and cR==p1L:#probably only need second half
            return 4 #child switch
        #cL matches neither p1, cR doesnt match p2L, cR matches p1L, cL matches p2R
        if cL!=p1L and cL!=p1R and cL!=p2L and cR==p1L and cL==p2R:
            return 5 #child switch, p2 switch
        #cL matches both p1, cR matches neither p2
        if cL==p1R and cL==p1L and cR!=p2L and cR!=p2R:
            return -1 #genotype error
        #cL doesnt match p1L, does match p1R.  cR matches p2L
        if cL!=p1L and cL==p1R and cR==p2L:
            return 2 #p1 switch
        #cL matches p1L, cR doesnt match p2L, matches p2R
        if cL==p1L and cR!=p2L and cR==p2R:
            return 1 #p2 switch
        #cL matches p1L, cR matches neither p2, cR matches p1R but not p1L, cL matches p2L
        if cL==p1L and cR!=p2L and cR!=p2R and cR==p1R and cR!=p1L and cL==p2L:
            return 6 #child switch, p1 switch

    return 7 #shouldnt happen
    
    
    
    


# In[52]:


def fixSwitchError(childGenotype,leftMatchingParentHaplotype,rightMatchingParentHaplotype,leftNonMatchingParentHaplotype,rightNonMatchingParentHaplotype,phaseP,phaseOP):
    if childGenotype[0]==leftMatchingParentHaplotype and childGenotype[2]==rightMatchingParentHaplotype: #matches correctly, return 
        return 0
    else:
        if childGenotype[0]==rightMatchingParentHaplotype and childGenotype[2]==leftMatchingParentHaplotype:  #switch in child, flip it
            return 1
        else:
            if childGenotype[0]!=leftMatchingParentHaplotype and childGenotype[0]==leftNonMatchingParentHaplotype and childGenotype[2]==rightMatchingParentHaplotype: #switch in leftParent only, fix parent
                if phaseP:
                    return 3
                else:
                    return 2
            elif childGenotype[2]!=rightMatchingParentHaplotype and childGenotype[2]==rightNonMatchingParentHaplotype and childGenotype[0]==leftMatchingParentHaplotype: #switch in rightParent only, fix it
                if phaseOP:
                    return 4
                else:
                    return 2
            else: #both switches, more likely genotype error, unphase it
                return 2


# In[53]:


def findInitialHap(lefts,rights,startingParentForEachChrom,startingParentHapForEachChrom):
    topHapsR=[]
    largestValR=0

    topHapsL=[]
    largestValL=0
    for i in range(4):
        #get all largest values incase there are multiple (very rare case but good to catch)
        if lefts[i]>largestValL:
            topHapsL=[i]
            largestValL=lefts[i]
        elif lefts[i]==largestValL:
            topHapsL.append(i)

        if rights[i]>largestValR:
            topHapsR=[i]
            largestValR=rights[i]
        elif rights[i]==largestValR:
            topHapsR.append(i)

    #detect the chromosomes
    if len(topHapsL)==1 and len(topHapsR)==1: #should be most of the time
        cL=topHapsL[0]
        cR=topHapsR[0]
        if (cL>=2 and cR>=2) or (cL<2 and cR<2):
            eprint("WARNING: "+chrom+ " has the same top haplotype.")

    elif len(topHapsL)>1 and len(topHapsR)==1: #resolve the chromosomes
        cR=topHapsR[0]
        if cR>=2:
            topHapsL=[x for x in topHapsL if x <2]
            if len(lefts)>2:
                eprint("WARNING: "+chrom+" child left haplotype matches both haplotypes on the same parent equally well.  Picking randomly")
            cL=topHapsL[0]
        else:
            lefts=[x for x in lefts if x >=2]
            if len(lefts)>2:
                eprint("WARNING: "+chrom+" child left haplotype matches both haplotypes on the same parent equally well.  Picking randomly")
            cL=topHapsL[0]

    elif len(topHapsL)>1 and len(topHapsR)==1: #resolve the chromosomes
        cL=topHapsL[0]
        if cL>=2:
            topHapsR=[x for x in topHapsR if x <2]
            if len(topHapsR)>2:
                eprint("WARNING: "+chrom+" child right haplotypes matches bothon the same parent equally well.  Picking randomly")
            cR=topHapsR[0]
        else:
            topHapsR=[x for x in topHapsR if x >=2]
            if len(topHapsR)>2:
                eprint("WARNING: "+chrom+" child right haplotypes matches bothon the same parent equally well.  Picking randomly")
            cR=topHapsR[0]

    else: #very rare, should basically never happen
        #multiple cases here, shouldn't happen, will code later if we choose to release the tool
        eprint("WARNING: indistinguishable chromosome.  Exiting")
        Tracer()()

    if cL==0 and cR==2:
        startingParentForEachChrom.append((0,1))
        startingParentHapForEachChrom.append((0,0))

    elif cL==1 and cR==2:
        startingParentForEachChrom.append((0,1))
        startingParentHapForEachChrom.append((1,0))

    elif cL==0 and cR==3:
        startingParentForEachChrom.append((0,1))
        startingParentHapForEachChrom.append((0,1))

    elif cL==1 and cR==3:
        startingParentForEachChrom.append((0,1))
        startingParentHapForEachChrom.append((1,1))

    elif cR==0 and cL==2:
        startingParentForEachChrom.append((1,0))
        startingParentHapForEachChrom.append((0,0))

    elif cR==1 and cL==2:
        startingParentForEachChrom.append((1,0))
        startingParentHapForEachChrom.append((0,1))

    elif cR==0 and cL==3:
        startingParentForEachChrom.append((1,0))
        startingParentHapForEachChrom.append((1,0))

    elif cR==1 and cL==3:
        startingParentForEachChrom.append((1,0))
        startingParentHapForEachChrom.append((1,1))                                 
    else:
        eprint("LOGIC ERROR.  EXITING")
        exit()
    return


# ### identify starting haplotypes

# In[55]:


startingParentForEachChrom=[] #list has chromosome number of entries, each tuples of (parentMatchingChildHap1,parentMatchingChildHap2)
startingParentHapForEachChrom=[] #list has chromosome number of entries, each tuples of (parentalHaplotype,parentalHaplotype)

lefts=[0,0,0,0]
rights=[0,0,0,0]  #both of these are p1L, p1R, p2L, p2R

chromOrder={}


previousBP=-1
switch=True
howManyChrs=-1
first=True

with open(vcfFile) as f:
    for line in tqdm(f):
        if line[0]=="#":
            if line[1]!="#": #get header so we know where to look for each sample
                line=line.strip().split()
                sampleToColumn={line[i]:i for i in range(len(line))}
            continue
        line=line.split()
        if first:
            chrom=line[0]
            first=False

        if chrom!=line[0]:
            chromOrder[chrom]=line[0]
            
            findInitialHap(lefts,rights,startingParentForEachChrom,startingParentHapForEachChrom)
            lefts=[0,0,0,0]
            rights=[0,0,0,0]
            
            chrom=line[0]
                    

        if line[sampleToColumn[person]][1]=="/" or line[sampleToColumn[family[person][0]]][1]=="/" or line[sampleToColumn[family[person][0]]][1]=="/":
            continue

        cL=line[sampleToColumn[person]][0]
        cR=line[sampleToColumn[person]][2]

        p1L=line[sampleToColumn[family[person][0]]][0]
        p1R=line[sampleToColumn[family[person][0]]][2]

        p2L=line[sampleToColumn[family[person][1]]][0]
        p2R=line[sampleToColumn[family[person][1]]][2]

        if cL==p1L:
            lefts[0]+=1
        if cL==p1R:
            lefts[1]+=1
        if cL==p2L:
            lefts[2]+=1
        if cL==p2R:
            lefts[3]+=1

        if cR==p1L:
            rights[0]+=1
        if cR==p1R:
            rights[1]+=1
        if cR==p2L:
            rights[2]+=1
        if cR==p2R:
            rights[3]+=1
            
findInitialHap(lefts,rights,startingParentForEachChrom,startingParentHapForEachChrom)    
    


# ### identify crossovers

# In[56]:


#read thru the vcfFile and stop at each chromsome boundary to chunk

blocks={x:[] for x in family} #holds the line numbers of the haplotype block splits
unphasedAtParent={x:set() for x in family}

for matchingHap in [0,1]: #0 for left, 1 for right.  Will need to loop this
    chromListInd=0
    matchingParent=startingParentForEachChrom[chromListInd][matchingHap] #0 for first parent, 1 for second
    matchingParentHap=0  #0 for left, 1 for right
    lineNumber=0

    chrom='chr1a'

    Match=[] #holds the match scores
    matchIndexes=[] #hold the indexes of each informative SNP
    dists=[]
    with open(vcfFile) as f:
        for line in tqdm(f):
            lineNumber+=1
            if line[0]=="#":
                if line[1]!="#": #get header so we know where to look for each sample
                    line=line.strip().split()
                    sampleToColumn={line[i]:i for i in range(len(line))}
                continue

            #find the peaks in the match score

            line=line.split()  

            if line[0]!=chrom: #when we advance chromsosomes, stop and process the previous chrom
                half2=sum(Match)
                half1=0
                dists=[half1-half2]
                for x in Match:
                    half1+=x
                    half2-=x
                    dists.append(half1-half2)
                            
                #add the extreme points and the parents that had them to the blocks list
                blocks[person].extend([(matchIndexes[x],matchingHap,'max',chrom) for x in fastFindMax(dists)])
                blocks[person].extend([(matchIndexes[x],matchingHap,'min',chrom) for x in fastFindMins(dists)])
                                

                #reset the lists and keep going
                Match=[]
                matchIndexes=[]
                dists=[]
                chromListInd+=1
                matchingParent=startingParentForEachChrom[chromListInd][matchingHap]
                chrom=chromOrder[chrom]


                #at the end, blocks has all the split points for a single haplotype in a single person in a single trio
            
            
            
            #check to make sure variant is phased in both parents:
            if line[sampleToColumn[family[person][matchingParent]]][1]!="|":
                unphasedAtParent[person].add(lineNumber)
                continue
                
            if line[sampleToColumn[person]][1]!="|":
                continue

            #check for match score.  if haplotypes are different, get +1 for match, -1 for mismatch, else get a 0
            p1=line[sampleToColumn[family[person][matchingParent]]].split("|")[matchingParentHap] #parental hap1
            p2=line[sampleToColumn[family[person][matchingParent]]].split("|")[(matchingParentHap+1)%2] #parental hap2
            c1=line[sampleToColumn[person]].split("|")[matchingHap] #child hap
            if p1==p2:
                pass
            elif p1==c1:
                Match.append(1)
                matchIndexes.append(lineNumber)
            elif p2==c1:
                Match.append(-1)
                matchIndexes.append(lineNumber)

    #need to process last line in file
    half2=sum(Match)
    half1=0
    dists=[half1-half2]
    for x in Match:
        half1+=x
        half2-=x
        dists.append(half1-half2)

    #add the extreme points and the parents that had them to the blocks list
    blocks[person].extend([(matchIndexes[x],matchingHap,'max',chrom) for x in fastFindMax(dists)])
    blocks[person].extend([(matchIndexes[x],matchingHap,'min',chrom) for x in fastFindMins(dists)])

    
t=[chromOrder[x] for x in chromOrder]
s=[x for x in chromOrder if x not in t][0]
i=0
sortOrder={}
while s in chromOrder:
    sortOrder[s]=i
    s=chromOrder[s]
    i+=1
sortOrder[s]=i
for x in blocks:
    blocks[x]=sorted(blocks[x],key=lambda x: (sortOrder[x[3]],x[0]))


# ### divide chromsome into "recombinant blocks"

# Take all crosover points in both chromosomes, and break VCF at each crossover point (ie P1 has 2 crossovers, P2 has 1 crossover, create 3 break points)

# In[ ]:


gentypeDict={}
for c in ["0|0","0|1","1|0","1|1","0/1","1/0"]:
    for p1 in ["0|0","0|1","1|0","1|1","0/1","1/0"]:
        for p2 in ["0|0","0|1","1|0","1|1","0/1","1/0"]:
            gentypeDict[c+","+p1+","+p2]=[]


# In[ ]:


concordDict={}


# In[ ]:


FLAG_NUM="a"

#these need to already be populated, all of the same size as "block"

#everytime you hit a block you flip haps on the hap the block break point came from

ReUnphasedVars=[]
Fixes=[]
parentFixes=[]

lineNumber=0
whichBlock=0

chrom='chr1a'
chromInd=0

p1hap,p2hap=getStartingHaps(chrom,blocks[person],chromInd,startingParentForEachChrom,startingParentHapForEachChrom)

if DEBUG:
    print (chrom+":"+str(p1hap)+","+str(p2hap))
with open(vcfFile) as f:
    with open(outFile,"w+") as g:
        for line in tqdm(f):
            lineNumber+=1
            if line[0]=="#":
                g.write(line)
            else:
                line=line.split()
                if line[0]!=chrom:
                    chromInd+=1
                    chrom=line[0]
                    p1hap,p2hap=getStartingHaps(chrom,blocks[person],chromInd,startingParentForEachChrom,startingParentHapForEachChrom)
                    if DEBUG:
                        print (chrom+":"+str(p1hap)+","+str(p2hap))
                    
                if lineNumber==blocks[person][whichBlock][0]:
                    if blocks[person][whichBlock][1]==0:
                        p1hap=(p1hap+1)%2
                    else:
                        p2hap=(p2hap+1)%2
                    if whichBlock!=len(blocks[person])-1:
                        if DEBUG:
                            print (chrom+":"+str(p1hap)+","+str(p2hap))
                        whichBlock+=1
                
                parInd=startingParentForEachChrom[chromInd][0]
                oparInd=startingParentForEachChrom[chromInd][1]        
                c=line[sampleToColumn[person]]
                p1=line[sampleToColumn[family[person][parInd]]] #matching parental hap
                p2=line[sampleToColumn[family[person][oparInd]]]#other parental matching hap
                
                
                concordKey=",".join([str(x) for x in (parInd,oparInd,p1hap,p2hap)])
                if concordKey in concordDict:
                    concordDict[concordKey][lineNumber]=c
                else:
                    concordDict[concordKey]={lineNumber:c}
                
                if line[1]==FLAG_NUM:
                    Tracer()()
                if c!='./.' and p1!='./.' and p2!='./.':

                    #phase unphased things.  To reset to old, just put a continue here if lineNumber in unphasedAtParents
                    p1=line[sampleToColumn[family[person][parInd]]] #matching parental hap
                    p2=line[sampleToColumn[family[person][oparInd]]]#other parental matching hap

                    p1L=p1[p1hap*2] #matching parental hap
                    p1R=p1[((p1hap+1)%2)*2] #nonmatching parental hap
                    p2L=p2[p2hap*2] #other parental matching hap
                    p2R=p2[((p2hap+1)%2)*2] #other parental nonmatching hap

                    if any([x=="2" for x in [p1L,p1R,p2L,p2R,c[0],c[2]]]):
                        parentFixes.append((chrom,int(line[1]),lineNumber,'Multiallelic'))
#                         g.write("\t".join(line)+"\n") #uncomment to write multiallelic lines
                        continue
                    
                    if c[0]==c[2]:
                        c=c[0]+"|"+c[2]
                    if p1[0]==p1[2]:
                        p1=p1[0]+"|"+p1[2]
                    if p2[0]==p2[2]:
                        p2=p2[0]+"|"+p2[2]
                    
                    gentypeDict[c+","+p1+","+p2].append(lineNumber)
                    entered=False
                    if any([x[1]=="/" for x in [c,p1,p2]]):
                        entered=True
                        np1,np2,nc=phaseUnphasedVar2(c,p1,p2,p1L,p2L,startingParentForEachChrom[chromInd][0],startingParentForEachChrom[chromInd][1])
                        Fixes.append(((chrom,int(line[1]),lineNumber,",".join((nc,np1,np2)))))
                        if person not in LOCKED_PEOPLE:
                            if nc!="^  ":
                                c=nc
                        if family[person][parInd] not in LOCKED_PEOPLE: 
                            if np1!="^  ":
                                p1=np1
                        if family[person][oparInd] not in LOCKED_PEOPLE:
                            if np2!="^  ":
                                p2=np2
                    
                    
                    if nc=="^  " and np1=="^  " and np2=="^  " and entered:
                        pass
                    elif True:
                        p1L=p1[p1hap*2] #matching parental hap
                        p1R=p1[((p1hap+1)%2)*2] #nonmatching parental hap
                        p2L=p2[p2hap*2] #other parental matching hap
                        p2R=p2[((p2hap+1)%2)*2] #other parental nonmatching hap

                        
                        cL=c[0]
                        cR=c[2]

                        if line[1]==FLAG_NUM:
                            Tracer()()
                        se=fixSwitchErrors2(cL,cR,p1L,p1R,p2L,p2R)
                        if line[1]==FLAG_NUM:
                            Tracer()()
                        if se==0:
                            pass
                        elif se==1:
                            parentFixes.append((chrom,int(line[1]),lineNumber,'P2S'))
                            if family[person][oparInd] not in LOCKED_PEOPLE:
                                p2=p2[2]+"|"+p2[0]
                            else:
                                p2=p2[2]+"/"+p2[0]
                        elif se==2:
                            parentFixes.append((chrom,int(line[1]),lineNumber,'P1S'))
                            if family[person][parInd] not in LOCKED_PEOPLE:
                                p1=p1[2]+"|"+p1[0]
                            else:
                                p1=p1[2]+"/"+p1[0]
                        elif se==-1:
                            parentFixes.append((chrom,int(line[1]),lineNumber,'GE'))
                            if person not in LOCKED_PEOPLE:
                                c=c[0]+"/"+c[2]
                        elif se==3:
                            parentFixes.append((chrom,int(line[1]),lineNumber,'P1S, P2S'))
                            if family[person][oparInd] not in LOCKED_PEOPLE and family[person][parInd] not in LOCKED_PEOPLE:
                                p2=p2[2]+"|"+p2[0]
                                p1=p1[2]+"|"+p1[0]
                            else:
                                p2=p2[2]+"/"+p2[0]
                                p1=p1[2]+"/"+p1[0]
                        elif se==4:
                            parentFixes.append((chrom,int(line[1]),lineNumber,'CS'))
                            if person not in LOCKED_PEOPLE:
                                c=c[2]+"|"+c[0]
                            else:
                                c=c[2]+"/"+c[0]

                        elif se==5:
                            parentFixes.append((chrom,int(line[1]),lineNumber,'CS P2S'))
                            if person not in LOCKED_PEOPLE and family[person][oparInd] not in LOCKED_PEOPLE:
                                c=c[2]+"|"+c[0]
                                p2=p2[2]+"|"+p2[0]
                            else:
                                c=c[2]+"/"+c[0]
                                p2=p2[2]+"/"+p2[0]
                        elif se==6:
                            parentFixes.append((chrom,int(line[1]),lineNumber,'CS P1S'))
                            if person not in LOCKED_PEOPLE and family[person][parInd] not in LOCKED_PEOPLE:
                                c=c[2]+"|"+c[0]
                                p1=p1[2]+"|"+p1[0]
                            else:
                                c=c[2]+"/"+c[0]
                                p1=p1[2]+"/"+p1[0]
                        elif se==7:
                            parentFixes.append((chrom,int(line[1]),lineNumber,'Error'))
                        else:
                            parentFixes.append((chrom,int(line[1]),lineNumber,'Null Return'))
                        
                        if line[1]==FLAG_NUM:
                            Tracer()()
                        
                        line[sampleToColumn[family[person][parInd]]]=p1
                        line[sampleToColumn[family[person][oparInd]]]=p2
                        line[sampleToColumn[person]]=c
                            
                    else:
                        parentFixes.append((chrom,int(line[1]),lineNumber,"Can't Phase"))
                
                    
                    
                    #make paternal on right, maternal on left
                    if startingParentForEachChrom[chromInd][0]==1:
                        line[sampleToColumn[person]]=line[sampleToColumn[person]][2]+"|"+line[sampleToColumn[person]][0]
                g.write("\t".join(line)+"\n")
if DEBUG:
    print (chrom+":"+str(p1hap)+","+str(p2hap))
                
                


