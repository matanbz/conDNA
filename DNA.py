
#!/usr/bin/python
# -*- coding: utf-8 -*-
#from glob import glob
import os, sys, csv, re
from numpy import array, loadtxt, savetxt,  amin, amax, uint8,arange,logical_or, mean
from numpy.random import random,shuffle
from time import time,sleep,localtime
#import pymorph as p
from PIL import Image
import matplotlib.cm as cm

#length of first generation sticky ends
fgsn=5
tvsn=8
typ='xj'
#typ='tw'
model="new"
testdir="/windows/D/NY/python/DNAmodule/"
if sys.platform=="win32":
    testdir="D:/NY/python/DNAmodule/"
templatedir=testdir+"templates/"
allstrandfile=testdir+"strands-new-cond-cnv-low8.csv"
if typ=='xj':
  tvsn=5
  fgsn=3


#bodylen=6
#selen=4

extract = lambda x, y: dict(zip(x, map(y.get, x)))
cols={10:"#008000;",11:"#aa0044;",12:"#d45500;",13:"#668000;",14:"#ff00ff;",15:"#00ff00;",16:"#008080;",17:"#000080;",18:"#800000;",20:"#800080;"}


def noRepetitions(lst,length, cgs=0):
    """takes a list of possible sequences and removes those containing sections of repeating bases of length 'length' and, optionally, those containing local CG accumulations of length 'cgs' (default: disabled with cgs=0).
    Input: string list for DNA, integer for length, returns reduced string list"""
    newitems=[]
    for item in lst:
        flag=0
        for base in ["G","A","T","C"]:
            if item.find(base*length)>-1:
                flag=1
                break
        if cgs>0:
            for pos in range(len(item)-cgs):
                if item[pos:pos+cgs].count("C")+item[pos:pos+cgs].count("G")==cgs:
                    flag=1
                    break

        if flag==0: newitems=newitems+[item]
    return newitems

def CGContent(lst,lower,upper,  verbose=True):
    """reduces a list of DNA sequences to ones with CG content in the range (lower,upper). Verbose prints the actual content"""
    newlist=[]
    for item in lst:
        cont=0
        for char in item:
            if char=="G" or char=="C":cont=cont+1
        cont=cont/float(len(item))
        if (lower<cont) & (cont<upper):
            newlist=newlist+[item]
            if verbose==True: print cont
    return newlist

#def findPalindromic(lseq):
 #   """takes a list of DNA sequences and returns the ones that are entirely palindromic. Why the hell do I need this function?"""
  #  ispal=[]
   # for seq in lseq:
	#seq = seq.upper()        
	#if seq[:len(seq)/2]==get_complement(seq[-(len(seq)/2):][::-1]): ispal=ispal+[seq]
#    return ispal


def genSingle(char,totallength):
    """generates a DNA duplex with maximised randomness (input random char length 'char' and intended length of strand)"""
    f=randBase(char)
    fcomp=get_complement(f)
    j=0
    while len(f)<totallength and j<20000:
        j=j+1
        breakflag=1
        bases=["G","A","T","C"]
        shuffle(bases)
        i=0
        for base in bases:
            i=i+1
            f=f+base
            fcomp=get_complement(f)
            if f[:-1].find(f[-char:])+f[:-1].find(f[-char:][::-1])+fcomp[:-1].find(f[-char:])+fcomp[:-1].find(f[-char:][::-1]) == -4:
                breakflag=0
                print i
                break
            f=f[:-1]
        if breakflag==1:
            f=f[:-1]
    print "iterations ", j
    return f

def findSortList(item,mylist,num):
    """looks for the occurrence of an item in a given list, starting from num, and removes it.
    Helper function for PossCompl"""
    for i in range(num, len(mylist)):
        if mylist[i]==item:
            mylist.pop(i)
            return i-1
            break

        #
def find_key(dic, val):
    #
    """return the key of dictionary dic given the value"""
    #
    return [k for k, v in dic.iteritems() if v == val][0]

def pylabPilPal(cmap_name,N):
    cmap = cm.get_cmap(cmap_name, N)
    return (cmap(arange(N))[:,0:3].reshape(1,-1)[0]*256).astype(uint8)


def uniq(seq, idfun=None):
    """return all unique items in a list (preserves list order)"""
    # order preserving
    if idfun is None:
        def idfun(x): return x
    seen = {}
    result = []
    for item in seq:
        marker = idfun(item)
        if marker in seen: continue
        seen[marker] = 1
        result.append(item)
    return result


def generate_possibles(length):
    b=['A','C','G','T']
    for i in range(2, length+1):
        a=['A'+s for s in b]
        c=['C'+s for s in b]
        g=['G'+s for s in b]
        t=['T'+s for s in b]
        b=a
        b.extend(c)
        b.extend(g)
        b.extend(t)
    return b

def sortedDictValues(adict):
    """sort a dictionary's values by the order of the keys in it."""
    keys = adict.keys()
    keys.sort()
    return map(adict.get, keys)


def OD2C(OD, lDNA,DNAvol=5., H2Ovol=500.):
  return (OD*3.5e-5)/(DNAvol*(1000/H2Ovol)*0.001*lDNA*330)

def randBase(length=1, bases=["G","A","T","C"]):
    result=""
    for i in range(length):
        a=int(random(1)*len(bases))
        result=result+bases[a]
    return result

def allstrands(allstrandfile,redux=False):
    """reads sequence identifiers and sequences from CSV into a dictionary."""
    strandsdict={}
    csvfile = open(allstrandfile)
    dialect = csv.Sniffer().sniff(csvfile.read(1024))
    csvfile.seek(0)
    data = csv.reader(csvfile, dialect)
    # ... process CSV file contents here ...
    for row in data:
        strandsdict[row[0]]=row[1]
    if redux:
        for key in strandsdict.keys():
            strandsdict[key]=strandsdict[key].replace(" ","")
            strandsdict[key]=strandsdict[key].replace("5'-","")
            strandsdict[key]=strandsdict[key].replace("-3'","")

    return strandsdict

def get_complement(string):
    """returns a complementary DNA sequence, e.g. GCTA for TAGC. Mind the sequence inversion!"""
    string=string.replace('G','1')
    string=string.replace('A','2')
    string=string.replace('T','3')
    string=string.replace('C','4')
    string=string.replace('1','C')
    string=string.replace('2','T')
    string=string.replace('3','A')
    string=string.replace('4','G')
    string=string.replace('X','Z')#for 'blank' sticky ends and cinnamate
    return string[::-1]

    
    
def test_strands(strands,charl=5,verbose=False):
  forb=0
  for n in range(len(strands)):
    for i in range(len(strands[n])-charl+1):
      test=get_complement(strands[n][i:i+charl])
      for j in range(n)+range(n+1,len(strands)):
	for k in range(len(strands[j])-charl+1):
	  if strands[j][k:k+charl]==test: 
	    if verbose: print 'forbidden: ', key, i, j
	    forb=forb+1
  return forb
    
    
    
def hamming(s1, s2):
    nc=0
    for i in range(len(s1)):
        if s2[-i]==get_complement(s1[i]): nc=nc+1
        if s2[i]==s1[i]:nc=nc+1
    return nc/float(len(s1))

def countercheck(sequence, strfile=allstrandfile,revtest=True, bmmatch=5,  emmatch=7):
    stdic=allstrands(strfile, redux=True)
    psdic=allstrands(testdir+"DNAspecs/pseudohelices.csv")
    out=""
    for key in stdic.keys():
        strand=stdic[key]
        rstrand=stdic[key][::-1]
        for i in range(len(strand)-len(sequence)+1):
            num=0
            rnum=0
            for j in range(len(sequence)):
                if not sequence[j]==get_complement(strand)[j+i]: num=num+1
                if not sequence[j]==get_complement(rstrand)[j+i]: rnum=rnum+1

            if revtest:
                if not i in [0, len(strand)-len(sequence)] and rnum<=bmmatch: out=out+ "Warning! insufficient mismatch, %s[%d]:%d\n"%(key,i,num)
                if i in [0, len(strand)-len(sequence)] and rnum<=emmatch: out=out+ "Warning! insufficient SE mismatch, %s[%d]:%d\n"%(key,i,num)
            if not i in [0, len(strand)-len(sequence)] and num<=bmmatch: out=out+ "Warning! insufficient reverse mismatch, %s[%d]:%d\n"%(key,i,rnum)
            if i in [0, len(strand)-len(sequence)] and num<=emmatch: out=out+ "Warning! insufficient SE reverse mismatch, %s[%d]:%d\n"%(key,i,rnum)
    for key in psdic.keys():
        strand=psdic[key]
        rstrand=psdic[key][::-1]
        for i in range(len(strand)-len(sequence)+1):
            num=0
            rnum=0
            for j in range(len(sequence)):
                if not sequence[j]==get_complement(strand)[j+i]: num=num+1
                if not sequence[j]==get_complement(rstrand)[j+i]: rnum=rnum+1
            if revtest:
                if num<=bmmatch: out=out+ "Warning! insufficient mismatch, %s[%d]:%d\n"%(key,i,num)
            if rnum<=bmmatch: out=out+ "Warning! insufficient reverse mismatch, %s[%d]:%d\n"%(key,i,rnum)
    return out

def TMelt(string, mg="12.5e-3", conc="5e-7", meltdir=testdir+'../MELTING5.0.3/executable', comp='', verbose=False):
    """passes a DNA sequence to the external MELTING Java app and returns
    the approximate melting temperature, entropy and enthalpy. Default values for Mg++
    and DNA concentration can be adjusted as well as the location of the melting5.jar
    executable. Note: This function assumes your system's decimal separator to be a dot '.', not a comma ','!"""
    thisdir=os.getcwd()
    #os.chdir(meltdir)
    #if os.path.exists('melttemp.txt'): os.remove('melttemp.txt')
    if len(comp) >0: os.system('melting -S%s -C%s -G%s -P%s -Hdnadna -Omelttemp.txt'%(string,comp[::-1],mg,conc))
    #os.system('java -jar melting5.jar -S %s -C %s -E Mg=%s -P %s -H dnadna -O melttemp.txt'%(string,comp[::-1],mg,conc))
    else:os.system('melting -S%s -G%s -P%s -Hdnadna -Omelttemp.txt'%(string,mg,conc))
    #os.system('java -jar melting5.jar -S %s -E Mg=%s -P %s -H dnadna -O melttemp.txt'%(string,mg,conc))
    f=open('melttemp.txt')
    dat=f.read()
    f.close()
    if verbose:
        print os.getcwd()
        print 'melting -S%s -G%s -P%s -Hdnadna -Omelttemp.txt'%(string,mg,conc)
        print dat
    dat=dat.replace(',','')
    if dat.find('SEVERE') > -1: res={'temp':-300,'S':-1,'H':-1}
    else:
        res={}
        se=re.finditer(': .*? J.mol', dat)
        res['H']=float(se.next().group(0)[2:-6])
        res['S']=float(se.next().group(0)[2:-6])
        res["temp"]=float(re.search('temperature: .*? deg C', dat).group(0)[13:-6])
    os.chdir(thisdir)
    if os.path.exists('melttemp.txt'): os.remove('melttemp.txt')
    return res



def TMeltSelect(lst, interv=(20, 40), comp=[], mg="12.5e-3", conc="5e-7", meltdir=testdir+'../MELTING5.0.3/executable',  verbose=False, threads=1):
    mdir={}
    thisdir=os.getcwd()
    #os.chdir(meltdir)
    for i in range(len(lst)/threads):
        if verbose: 
	  if len(comp)>0: print ';'.join(["melting -S'%s' -C'%s' -G%s -P%s -Hdnadna -Omelttemp%d.txt"%(lst[threads*i+j],comp[threads*i+j],mg,conc, j) for j in range(threads)])
	  else: print ';'.join(["melting -S'%s' -G%s -P%s -Hdnadna -Omelttemp%d.txt"%(lst[threads*i+j],mg,conc, j) for j in range(threads)])
        if len(comp) >0: os.system(';'.join(["melting -S'%s' -C'%s' -G%s -P%s -Hdnadna -Omelttemp%d.txt"%(lst[threads*i+j],comp[threads*i+j],mg,conc, j) for j in range(threads)]))
        else:os.system(';'.join(["melting -S'%s' -G%s -P%s -Hdnadna -Omelttemp%d.txt"%(lst[threads*i+j],mg,conc, j) for j in range(threads)]))
        for j in range(threads):
	  f=open('melttemp%d.txt'%j)
	  dat=f.read()
	  f.close()
	  dat=dat.replace(',','')
	  if dat.find('SEVERE') > -1: res={'temp':-300,'S':-1,'H':-1}
	  else:
	    res={}
	    se=re.finditer(': .*? J.mol', dat)
	    res['H']=float(se.next().group(0)[2:-6])
	    res['S']=float(se.next().group(0)[2:-6])
	    res["temp"]=float(re.search('temperature: .*? deg C', dat).group(0)[13:-6])
            if interv[0]<res['temp']<interv[1]: 
                mdir[lst[i*threads+j]]=res['temp']
                if verbose: print "%d: s.e. %s, T melt. %.02f"%(i*threads+j,  lst[i*threads+j],  res['temp'])
    os.system('rm melttemp?.txt')
    os.chdir(thisdir)
    return mdir

def ExchangeSticky(trange=(10, 45)):
    sevensticky=["tw-ICN4-11","tw-ICN4-13", "tw-AA2-11", "tw-AA2-13","tw-BB2-11", "tw-BB2-13","tw-ICN4-03","tw-ICN4-05", "tw-AA2-03","tw-AA2-05","tw-u3B-03","tw-u3B-05"]
    replacedir={}
    for end in sevensticky: replacedir[end]=""

class seedtile:
    """Definitions and helper functions for Tong's seed tiles."""
    def __init__(self,newID):
        self.seedID=newID
        self.type=self.seedID[0].upper()
        try: self.num=int(self.seedID[1])
        except ValueError: self.num=0
        self.stringIDs={}
        self.stringSequences={}
        self.subsequences={}
        self.cand=[]
        self.matches={"10B":"15C","10C":"11B","10D":"12E","10F":"13B","10G":"14B","15B":"16H","11C":"16G", "11D":"16F","12C":"16E","13C":"16D","13D":"16C","18C":"16B","17B":"20A","17C":"11E","17D":"12B", "17E":"12A","17F":"13E","17G":"18B",'11E': '17C', '12A': '17E', '12B': '17D', '16H': '15B', '12E': '10D', '13E': '17F', '11B': '10C', '16D': '13C', '16E': '12C', '16F': '11D', '16G': '11C', '16B': '18C', '16C': '13D', '20A': '17B', '15C': '10B', '14B': '10G', '18B': '17G', '12D': '10E'}
        self.semfile=testdir+'sematches.txt'
        if model=='new':self.semfile=testdir+'sematches-new.txt'
        self.connect=[("10D","10E"),("12D","12E"),("11C","11D"),("16G","16F"),("13C","13D"),("16C","16E"),("17D","17E"),("12A","12B")]
        self.sticky=["11F","13F","11A","13A"]
        self.templatefile="SEEDmidtemplate%d.svg"%tvsn
        if self.num==1: self.templatefile="A1template%d.svg"%tvsn
        if self.num==7: self.templatefile="B7template%d.svg"%tvsn
        if typ=='xj': 
	  self.templatefile="SEEDmidtemplatexj.svg"
          if self.num==1: self.templatefile="A1templatexj.svg"
          if self.num==7: self.templatefile="B7templatexj.svg"
        self.templateStrands={}
        self.stringNums=[10,11,12,13,14,15,16,17,18,20]

    def init(self, newID):
        """initialises the tile to a new template"""
        self.seedID=newID
        self.type=self.seedID[0].upper()        
        if self.num==7: self.templatefile="B7template%d.svg"%tvsn
        if typ=='xj': 
	  self.templatefile="SEEDmidtemplatexj.svg"
          if self.num==1: self.templatefile="A1templatexj.svg"
          if self.num==7: self.templatefile="B7templatexj.svg"
        self.num=int(self.seedID[1])
        self.stringIDs={}
        self.stringSequences={}

    def readSVGtemplate(self, verbose=False):
        o=open(templatedir+self.templatefile)
        svg=o.read()
        o.close()
        lengths={}
        for num in self.stringNums:
            index=0
            for alph in "ABCDEFGHIJKLMNOPQRSTUVWXYZ": 
                n=re.search("%02d%s.*?<"%(num,alph), svg)
                if type(n).__name__!='SRE_Match':
                    pass
                else:
                    n=n.group(0)[:-1]
                    aug=len(n)
                    if n[-1]=="R":  self.templateStrands["%02d%s"%(num,alph)]=(index+aug,index)
                    if n[-1]=="N":  self.templateStrands["%02d%s"%(num,alph)]=(index,index+aug)
                    if n[-1]=="T":
                        aug=int(n[-2])
                        self.templateStrands["%02d%s"%(num,alph)]=(index,index+aug)
                        
                    if aug<4: #infer from 'B' strand
                        if re.search("%02dB.*?<"%num, svg).group(0)[-1]=='R':self.templateStrands["%02d%s"%(num,alph)]=(index+aug,index)
                        else:self.templateStrands["%02d%s"%(num,alph)]=(index,index+aug)

                    if verbose:
                        print "%s, %02d"%(n,aug)
                    index=index+aug
            if verbose: print "Strand %d, length %d"%(num, index)
            lengths[num]=index
        return lengths    

    def replaceStrand(self, s1, s2):
        num=find_key(self.stringIDs, s1)
        self.stringIDs[num]=s2
        self.getStrands()


    def getIDs(self):
        """identifies SEED tile strings according to the seed ID and
        writes them into the stringIDs dictionary"""
        self.stringIDs[12]="tw-ICEN-12"
        self.stringIDs[20]="tw-ICEN-20"
        a=[10,14,15,16,17,18]
        for i in a:
            self.stringIDs[i]="tw-%s-%02d"%(self.seedID,i)
        a=[11,13]
        for i in a:
            self.stringIDs[i]="tw-%s%d-%02d"%(2*self.seedID[0],tvsn, i)
        if self.seedID=="A1":
            self.stringIDs[10]="tw-ICEN-10"
            self.stringIDs[11]="tw-ICN%d-11"%tvsn
            self.stringIDs[13]="tw-ICN%d-13"%tvsn
            self.stringIDs[17]="tw-ICEN-17"
        if self.seedID=="B7":
            self.stringIDs[16]="tw-ICEN-16"
        if typ=='xj':
	  for i in [10,11,12,13,14,17,20]:self.stringIDs[i]=self.stringIDs[i]+'xj'
        #print self.stringIDs

    def getStrands(self, verbose=False):
        """reads DNA sequences from CSV file and pairs them with identifiers (dictionary)"""
        self.getIDs()
        strandsd=allstrands(allstrandfile)
        self.stringSequences=extract(self.stringIDs.values(),strandsd)
        for key in self.stringSequences.keys():
            if verbose: print key
            self.stringSequences[key]=self.stringSequences[key].replace(" ","")
            self.stringSequences[key]=self.stringSequences[key].replace("5'-","")
            self.stringSequences[key]=self.stringSequences[key].replace("-3'","")




    def getSubseq(self,num, comp=False):
        """splits all sequences into possible subsequences of length num (for num=5, as in 12345, 23456, 34567 ...).
        Identifiers AABBB, with AA sequence number, BBB position of first base in string.
        Negative identifiers denote reversed sequences.
        If comp is set, the function returns complementary subsequences instead."""
        ##add reverses!!!
        self.getStrands()
        for  key in self.stringSequences.keys():
            string=self.stringSequences[key].replace("5'-","")
            string=string.replace("-3'","")
            string=string.replace(" ","")
            for i in range(len(string)+1-num):
                if comp:
                    self.subsequences[get_complement(string[i:i+num])]=int(key[-2:])*1000+i
                    self.subsequences[get_complement(string[i:i+num][::-1])]=-int(key[-2:])*1000+i
                else:
                    self.subsequences[string[i:i+num]]=int(key[-2:])*1000+i
                    self.subsequences[string[i:i+num][::-1]]=-int(key[-2:])*1000+i


    def PossCompl(self, num):
        self.getSubseq(num,comp=True)
        vals=self.subsequences.keys()
        vals.sort()
        vals=uniq(vals)
        self.cand=generate_possibles(num)
        print "before: %d"%len(self.cand)
        a=0
        for string in vals:
            a=findSortList(string,self.cand,a)
        print "after: %d"%len(self.cand)


    def fetchAllSticky(self,charlen):
        allsticky={}
        strandsd=allstrands(allstrandfile)
        for key in strandsd.keys():
            strand=strandsd[key]
            strand=strand.replace(" ","")
            strand=strand.replace("5'-","")
            strand=strand.replace("-3'","")
            strandsd[key]=strand

        for tile in ["A1", "B2","B3", "A4", "B5", "A6", "B7"]:
            for i in [10,14,15,16,17,18]:
                thisid="tw-%s-%02d"%(tile,i)
                try:
                    allsticky[thisid]=strandsd[thisid][:16+charlen-1]
                    if thisid=="tw-B7-18":
                        allsticky[thisid+"e"]=strandsd[thisid][-(6+charlen-1):][::-1]
                    if allsticky[thisid][:6]=="TTTTTT":
                        allsticky[thisid]=allsticky[thisid][:6+charlen-1]
                except KeyError:
                    pass
            if model=='new':#why doesn't this register?
                for i in [10,14,15,16,17,18]:
                    thisid="tw-%s-%02dn"%(tile,i)#use a newer version if possible (what to do about dated strands?)
                    try:
                        allsticky[thisid]=strandsd[thisid][:16+charlen-1]
                        dummy=allsticky.pop(thisid[:-1])
                        if thisid=="tw-B7-18":
                            allsticky[thisid+"e"]=strandsd[thisid][-(6+charlen-1):][::-1]
                        if allsticky[thisid][:6]=="TTTTTT":
                            allsticky[thisid]=allsticky[thisid][:6+charlen-1]
                    except KeyError:
                        pass

        for thisid in ["tw-ICN%d-11"%tvsn,"tw-ICN%d-13"%tvsn, "tw-AA%d-11"%tvsn, "tw-AA%d-13"%tvsn,"tw-BB%d-11"%tvsn, "tw-BB%d-13"%tvsn,"tw-ICN%d-03"%tvsn,"tw-ICN%d-05"%tvsn, "tw-AA%d-03"%tvsn,"tw-AA%d-05"%tvsn,"tw-u3B%d-03"%tvsn,"tw-u3B%d-05"%tvsn,"tw-BB%d-03"%tvsn,"tw-BB%d-05"%tvsn]:
            allsticky[thisid+"b"]=strandsd[thisid][:tvsn+charlen-1]
            allsticky[thisid+"e"]=strandsd[thisid][-(tvsn+charlen-1):][::-1]


        if fgsn==0:
            for thisid in ["tw-u2AB-02","tw-u2AB-06","tw-u2AB-09", "tw-u2AB-10", "tw-u2AB-16", "tw-u2AB-17"]:#find solution for new model
                allsticky[thisid+"b"]=strandsd[thisid][:9+charlen-1]
                allsticky[thisid+"e"]=strandsd[thisid][-(9+charlen-1):][::-1]
        else:
            for i in [1, 2, 6, 7, 8, 9, 10,14,15,16,17,18, 19]:
                thisid="tw-u2AB%d-%02d"%(fgsn,i)#use a newer version if possible (what to do about dated strands?)
                allsticky[thisid]=strandsd[thisid][:fgsn+charlen-1]
        return allsticky

    def fetchBody(self):
        body={}
        bfile=open(testdir+"blunt-seed.txt")
        for line in bfile:
            body[line.split("\t")[0]]=line.split("\t")[1].replace('\n','')

        return body

    def fetchMatches(self):
        body={}
        bfile=open(self.semfile)
        for line in bfile:
            line=line.replace('#t','%d'%tvsn)
            body[line.split("\t")[0]]=line.split("\t")[1][:-1]
            body[line.split("\t")[1][:-1]]=line.split("\t")[0]
        return body

    def ReplaceSticky(self,seqID,selen,bm=6,sm=5,palin=0, revtest=False,  replacements={},tlim=300,  nlim=20000, verbose=False, constraints=''):
        """argument: strand ID like 'tw-B7-10'. Sticky ends are _always_ at 5' end!
         options: bm=6, maximum number of random body matches in a row
              sm=5, maximum number of random sticky end matches in a row
              palin=0, number of palindromic bases (for hairpins) at each end. This results           in a shortened path matrix
              revtest=False, check for reverse compelentarity (not necessary, parallel DNA binding is unfavourable
              replacements: dictionary of already generated replacement SE's {'tw-A6-17':'AGTAGGC...'}
              tlim=300, time limit in seconds
              nlim=20000, max. number of generated possibilities
              verbose=False, prints all matches (A LOT!!!!)
              constraints='', string of pre-determined bases like 'UYAXXXGC', 'S' and 'W' for G/C or A/T, 'X' for free choice, Replaces rows in the path matrix. """
        t1=time()
        self.allsticky=self.fetchAllSticky(max(bm,sm))
        if len(replacements)>0:
            for key in replacements.keys():
                if verbose==True: print key,":",self.allsticky[key]
                self.allsticky[key]=replacements[key]+self.allsticky[key][-max(bm,sm)+1:]
                if verbose==True: print "replaced by", self.allsticky[key]
                self.getStrands()
        body=self.fetchBody()
        if self.seedID in ["A1","B2", "B3", "A4", "B5", "A6", "B7", "A''", "B''" , "I''"]:
            #comptile=firstgentile(self.type+"'")
            comptile=newfgtile(self.type+"'")
        if self.seedID in ["A'", "B'", "I'"]:
            #comptile=secgentile(self.type+"''")
            comptile=newsgtile(self.type+"''")
        compbody=comptile.fetchBody()
        matches=self.fetchMatches()
        # identify sticky end
        stubseq=self.allsticky[seqID][-max(bm,sm)+1:]
        print len(stubseq)
        #take it from allsticky... and delete
        del self.allsticky[seqID]

        try: del self.allsticky[matches[seqID]]
        except KeyError: pass

        self.basematrix=[]
        for i in range(0,selen-palin):
            bases=["G","A","T","C"]
            shuffle(bases)
            self.basematrix=self.basematrix+bases
        self.basematrix=array(self.basematrix).reshape(selen-palin,4)
        if len(constraints)>0:
            for i in range(0,selen-palin):
                try:
                    if constraints[i]=='S': self.basematrix[i,:]=['C', 'C', 'G', 'G']
                    if constraints[i]=='W': self.basematrix[i,:]=['A', 'A', 'T', 'T']
                    if constraints[i] in ['G','A','T','C']: self.basematrix[i,:]=[constraints[i]]*4
                except IndexError: break

        possibles=[]
        baseind=0
        seind=0
        count=0
        path=[0]
        alloverflag=0

        while(alloverflag==0):
            matchflag=0
            testseq=stubseq
            for i in range(len(path)):
                testseq=self.basematrix[i,path[i]]+testseq
            seseq=get_complement(testseq[:sm])
            for m,n in self.allsticky.iteritems():
                count=count+1
                if revtest:
                    rev=n.find(seseq[::-1])
                else: rev=-1
                if rev+n.find(seseq)>-2:
                    if verbose: print "s.e.",count,m,rev,n.find(seseq),path, len(possibles)
                    matchflag=1
                    break
            bodyseq=get_complement(testseq[:bm])
            for m,n in body.iteritems():
                count=count+1
                if count%10000==0: print len(possibles)
                if revtest:
                    rev=n.find(bodyseq[::-1])
                else: rev=-1
                if rev+n.find(bodyseq)>-2:
                    if verbose: print "body",count,m,rev,n.find(bodyseq),path, len(possibles)
                    matchflag=1
                    break

        #-----------conditions for aborting-----------------------
            if mean(path)==3 and len(path)==selen-palin-1:
                print len(path), baseind,count
                print "no combinations left, aborting!"
                alloverflag=1
                break
            else:
                if len(possibles)>nlim:
                    print "%d possible sequences, aborting!"%nlim
                    alloverflag=1
                    break
                if time()-t1>tlim:
                    print "time limit of %d seconds exceeded, aborting!"%tlim
                    alloverflag=1
                    break
                if matchflag==0 and seind==selen-palin-1: #base passed testss
                    nmatchflag=0
                    if palin>0:
                        # make sure that middle bit _not_ palindromic
                        for i in range((selen-2*palin)/2):
                            if get_complement(testseq[:selen-2*palin][-i-1])==testseq[:selen-2*palin][i]:
                                nmatchflag=1 #
                                break #
                        testseq=get_complement(testseq[selen-2*palin:selen-palin])+testseq #append the palindromic part's complement
                        for c in range(palin):
                            seseq=get_complement(testseq[c:c+sm])
                            for m,n in self.allsticky.iteritems():
                                count=count+1
                                if revtest:
                                    rev=n.find(seseq[::-1])
                                else: rev=-1
                                if rev+n.find(seseq)>-2:
                                    if verbose: print "Palindrome SE match",count,m,rev,n.find(seseq)
                                    nmatchflag=1
                                    break
                            if nmatchflag==1: break
                            seseq=get_complement(testseq[c:c+bm])
                            for m,n in body.iteritems():
                                count=count+1
                                if revtest:
                                    rev=n.find(seseq[::-1])
                                else: rev=-1
                                if rev+n.find(seseq)>-2:
                                    if verbose: print "Palindrome body match",count,m,rev,n.find(seseq) #TODO: probably only [::-1] case necessary
                                    nmatchflag=1 #f
                                    break
                            if nmatchflag==1: break
                        #TODO: where is the countercheck routine?
                    try:
                        match=matches[seqID]
                        match=get_complement(testseq[:selen])+match[selen:]
                    except KeyError:
                        match=get_complement(testseq[:selen])
                    for a in range(len(match)-sm):
                        seseq=get_complement(match[a:a+sm])
                        for m,n in self.allsticky.iteritems():
                            count=count+1
                            if n.find(seseq)>-1:
                                if verbose: print "comp s.e.",count,m,n.find(seseq)#,n.find(seseq[::-1])
                                nmatchflag=1
                                break
                        if nmatchflag==1: break
                    for a in range(len(match)-bm):
                        seseq=match[a:a+bm]
                        for m,n in body.iteritems():
                            count=count+1
                            if n.find(seseq)>-1:
                                if verbose: print "comp body",count,m,n.find(seseq)#,n.find(seseq[::-1])
                                nmatchflag=1
                                break
                        if nmatchflag==1: break
                    if nmatchflag==0: possibles=possibles+[testseq[:selen]]
                    else:
                        if verbose==True: print "complement match!"

                if matchflag==1 or seind==selen-palin-1:
                    while(baseind==3):
                        testseq=testseq[1:]#when?
                        seind=seind-1
                        path=path[:-1]
                        try: baseind=path[-1]
                        except IndexError:
                            print "End of line!"
                            alloverflag=1
                            break
                    if alloverflag==1: break
                    baseind=baseind+1
                    testseq=bases[baseind]+testseq[1:]
                    path[-1]=baseind

                if matchflag==0 and seind<selen-palin-1:
                    seind=seind+1
                    baseind=0
                    testseq=bases[baseind]+testseq
                    path=path+[baseind]


        print "Finished in %.2f minutes, resulting in %d possible sequences from %d string comparisons"%((time()-t1)/60.0,len(possibles),count)
        return possibles


    def testStructure(self):
        self.getStrands()
        self.readSVGtemplate()
        failcount=0
        stickcount=0
        for key in self.matches.keys():
            si=self.templateStrands[key]
            try:
                sii=self.templateStrands[self.matches[key]]
                if si[1]<si[0]:
                    seqi=self.stringSequences[self.stringIDs[int(key[0:2])]][si[1]:si[0]][::-1]
                else:
                    seqi=self.stringSequences[self.stringIDs[int(key[0:2])]][si[0]:si[1]]
                if sii[1]<sii[0]:
                    seqii=self.stringSequences[self.stringIDs[int(self.matches[key][0:2])]][sii[1]:sii[0]][::-1]
                else:
                    seqii=self.stringSequences[self.stringIDs[int(self.matches[key][0:2])]][sii[0]:sii[1]]
                if seqi != get_complement(seqii):
                    print "Match? %s:%s %s:%s"%(key,seqi,self.matches[key],get_complement(seqii))
                    failcount+=1
            except KeyError: print "Key Error: %s,%s"%(key,self.matches[key])
        print "Missed matches: %s"%failcount
        m=firstgentile("%s'"%self.type)
        try: dummy= self.num
        except AttributeError:
            print "First gen?"
            return
        if self.num==1: m=firstgentile("I'")
        m.getStrands()
        m.readSVGtemplate()
        for i in range(4):
            #print i
            si=self.templateStrands[self.sticky[i]]
            sii=m.templateStrands[m.sticky[(i+2)%4]]
            if si[1]<si[0]:
                seqi=self.stringSequences[self.stringIDs[int(self.sticky[i][0:2])]][si[1]:si[0]][::-1]
            else:
                seqi=self.stringSequences[self.stringIDs[int(self.sticky[i][0:2])]][si[0]:si[1]]
            if sii[1]<sii[0]:
                seqii=m.stringSequences[m.stringIDs[int(m.sticky[(i+2)%4][0:2])]][sii[1]:sii[0]][::-1]
            else:
                seqii=m.stringSequences[m.stringIDs[int(m.sticky[(i+2)%4][0:2])]][sii[0]:sii[1]]
            if seqi != get_complement(seqii):
                print "Sticky ends? %s:%s %s:%s"%(self.sticky[i],seqi,m.sticky[(i+2)%4],get_complement(seqii))
                stickcount+=1
        return (failcount, stickcount)

    def wikiOutput(self):
        """format a tile image and its sequences for Tiddly Wiki copy and paste output."""
        self.getStrands()
        output='[img(100%%+,+)[img/col_%s.png]]\n'%(self.seedID)
        for i in self.stringIDs.keys():
            iID=''
            key=self.stringIDs[i]
            out="%02d. %s. %s"%(i, key,self.stringSequences[key])
            if len(out)>80:
                output=output+"@@color:"+cols[i]+"font-size:0.92em;"+out[:int(len(out)/2)+8]+"@@\n@@color:"+cols[i]+"font-size:0.92em;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;"+out[int(len(out)/2)+8:]+"@@\n"
            else:
                output=output+"@@color:"+cols[i]+"font-size:0.92em;"+out+"@@\n"
        print output

    def fillTemplate(self, fname, verbose=False):
	if self.stringSequences=={} or self.stringIDs=={}:  self.getStrands()
        ls=self.readSVGtemplate()
        o=open(templatedir+self.templatefile)
        svg=o.read()
        o.close()
        for num in self.stringNums:
            index=0
            count=0
            if len(self.stringSequences[self.stringIDs[num]])!=ls[num]: 
	      print "warning: strand %s doesn't match template!"%self.stringIDs[num]
            for alph in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
                try:
                    si=self.templateStrands["%02d%s"%(num,alph)]
                    if si[1]<si[0]:
                        seq=self.stringSequences[self.stringIDs[num]][si[1]:si[0]][::-1]
                        if si[0]-si[1]>3: plac="%02d%s"%(num,alph+"X"*(si[0]-si[1]-4)+"R")
                        else: plac="%02d%s"%(num,alph)
                    else:
                        seq=self.stringSequences[self.stringIDs[num]][si[0]:si[1]]
                        if si[1]-si[0]>3: plac="%02d%s"%(num,alph+"X"*(si[1]-si[0]-4)+"N")
                        else: plac="%02d%s"%(num,alph)
                    if verbose: print "%s %s"%(plac,seq)
                    svg=svg.replace(plac,seq)
                except KeyError: pass
        svg=svg.replace("TITLE",self.seedID)
        o=open(fname, "w")
        o.write(svg)
        o.close()
        os.system('inkscape -z -A "'+fname[:-4]+'.pdf" "'+fname+'"')
        print fname




class firstgentile(seedtile):
    def __init__(self,newID):
        seedtile.__init__(self,newID)
        del self.num
        if self.type=="B":
            self.matches={"01A":"02G","03I":"02F", "04A":"O2E","04F":"02D", "05I":"02C","19A":"02B","06B":"01B","06C":"03H","06D":"03C","06E":"04G","06F":"05H","06G":"05C","06H":"07A","08A":"09G","03B":"09F","04H":"09E","04M":"09D","05B":"09C","07B":"09B","03F":"03D","04B":"04D","04I":"04K","05F":"05D"}
            self.connect=[]
            self.sticky=["03A","05A","03J","05J"]
            self.templatefile="B'template.svg"
        else:
            self.matches={"01A":"02G","03E":"02F", "04A":"O2E","04B":"02D", "05E":"02C","19A":"02B","06B":"01B","06C":"03D","06D":"03C","06E":"04C","06F":"05D","06G":"05C","06H":"07A","08A":"09G","03B":"09F","04H":"09E","04M":"09D","05B":"09C","07B":"09B"}
            self.connect=[("04A","04B"),("02D","02E"),("06C","06D"),("03C","03D"),("06F","06G"),("05C","05D"),("04D","04E"),("09D","09E")]
            self.sticky=["03A","05A","03F","05F"]
            if self.type=="I": self.templatefile="I'template.svg"
            else: self.templatefile="A'template.svg"
        self.stringNums=[1,2,3,4,5,6,7,8,9,19]
        
	  

        

    def getIDs(self):
        """identifies SEED tile strings according to the seed ID and
        writes them into the stringIDs dictionary"""
        for i in [1, 2, 6, 7, 8, 9, 19]:
            self.stringIDs[i]="tw-u2AB-%02d"%i
        self.stringIDs[3]="tw-AA%d-03"%tvsn
        self.stringIDs[4]="tw-ICEN-04"
        self.stringIDs[5]="tw-AA%d-05"%tvsn
        
        if self.type=="I":
            self.stringIDs[1]="tw-ICEN-01"
            self.stringIDs[2]="tw-u1I-02"
            self.stringIDs[3]="tw-ICN%d-03"%tvsn
            self.stringIDs[4]="tw-ICEN-04"
            self.stringIDs[5]="tw-ICN%d-05"%tvsn
            self.stringIDs[6]="tw-u1I-06"
            self.stringIDs[8]="tw-ICEN-08"
            self.stringIDs[9]="tw-u1I-09"
            self.stringIDs[19]="tw-ICEN-19"

        if self.type=="B":
            self.stringIDs[3]="tw-u3B%d-03"%tvsn
            self.stringIDs[4]="tw-u3AB-04"
            self.stringIDs[5]="tw-u3B%d-05"%tvsn

    def fetchBody(self):
        body={}
        if self.type=="B":
            bfile=open(testdir+"blunt-fgen-hp.txt")
        else:
            bfile=open(testdir+"blunt-fgen.txt")
        for line in bfile:
            body[line.split("\t")[0]]=line.split("\t")[1][:-1]
        return body


class newfgtile(firstgentile):
    def __init__(self,newID):
        seedtile.__init__(self,newID)
        del self.num
        if self.type=="B":
            self.matches={"01B":"02G","03I":"02F", "04A":"O2E","04F":"02D", "05I":"02C","19A":"02B","06B":"01C","06C":"03H","06D":"03C","06E":"04G","06F":"05H","06G":"05C","06H":"07B","08B":"09G","03B":"09F","04H":"09E","04M":"09D","05B":"09C","07C":"09B","03F":"03D","04B":"04D","04I":"04K","05F":"05D"}
            self.connect=[]
            self.sticky=["03A","05A","03J","05J"]
            self.templatefile="B'newtemplate%d-%d.svg"%(fgsn,tvsn)
        else:
            self.matches={"01B":"02G","03E":"02F", "04A":"O2E","04B":"02D", "05E":"02C","19A":"02B","06B":"01C","06C":"03D","06D":"03C","06E":"04C","06F":"05D","06G":"05C","06H":"07B","08B":"09G","03B":"09F","04D":"09E","04E":"09D","05B":"09C","07C":"09B"}
            self.connect=[("04A","04B"),("02D","02E"),("06C","06D"),("03C","03D"),("06F","06G"),("05C","05D"),("04D","04E"),("09D","09E")]
            self.sticky=["03A","05A","03F","05F"]
            if self.type=="I":
                self.templatefile="I'newtemplate%d-%d.svg"%(fgsn,tvsn)
            else:
                self.templatefile="A'newtemplate%d-%d.svg"%(fgsn,tvsn)
        self.stringNums=[1,2,3,4,5,6,7,8,9,19]
        if typ=='xj': 
	  self.templatefile="A'newtemplate%dxj.svg"%fgsn

        

    def getIDs(self):#TODO:what about I???
        """identifies SEED tile strings according to the seed ID and
        writes them into the stringIDs dictionary"""
        for i in [1, 2, 6, 7, 8, 9, 19]:
            self.stringIDs[i]="tw-u2AB%d-%02d"%(fgsn, i)
        self.stringIDs[3]="tw-%s%d-03"%(2*self.type, tvsn)
        self.stringIDs[4]="tw-ICEN-04"
        self.stringIDs[5]="tw-%s%d-05"%(2*self.type, tvsn)

        if self.type=="I":
            for i in [1, 6, 8]:
                #self.stringIDs[i]="tw-u2AB%dB-%02d"%(fgsn, i) #TODO: missing!!!!!
                self.stringIDs[i]="tw-u2AB%d-%02d"%(fgsn, i) #TODO: missing!!!!!
            self.stringIDs[3]="tw-ICN%d-03"%tvsn
            self.stringIDs[5]="tw-ICN%d-05"%tvsn
        
        if self.type=="B":
            #self.stringIDs[3]="tw-u3B%d-03"%tvsn
            #self.stringIDs[4]="tw-u3AB-04"
            #self.stringIDs[5]="tw-u3B%d-05"%tvsn            
	    self.stringIDs[3]="tw-%s%d-03"%(2*self.type, tvsn)
	    self.stringIDs[4]="tw-ICEN-04"
	    self.stringIDs[5]="tw-%s%d-05"%(2*self.type, tvsn)
        if typ=='xj': 
	  for i in [3,4,5,8,19]:self.stringIDs[i]=self.stringIDs[i]+'xj'
	  self.stringIDs[2]="tw-u2%s%d-02xj"%(2*self.type, fgsn)
	  self.stringIDs[9]="tw-u2%s%d-09xj"%(2*self.type, fgsn)

    def fetchBody(self):
        body={}
        if self.type=="B":
            bfile=open(testdir+"blunt-fgen-hp-new%d.txt"%fgsn)
        else:
            bfile=open(testdir+"blunt-fgen-new%d.txt"%fgsn)
        for line in bfile:
            body[line.split("\t")[0]]=line.split("\t")[1].replace('\n','')
        bfile.close()
        return body


    def fetchAllSticky(self,charlen):
        allsticky={}
        strandsd=allstrands(allstrandfile)
        for key in strandsd.keys():
            strand=strandsd[key]
            strand=strand.replace(" ","")
            strand=strand.replace("5'-","")
            strand=strand.replace("-3'","")
            strandsd[key]=strand

        for tile in ["A1", "B2","B3", "A4", "B5", "A6", "B7"]:
            for i in [10,14,15,16,17,18]:
                thisid="tw-%s-%02d"%(tile,i)
                try:
                    allsticky[thisid]=strandsd[thisid][:16+charlen-1]
                    if thisid=="tw-B7-18":
                        allsticky[thisid+"e"]=strandsd[thisid][-(6+charlen-1):][::-1]
                    if allsticky[thisid][:6]=="TTTTTT":
                        allsticky[thisid]=allsticky[thisid][:6+charlen-1]
                except KeyError:
                    pass

            if model=='new':#why doesn't this register?
                for i in [10,14,15,16,17,18]:
                    thisid="tw-%s-%02dn"%(tile,i)#use a newer version if possible (what to do about dated strands?)
                    try:
                        allsticky[thisid]=strandsd[thisid][:16+charlen-1]
                        dummy=allsticky.pop(thisid[:-1])
                        if thisid=="tw-B7-18":
                            allsticky[thisid+"e"]=strandsd[thisid][-(6+charlen-1):][::-1]
                        if allsticky[thisid][:6]=="TTTTTT":
                            allsticky[thisid]=allsticky[thisid][:6+charlen-1]
                    except KeyError:
                        pass

            for thisid in ["tw-ICN%d-11"%tvsn,"tw-ICN%d-13"%tvsn, "tw-AA%d-11"%tvsn, "tw-AA%d-13"%tvsn,"tw-BB%d-11"%tvsn, "tw-BB%d-13"%tvsn,"tw-ICN%d-03"%tvsn,"tw-ICN%d-05"%tvsn, "tw-AA%d-03"%tvsn,"tw-AA%d-05"%tvsn,"tw-u3B%d-03"%tvsn,"tw-u3B%d-05"%tvsn, "tw-BB%d-03"%tvsn, "tw-BB%d-05"%tvsn ]:
                allsticky[thisid+"b"]=strandsd[thisid][:7+charlen-1]
                allsticky[thisid+"e"]=strandsd[thisid][-(7+charlen-1):][::-1]

            for thisid in ["tw-u2AB%d-%02d"%(fgsn, i) for i in [1,2,6,7,8,9,10,14,15,16,17,18]]:
		try: allsticky[thisid]=strandsd[thisid][:fgsn+charlen-1]
		except KeyError: print 'key not found: '+thisid

        return allsticky


class secgentile(seedtile):
    def __init__(self,newID):
        seedtile.__init__(self,newID)
        del self.num
        if self.type=="B":
            self.matches={"10B":"15B", "10C":"11B", "10D":"12E", "10E":"12D", "10F":"13B", "10G":"14A", "15A":"16H", "11C":"16G", "11D":"16F", "12C":"16E", "13C":"16D", "13D":"16C", "18B":"16B", "17B":"20A", "17C":"11E",  "17D":"12B", "17E":"12A", "17F":"13E",  "17G":"18A"}
            self.connect=[("10D","10E"),("12D","12E"), ("11C","11D"),("16F","16G"), ("13C","13D"),("16C","16D"),("17D","17E"),("12A","12B")]
            self.sticky=["11A","13A","11J","13J"]
            self.templatefile="B''template.svg"
        else:
            self.matches={"10B":"15B", "10C":"11B", "10D":"12M", "10E":"12H", "10F":"13B", "10G":"14A", "15A":"16H", "11C":"16G", "11H":"16F", "12G":"16E", "13C":"16D", "13H":"16C", "18B":"16B", "17B":"20A", "17C":"11I", "17D":"12F", "17E":"12A", "17F":"13I", "17G":"18A", "11D":"11F", "12I":"12K", "13D":"13F", "12B":"12D"}
            self.sticky=["11A","13A","11J","13J"]
            self.connect=[]
            if self.type=="I": self.templatefile="I''template.svg"
            else: self.templatefile="A''template.svg"
        self.stringNums=[10,11,12,13,14,15,16,17,18,20]

    def getIDs(self):
        """identifies SEED tile strings according to the seed ID and
        writes them into the stringIDs dictionary. This had better work"""
        for i in [10, 14, 15, 16, 17, 18, 20]:
            self.stringIDs[10]="tw-u2AB-%02d"%i
            
        self.stringIDs[11]="tw-u3A%d-11"%tvsn
        self.stringIDs[12]="tw-u3A-12"
        self.stringIDs[13]="tw-u3A%d-13"%tvsn

        if self.type=="I":
            self.stringIDs[10]="tw-u2I-10"
            self.stringIDs[11]="tw-u3I%d-11"%tvsn
            self.stringIDs[12]="tw-u3AB-12"
            self.stringIDs[13]="tw-u3I%d-13"%tvsn
            self.stringIDs[16]="tw-u2I-16"
            self.stringIDs[17]="tw-u2I-17"
            self.stringIDs[20]="tw-ICEN-20"

        if self.type=="B":
            self.stringIDs[11]="tw-BB%d-11"%tvsn
            self.stringIDs[12]="tw-ICEN-12"
            self.stringIDs[13]="tw-BB%d-13"%tvsn

    def fetchBody(self):
        body={}
        if self.type=="B":
            bfile=open(testdir+"blunt-sg.txt")
        else:
            bfile=open(testdir+"blunt-sg-hp.txt")
        for line in bfile:
            body[line.split("\t")[0]]=line.split("\t")[1].replace('\n','')
        return body

class newsgtile(seedtile):
    def __init__(self,newID):
        seedtile.__init__(self,newID)
        del self.num
        if self.type=="B":
            self.matches={"10B":"15C", "10C":"11B", "10D":"12E", "10E":"12D", "10F":"13B", "10G":"14B", "15B":"16H", "11C":"16G", "11D":"16F", "12C":"16E", "13C":"16D", "13D":"16C", "18C":"16B", "17B":"20A", "17C":"11E",  "17D":"12B", "17E":"12A", "17F":"13E",  "17G":"18B"}
            self.connect=[("10D","10E"),("12D","12E"), ("11C","11D"),("16F","16G"), ("13C","13D"),("16C","16D"),("17D","17E"),("12A","12B")]
            self.sticky=["11A","13A","11J","13J"]
            self.templatefile="B''newtemplate%d-%d.svg"%(fgsn,tvsn)
            #if fgsn==9: self.templatefile="B''newtemplate9.svg"
        else:
            self.matches={"10B":"15C", "10C":"11B", "10D":"12M", "10E":"12H", "10F":"13B", "10G":"14B", "15B":"16H", "11C":"16G", "11H":"16F", "12G":"16E", "13C":"16D", "13H":"16C", "18C":"16B", "17B":"20A", "17C":"11I", "17D":"12F", "17E":"12A", "17F":"13I", "17G":"18B", "11D":"11F", "12I":"12K", "13D":"13F", "12B":"12D"}
            self.sticky=["11A","13A","11J","13J"]
            self.connect=[]
            if self.type=="I":
                self.templatefile="I''newtemplate%d-%d.svg"%(fgsn,tvsn)
                #if fgsn==9: self.templatefile="I''newtemplate9.svg"
            else:
                self.templatefile="A''newtemplate%d-%d.svg"%(fgsn,tvsn)
                #if fgsn==9: self.templatefile="A''newtemplate9.svg"
        self.stringNums=[10,11,12,13,14,15,16,17,18,20]

    def getIDs(self):
        """identifies SEED tile strings according to the seed ID and
        writes them into the stringIDs dictionary. This had better work"""
        
        for i in [10, 14, 15, 16, 17, 18, 20]:
            self.stringIDs[i]="tw-u2AB%d-%02d"%(fgsn, i)
        
        self.stringIDs[11]="tw-u3A%d-11"%tvsn
        self.stringIDs[12]="tw-u3AB-12"
        self.stringIDs[13]="tw-u3A%d-13"%tvsn

        if self.type=="B":
            self.stringIDs[11]="tw-BB%d-11"%tvsn
            self.stringIDs[12]="tw-ICEN-12"
            self.stringIDs[13]="tw-BB%d-13"%tvsn

        if self.type=="I":
            self.stringIDs[10]="tw-u2I%d-10"%fgsn
            self.stringIDs[11]="tw-u3I%d-11"%tvsn
            self.stringIDs[12]="tw-u3AB-12"
            self.stringIDs[13]="tw-u3I%d-13"%tvsn
            self.stringIDs[15]="tw-u2I%d-15"%fgsn
            self.stringIDs[17]="tw-u2I%d-17"%fgsn

    def fetchBody(self):
        body={}
        if self.type=="B":
            bfile=open(testdir+"blunt-sgen-hp-new%d.txt"%fgsn)
        else:
            bfile=open(testdir+"blunt-sgen-new%d.txt"%fgsn)
        for line in bfile:
            if len(line)>2: body[line.split("\t")[0]]=line.split("\t")[1].replace('\n','')
        return body

    def fetchAllSticky(self,charlen):
        allsticky={}
        strandsd=allstrands(allstrandfile)
        for key in strandsd.keys():
            strand=strandsd[key]
            strand=strand.replace(" ","")
            strand=strand.replace("5'-","")
            strand=strand.replace("-3'","")
            strandsd[key]=strand

        for tile in ["A1", "B2","B3", "A4", "B5", "A6", "B7"]:
            for i in [10,14,15,16,17,18]:
                thisid="tw-%s-%02d"%(tile,i)
                try:
                    allsticky[thisid]=strandsd[thisid][:16+charlen-1]
                    if thisid=="tw-B7-18":
                        allsticky[thisid+"e"]=strandsd[thisid][-(6+charlen-1):][::-1]
                    if allsticky[thisid][:6]=="TTTTTT":
                        allsticky[thisid]=allsticky[thisid][:6+charlen-1]
                except KeyError:
                    pass

            for thisid in ["tw-ICN4-11","tw-ICN4-13", "tw-AA2-11", "tw-AA2-13","tw-BB2-11", "tw-BB2-13","tw-ICN4-03","tw-ICN4-05", "tw-AA2-03","tw-AA2-05","tw-u3B-03","tw-u3B-05"]:
                allsticky[thisid+"b"]=strandsd[thisid][:tvsn+charlen-1]
                allsticky[thisid+"e"]=strandsd[thisid][-(tvsn+charlen-1):][::-1]

            for thisid in ["tw-u2ABn-01","tw-u2ABn-02","tw-u2ABn-06","tw-u2ABn-07","tw-u2ABn-08","tw-u2ABn-09","tw-u2ABn-10","tw-u2ABn-14","tw-u2ABn-15","tw-u2ABn-16","tw-u2ABn-17","tw-u2ABn-18"]:
                allsticky[thisid]=strandsd[thisid][:fgsn+charlen-1]

        return allsticky




if __name__ == "__main__":
	pass

#    repl={"tw-BB8-11b":"AGTCAAGC", "tw-BB8-11e":"GCAGAGTT", "tw-BB8-13b":'GAGGTTCC', "tw-BB8-13e":'ACCTGTCT', "tw-AA8-11b":'CTCCTACG', 'tw-AA8-11e':'GGTCCTTC', 'tw-AA8-13b':'TCTCTGCT', 'tw-AA8-13e':'TCCACCTT' }
#    repl={"tw-BB8-11b":"AGTCAAGC", "tw-BB8-11e":"GCAGAGTT", "tw-BB8-13b":'GAGGTTCC', 'tw-BB8-13e':'CTTCACGA',"tw-AA8-11b":'CTCCTACG', 'tw-AA8-11e':'GGTCCTTC', 'tw-AA8-13b':'TCTCTGCT', 'tw-AA8-13e':'TCCACCTT', 'tw-ICN8-11b':'ATGACAGC',  'tw-ICN8-11e':'GACTCCAG', 'tw-ICN8-13b':'TCAAGACG', 'tw-ICN8-13e':'AGAGCAGA'}melting -Hdnadna -SAGTCAAGC -CTCAGTTCG -G12.5e-3 -P0.5e-7

    ##repl={'tw-u2AB7-01':'CGTCTTA','tw-u2AB7-02':'TAAXACG', 'tw-u2AB7-06':'TATCACG', 'tw-u2AB7-07':'CGTXATA', 'tw-u2AB7-08':'ATTGGAC', 'tw-u2AB7-09':'GTCXAAT'}#1:14.81, 6:14.99, 8:14.38#wrong ones
    #repl={'tw-u2AB7-01':'GTTCTCA','tw-u2AB7-02':'TCAXAAC', 'tw-u2AB7-06':'TCACTTC', 'tw-u2AB7-07':'GAAXTGA', 'tw-u2AB7-08':'TTCCTTG', 'tw-u2AB7-09':'CAAXGAA'}#1:14.38, 6:14.38, 8:14.48
    #repl={'tw-u2AB5-01':'TCGTG','tw-u2AB5-02':'CAXGA', 'tw-u2AB5-06':'CGGTT', 'tw-u2AB5-07':'AAXCG', 'tw-u2AB5-08':'ACGTC', 'tw-u2AB5-09':'GAXGT'}#1:-1.28, 6:-1.73, 8:14.38
    ###repl={'tw-u2AB5-01':'TTCGG','tw-u2AB5-02':'CCXAA', 'tw-u2AB5-06':'GTCGA', 'tw-u2AB5-07':'TCXAC', 'tw-u2AB5-08':'ATGGC', 'tw-u2AB5-09':'GCXAT'}#1:-2.98,6:-2.67, 8:-3.97#wrong ones
    #fgsn=7

#    repl={"tw-A6-15":"ATGCCTCGCTCAATTG",
#      "tw-B5-16":"CAATTGAGCGAGGCAT",
#      "tw-A6-16":"ATCGTCATGTAGCACC",#output100922c.txt
#      "tw-B7-15":"GGTGCTACATGACGAT",
#      "tw-A1-16":"TCGTGTTCACAACACC",#output100922d.txt
#      "tw-B2-15":"GGTGTTGTGAACACGA",
#      "tw-B3-16":"AACTTGGTTCGGAACC",#output100922e.txt
#      "tw-A4-15":"GGTTCCGAACCAAGTT",
#      "tw-B3-18":"AGGTGGCAGAACAATG",#output100922f.txt
#      "tw-A4-17":"CATTGTTCTGCCACCT",
#      "tw-A6-14":"AATACTTGCGTGCGTG",#output100922g.txt
#      "tw-B7-10":"CACGCACGCAAGTATT"

#    repl={"tw-A6-15":"ATGCCTCGCTCAATTG",
#    "tw-B5-16":"CAATTGAGCGAGGCAT",
#     "tw-A6-16":"ATCGTCATGTAGCACC",#output100922c.txt
#     "tw-B7-15":"GGTGCTACATGACGAT",
#     "tw-A1-16":"TCGTGTTCACAACACC",#output100922d.txt
#     "tw-B2-15":"GGTGTTGTGAACACGA",
#     "tw-B3-16":"AACTTGGTTCGGAACC",#output100922e.txt
#     "tw-A4-15":"GGTTCCGAACCAAGTT",
#     "tw-B3-18":"AGGTGGCAGAACAATG",#output100922f.txt
#     "tw-A4-17":"CATTGTTCTGCCACCT",
#     "tw-A6-14":"AATACTTGCGTGCGTG",#output100922g.txt
#     "tw-B7-10":"CACGCACGCAAGTATT"}

  #repl={
  #"tw-BB8-11b":"AGTAATCT", "tw-BB8-13b":"GAGAATAA",
  #"tw-ICN8-11b":"TACCTTAA", "tw-ICN8-13b":"TATAATGC",
  #"tw-AA8-11b":"GTAATTAC", "tw-AA8-13b":"TAGATTGA",
  #"tw-BB8-03b":"AGATTACT", "tw-BB8-05b":"TTATTCTC",
  #"tw-AA8-03b":"GTAATTAC", "tw-AA8-05b":"TCAATCTA",
  #"tw-ICN8-03b":"TTAAGGTA", "tw-ICN8-05b":"GCATTATA",
  #"tw-BB8-11e":"ATTCCATA", "tw-BB8-13e":"ACTAAGTA",
  #"tw-ICN8-11e":"TATGTGTA", "tw-ICN8-13e":"AATAGTCA",
  #"tw-AA8-11e":"TATTGGAA", "tw-AA8-13e":"CAAGTTAA",
  #"tw-BB8-03e":"TATGGAAT", "tw-BB8-05e":"TACTTAGT",
  #"tw-AA8-03e":"TTCCAATA", "tw-AA8-05e":"TTAACTTG",
  #"tw-ICN8-03e":"TACACATA", "tw-ICN8-05e":"TGACTATT"
  #}
  #thisrep="tw-ICN8-05e"
  #fgsn=7
  #tvsn=8
  #apr=newfgtile("B'") 
  #cand=[]
  #acand={}
  #for i in range(100):
    #new1=noRepetitions(apr.ReplaceSticky(thisrep,8, bm=7, sm=6,  nlim=80, tlim=80, verbose=False, replacements=repl),3)
  #if len(new1)>0:  
		#cand=TMeltSelect(new1,(0,20),verbose=True, threads=4)
		#acand=dict(acand, **cand)
  #else: print "no strands found"
  #f=localtime()
  #fil=open('/windows/D/NY/randomgel/sticky%s-%02d-%02d-%02d.%02d.log'%(str(f.tm_year)[-2:],f.tm_mon,f.tm_mday,f.tm_hour,f.tm_min),'w')
  #f=open('scripts/cand%s-%02d-%02d-%02d.%02d.%02d.log'%(str(f.tm_year)[-2:],f.tm_mon,f.tm_mday,f.tm_hour,f.tm_min,f.tm_sec),'w')
  #for item in repl.keys(): f.write('%s\t%s\n'%(item,repl[item]))
  #f.write('\n################%s################\n\n'%thisrep)
  #for item in acand.keys(): f.write('%s\t%.02f\n'%(item,acand[item]))
  #f.close()


#########################stuff that shouldn't be in this module################

def nobranches(grey,cut,linelen):
    bbin=grey>cut
    blab=p.label(1-bbin)
    blobs=p.blob(blab,'boundingbox','data')
    nums={}
    Image.fromarray(blab.astype(uint8)*(255/amax(blab))).show()
    for num in arange(1,amax(blab)+1):
        nums[num]=[]
    for ang in arange(6)*30:
        f=p.closeth(bbin,p.seline(linelen,ang))
        ad= uniq(blab[(f==bbin).nonzero()])
        for key in nums.keys():
            if key in ad: nums[key].append(ang)
    for key in nums.keys():
        if logical_or(len(nums[key])==0, len(nums[key])>2):
            blab[(blab==key).nonzero()]=0
            print "Blob removed: %d"%key
    Image.fromarray(blab.astype(uint8)*(255/amax(blab))).show()
    bbin=blab>0
    nblab=p.label(bbin)
    blobs=p.blob(nblab,'boundingbox','data')
    return (blobs,nums,blab)

class AFMpic:
    """Does things with AFM ascii exports"""
    def __init__(self,fname):
        self.file=fname
        self.size=0
        self.data=array([])
        self.grey=array([])

    def readData(self):
        f=open(self.file)
        b=f.read()
        f.close()
        c=b.find('\Scan size:')
        d=b.find(' nm',c)
        self.size=float(b[c+11:d])
        c=b.find('\Exported image units: nm')
        if c>-1:
            temp=b[c+26:]
            self.lengths=True
        else:
            c=b.find('\*File list end\r\n')
            temp=b[c+16:]
            if verbose: print key
            self.lengths=False
        f=open('temp','w')
        f.write(temp)
        f.close()
        self.data=loadtxt('temp')
        self.grey=((self.data-amin(self.data))*255/(amax(self.data)-amin(self.data))).astype(uint8)
        #u=Image.fromarray(self.grey)
        #u.show()

