#!/usr/bin/python

"""DNA structure object and helper functions.
Can be used as a standalone program to populate a SVG template."""

#Last Update 2017 10 10 by Matan Yah Ben Zion

import re,os,sys,time
from random import shuffle, random
from argparse import ArgumentParser,RawTextHelpFormatter
from numpy import array,log,ceil,arange, append
import string
from itertools import permutations

strandfileDefault='Tri2TUstrands.csv' #CSV file name with DNA sequence dictionary
strandnumsDefault='Tri2TUstrandNums.csv' #CSV file name with DNA sequence/SVG strand number map
templateDefault='Tri2TUtemplate.svg' #SVG template file name
twistfile='twists.csv' #CSV file with next-neighbour twist data


alphabet='ABCDEFGHIJKLMNOPQRSTUVWXYZ'

# Path to inkscape executable (Win/Mac only)
if sys.platform=="win32":
    inkscape_path='C:\\Program Files\\Inkscape\\'
if sys.platform=="darwin":
    inkscape_path="/Applications/Inkscape.app/Contents/Resources/bin/"
if sys.platform.startswith("linux"):
    inkscape_path=""

def complement(dna, reverse=False, verbose=True):
    """Returns the complementary string to a 'GATC' type sequence string.
    Alien bases are left unchanged, the reverse complement is returned with the
    reverse=True keyword."""
    cdir={'G':'C','C':'G','A':'T','T':'A'}
    comp=''
    for c in dna.upper():
        try: comp=comp+cdir[c]
        except KeyError:
            if verbose:
                print "Warning: non GATC-base!"
            comp=comp+c
    if reverse: comp=comp[::-1]
    return comp
  
def CG(seq):
  """CG content of a given sequence, 0<=CG<=1. Input string, returns float"""
  return float(seq.count('C')+seq.count('G'))/len(seq)

def twist(seq,twistfile=twistfile, periodic=False):
    """Returns total twist (in degrees) for a DNA sequence string.
    Change twist data file with twistfile keyword; the periodic keyword
    appends the twist between the last and first base for periodic structures.
    Default twist values from PNAS 107(35):15421-6 (2010), SI."""
    twists=dictfromfile(twistfile, valtype='float')
    totaltwist=0
    for i in range(len(seq)-1):
        pair=seq[i:i+2]
        totaltwist=totaltwist+twists[pair]
    if periodic:
        totaltwist=totaltwist++twists[seq[-1]+seq[0]]
    return totaltwist
    
def randCGconstr(fname, prob=0.6, oname=None):
    if not oname: oname=fname
    if prob < 0 or prob>1:
        print 'Please enter a CG content between 0 and 1.'
    else:
        if prob<.5: 
            base='W'
            prob=2-2*prob
        if prob>.5: 
            prob=2*prob-1
            base='S'
        with open(fname,'r') as f:
            pos=0
            txt=f.read()
        while pos>-1:
            r=random()
            pos=txt.find('X',pos+1)
            if r<prob: txt=txt[:pos]+base+txt[pos+1:]
        with open(oname,'w') as f: f.write(txt)

                     

def dictfromfile(fname, keytype='str', valtype='str',inds=(0,-1), quotech='',verbose=False):
    """Returns a dictionary from a comma/whitespace separated text file.
    If keytype and valtype keywords are present (values int|str|float),
    string values will be converted to a given type.
    If key/value pairs don't correspond to the first and last columns,
    change the inds keyword from (0,-1).
    For spread sheet exporting text in quotes, strip them by setting the quotech keyword.
    verbose output prints the type conversion syntax"""
    fdict={} #initialise empty dictionary
    with open(fname,'r') as f:
        for line in f:
            try:
            #split line by commas/whitespace, use predefined colums (inds tuple) and type convert
            #via exec statement using type string. Strip quotation marks first.
                n="a=%s('%s')"%(keytype, line.split()[inds[0]].replace(quotech,'')); exec n
                m="b=%s('%s')"%(valtype, line.split()[inds[-1]].replace(quotech,'')); exec m
                if verbose: print 'key: ',n; print 'value: ',m
                fdict[a]=b #setting new dictionary entry
            except: pass #if it fails, move on (useful if file ends with empty line)
    return fdict

def reformat(s,c, colour="#000000", weight="bold"):
    """Helper function to highlight out non-matching bases.
    Set colour or weight to None if you don't want to use them for highlighting.
    Returns the complete string if lengths are inconsistent, but mismatch highlighting will probably be incorrect."""
    c=complement(c, reverse=True)
    news=''
    if colour and not weight: style="fill:%s"%colour
    if weight and not colour: style="font-weight:%s"%weight
    if weight and colour: style="font-weight:%s;fill:%s"%(weight, colour)
    if not weight and not colour: return s
    for i in range(len(s)):
        try:
            if c[i]!=s[i]: news=news+'<tspan style="%s">%s</tspan>'%(style, s[i])
            else: news=news+s[i]
        except IndexError:
            news=news+s[i:]
            break
    return news

def d2b(num, base, dd=False):
    if not 2 <= base <= 36:
        raise ValueError, 'The base number must be between 2 and 36.'
    if not dd:
        dd = dict(zip(range(36), list(string.digits+string.ascii_lowercase)))
    if num == 0: return ''
    num, rem = divmod(num, base)
    return d2b(num, base, dd)+dd[rem]
    
def pseqAug(pseq,bpath=None):
    seq=''.join(pseq.astype(str))
    l=len(seq)
    z=int(seq,4)+1
    if z<4**l:
        new=array(list(d2b(z,4).zfill(l))).astype(int)
        new=new[:(pseq==new).tolist().index(False)+1]
        
        try:
            bs=(bpath[:len(new),0]==bpath[:len(new),1]).nonzero()[0].tolist()
            for i in bs:
                seq=''.join(pseq.astype(str))
                if new[i]!=seq[i]:
                    j=i
                    while j in bs:
                        bs.remove(j)
                        j-=1
                    if len(bs)>0:
                        j=bs[-1]
                        z=int(''.join(new[:bs[-1]]),4)+1
                        if z>4**bs[-1]:
                            print "End of line!"
                            return -1
                        else:
                            new=array(list(d2b(z,4).zfill(l))).astype(int)+new[bs[-1]:]
        except: pass
    else:
        print "End of line!"
        return -1
    return new
  

                
                

def word_check(seq1,seq2,word):
  """Returns False and aborts if seq2 contains a substring of seq1 of length word. Returns True otherwise"""
  for i in range(len(seq1)-word+1):
    if seq2.find(seq1[i:i+word])>-1: return seq2.find(seq1[i:i+word])
  return -1

def noGCrep(seq,rep=4):
  #checks for CG repetitions (CCC,GGG,CGC etc.) etc of length rep
  seq = seq.upper()
  rx = r'[GC]{'+str(rep)+',}'
  return re.search(rx, seq) == None
  

"""def noCGrep(seq,rep):
  #checks for CG repetitions (CCC,GGG,CGC etc.) etc of length rep
  cand=[]
  seq = seq.upper()
  for i in range(rep+1):
    cand=cand+list(set(map(''.join, permutations((rep-i)*'C'+i*'G'))))
  for item in cand:
    if seq.find(item)>-1: return False
  return True """
    
  
def palin_check(seq, num):
  seq=seq.upper()
  for i in range(len(seq)+1-num):
    b=complement(seq[i:i+num], reverse=True)
    f=seq.find(b)
    if seq.find(b)>-1 and f!=i: return i
  return -1
    
    
def generate_path(length, constraints=False):
    path=array([['G','A','T','C']]*length)
    for s in path: shuffle(s)
    if constraints:
        for i in range(len(constraints)):
            try:
                if constraints[i]=='S': f=['C', 'C', 'G', 'G']#influence CG content
                if constraints[i]=='W': f=['A', 'A', 'T', 'T']
                #add contraints for purines and pyrimidines
                if constraints[i] in ['G','A','T','C']: f=[constraints[i]]*4
                else: f=['G','A','T','C']
                shuffle(f)
                path[i,:]=f
            except IndexError: break
    return path
  
def rand_seq(length, cg=(.4,.6),tmax=100):
  seq=''
  bases=['G','A','T','C']
  t0=time.time()
  while len(seq)<length and (time.time()-t0)<tmax:
    shuffle(bases)
    seq=seq+bases[0]
    if len(seq)==length and not (cg[0] <CG(seq)<cg[1]): seq=''
  return seq

def reptest(string,num):
    string = string.upper()
    if string.find('A'*num)>=0 or string.find('C'*num)>=0 or string.find('G'*num)>=0 or string.find('T'*num)>=0: 
        return True
    else: return False

    
class MyError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)    

class DNAstrand():
    def __init__(self, ident, length, seq='', constraints=False):
        self.id=ident
        self.seq=seq
        self.bpath=generate_path(length, constraints)        
        self.conn={}
        self.chunks={}
        self.ivchunks={}
       
       
    
    
class DNAStructure():
    """DNA structure object. Contains strand dictionary, strand map, SVG structural template.
    Set template file, strand dictionary file and SVG strand map file with temp, sts and nums keywords."""
    def __init__(self, temp=templateDefault, sts=strandfileDefault, nums=strandnumsDefault, stcols=(0,-1)):
        self.strands=dictfromfile(sts,inds=stcols)
        self.strandnums=dictfromfile(nums, keytype='int')
        self.strandfile=sts
        self.numfile=nums
        self.template=temp
        self.matchmap=self.readMatchmap()
        self.weight=sum([len(i) for i in self.strands.values()])
        self.strobs={}
    
    def lookCbase(self,sid,idx):
        try: 
            b=self.strobs[sid].conn[idx]
            return self.strobs[b[0]].seq[b[1]]
        except: 
            return -1
        
    def readTemplate(self,  init=True, verbose=False):
        """reconstructs strand structure and connectivity from template file"""
        if not init: self.strandfile=self.template[:-4]+'-seq.csv'
        if not init: self.numfile=self.template[:-4]+'-num.csv'
        mmap=array(self.readMatchmap())
        self.mmap=mmap
        if not init: self.strands={}
        if not init: self.strandnums={}
        with open(self.template, 'r') as f:
            svg=f.read()
        if not init: r=range(1, 101)
        else: r=self.strandnums.keys()
        for i in r:
            if not init: id='strand%02d'%i
            else: id=self.strandnums[i]
            chunks=re.findall('(?<=>)%02d[A-Z]X*?(?=<)'%i,svg)
            if verbose: print chunks
            if len(chunks)>0:
                chunks.sort()
                seq=''.join(chunks)
                self.strobs[id]=DNAstrand(id,len(seq))
                self.strobs[id].seq='X'*len(seq)
                self.strandnums[i]=id
                self.strands[id]='X'*len(seq)
                count=0
                for chunk in chunks:
                    self.strobs[id].chunks[chunk[:3]]=(count, len(chunk))
                    self.strobs[id].ivchunks[count]=chunk[:3]
                    count=count+len(chunk)
            else: pass
        for sid in self.strobs:
            for ch in self.strobs[sid].chunks.keys():
                st,l=self.strobs[sid].chunks[ch]
                idx=(mmap==ch[:3]).nonzero()
                if len(idx[1])==0: pass                
                elif idx[1][0]==0:
                    cch=mmap[idx[0][0],1] 
                    if not init: id='strand%02d'%cch[:2]
                    else: id=self.strandnums[int(cch[:2])]
                    cst,cl=self.strobs[id].chunks[cch]
                    for i in range(l):
                        self.strobs[sid].conn[st+i]=(id,cst+min(cl,mmap[idx[0][0],3])-i-1)
                elif idx[1][0]==1:
                    cch=mmap[idx[0][0],0]
                    if not init: id='strand%02d'%cch[:2]
                    else: id=self.strandnums[int(cch[:2])]
                    cst,cl=self.strobs[id].chunks[cch]
                    for i in range(len(('x'*l)[int(mmap[idx[0][0],2]):int(mmap[idx[0][0],3])])):
                        self.strobs[sid].conn[st+int(mmap[idx[0][0],2])+i]=(id,cst+cl-i-1)
        self.weight=sum([len(i) for i in self.strands.values()])
                
            
    def populate(self, word=0, repeat=4, constrfile=None, lim_t=300000, cg=.5, verbose=False):
        sts=self.strandfile
        nums=self.numfile
        self.readTemplate(init=True)
        if sts!='': self.strandfile=sts
        if nums!='': self.numfile=nums
        if verbose: print 'hello'
        t1=time.time()
        if word==0: 
            word=int(ceil(log(self.weight)/log(4))+1)
            print 'word length: %d'%word
        if constrfile:
            if cg!=.5:
                randCGconstr(constrfile,cg,oname='temp.txt')
                const=dictfromfile('temp.txt')
                os.remove('temp.txt')
            else: const=dictfromfile(constrfile)
        else: 
            const={key:'X'*len(self.strands[key])}# for key in self.strands.keys()}
        for key in self.strands.keys():
            self.strobs[key].seq=self.strands[key]
            self.strandnums=dictfromfile(self.numfile,keytype='int')
            self.strands=dictfromfile(self.strandfile)
        for key in const.keys():
            self.strobs[key].bpath=generate_path(len(self.strobs[key].seq), constraints=const[key])
            self.strands[key]=const[key]+self.strands[key][len(const[key]):]
            self.strobs[key].seq=self.strands[key]
        dnas=self.strobs.values()
        while len(dnas)>0:
            shuffle(dnas)
            startover=False
            for dna in dnas:
                if dna.seq.find('X')==-1:
                    if verbose: print dna.id,  dna.seq,  'excluded'
                    dnas.remove(dna)
                    break
                pseq=array([0]*word)
                for i in range(len(dna.seq)):
                    b=self.lookCbase(dna.id,i)
                    if b in ['G','A','T','C']:
                        dna.bpath[i,:]=[complement(b)]*4
                        if i<word: pseq[i]=3
                test=''.join(dna.bpath[range(len(pseq)),pseq])
                while reptest(test,repeat) and type(pseq).__name__=='ndarray':
                    pseq=pseqAug(pseq,dna.bpath)
                test=''.join(dna.bpath[range(len(pseq)),pseq])
                dna.seq=test+dna.seq[word:]
                self.strands[dna.id]=dna.seq
                for i in range(word):
                    try:
                        s=self.strobs[dna.conn[i][0]].seq 
                        j=self.strobs[dna.conn[i][1]]
                        self.strobs[dna.conn[i][0]].seq=s[:j]+complement(test[i])+s[j+1:]
                        self.strands[dna.conn[i][0]]=self.strobs[dna.conn[i][0]].seq
                    except KeyError: pass
                if verbose: print dna.seq
                while set(dna.seq)!=set(['A','T','G','C']):
                    for key in self.strands.keys():
                        find=self.strands[key].find(complement(test[-word:], reverse=True))
                        if time.time()-t1 > lim_t: 
                                    print "Time limit exceeded!"
                                    return -1
                        if find>0:
                            try: 
                                if not (dna.conn[len(test)-word]==(key, find+word-1) and dna.conn[len(test)-1]==(key, find) 
                                and self.strands[key].find(complement(test[-word:], reverse=True),find+1)<0): #random match?
                                    raise KeyError #trick to avoid KeyError, case is handled in exception
                                else: pass #if complement is allowed, move on.
                                
                            except KeyError:
                                pseq=pseqAug(pseq, dna.bpath) #go to next possible path
                                while reptest(test,repeat) and type(pseq).__name__=='ndarray':
                                    pseq=pseqAug(pseq,dna.bpath)
                                    try: test=''.join(dna.bpath[range(len(pseq)),pseq])
                                    except KeyError:
                                      print 'KeyError!'
                                      return -1
                                if type(pseq).__name__!='ndarray':
                                    startover=True #all possible paths explored, discard this strand and one complement! 
                                try: test=''.join(dna.bpath[range(len(pseq)),pseq])
                                except TypeError:
                                  print 'Type error', sys.exc_info()[0]
                                  return -1
                                break #don't check the other sequences if base doesn't work
                                if time.time()-t1 > lim_t: 
                                    print "Time limit exceeded!"
                                    return -1
                                           
                    if startover: 
                        for v in dna.conn.values(): #put one of the complement strands back in the pot
                            if not self.strobs[v[0]] in dnas:
                                dnas.append(self.strobs[v[0]])
                                self.strobs[v[0]].seq='X'*len(self.strobs[v[0]].seq) #reset sequences for this and the complement
                                if v[0] in const:
                                    self.strobs[v[0]].seq=const[v[0]]+self.strobs[v[0]].seq[len(const[v[0]]):]
                                self.strands[v[0]]=self.strobs[v[0]].seq
                                dna.seq='X'*len(dna.seq)
                                if dna.id in const:
                                    dna.seq=const[dna.id]+dna.seq[len(const[dna.id]):]
                                self.strands[dna.id]=dna.seq
                                break
                        break
                    dna.seq=test+dna.seq[len(test):]
                    self.strands[dna.id]=dna.seq
                    if verbose: print dna.seq
                    try:
                        b=dna.conn[len(test)-1]
                        self.strobs[b[0]].seq=self.strobs[b[0]].seq[:b[1]]+complement(test[-1])+self.strobs[b[0]].seq[b[1]+1:]
                    except KeyError: pass
                    if len(pseq)<len(dna.bpath[:,0]):    
                        pseq=append(pseq,0)
                        test=''.join(dna.bpath[range(len(pseq)),pseq])
                        if verbose: print test
                        if verbose: print pseq
                    else:
                        if verbose: print dna.id,  dna.seq,  'path'
                        dnas.remove(dna)
                        break    
                if time.time()-t1 > lim_t: 
                    print "Time limit exceeded!"
                    return -1
                if startover: break
        for d in self.strobs.values():
            self.strands[d.id]=d.seq
        print "Populated! Elapsed time %s h."%time.strftime('%H:%M:%S', time.gmtime(time.time()-t1))
        return 0
      
      
    def writeCSVs(self, fname=''):
      if fname=='':
        if os.access(self.strandfile, os.F_OK):
            os.rename(self.strandfile, self.strandfile+'.bak')
      else: fname==self.strandsfile
      with open(fname, 'w') as f:
          for s in self.strands:
              f.write('%s %s'%(s, self.strands[s]))
      if os.access(self.numfile, os.F_OK):
          os.rename(self.numfile, self.numfile+'.bak')
      with open(self.numfile, 'w') as f:
          for s in self.strandnums:
              f.write('%s %s'%(s, self.strandnums[s]))
                
    def readMatchmap(self):
        """extracts a list mapping complementarity information for different strand segments
        from a special comment in the SVG template."""
        matches=[] #returns empty list if not successful
        try:
            with open(self.template,'r') as f: #find special comment
                a=re.search('(?<=<!--matchmap).*?(?=-->)',f.read(),re.DOTALL)
                exec "matches=%s"%a.group()
                #execute comment string data (already in map form)
        except: print "No match map found in %s!"%self.template
        return matches

    def fillTemplate(self, fname, verbose=False, twistfile=twistfile, inkpath=inkscape_path, label='label'):
        """Populates SVG template file with DNA sequences with the target filename as argument.
        Checks for strand length consistency.
        If the template contains a match map, also checks for complementarity errors.
        If the template contains twist markers, computes twist values for certain segments.
        If inkscape is found, produces additional PDF file.
        Keywords: twistfile - change file containing next neighbour pair twist values;
        inkpath - path to inkscape executable."""
        try:
            with open(self.template,'r') as f: svg=f.read()
        except IOError:
            print 'Could not read %s'%self.template
            return -1 #if svg reading fails, abort
        for num in self.strandnums.keys():
            #walk through all strands in map dictionary
            end=0
            strand=self.strands[self.strandnums[num]] #strand sequence
            for alph in alphabet:
                #walk through template chunks, search for '01AXXXX' etc.
                n=re.search(">%02d%sX*<"%(num,alph), svg)
                if type(n).__name__!='SRE_Match':
                    pass #ignore [A-Z] letters not in template (search returns None)
                else:
                    n=n.group(0) # '>01AXXXX<'-like chunk
                    beg=end #move DNA string indices to next SVG segment
                    end=beg+len(n)-2
                    if svg.find('twist%02d%s'%(num,alph))>0:
                        #if template contains twist note for segment, calculate and replace
                        svg=svg.replace('twist%02d%s'%(num,alph), 'twist %.02f'%twist(strand[beg:end+1]))
                    rep=strand[beg:end]
                    for m in self.matchmap:
                    #complement check
                        if m[0]=="%02d%s"%(num,alph):
                            m[4]=strand[beg:end] #enter segment in match map
                            if type(m[5]).__name__=='str': #if complement has been entered as well, proceed
                                if m[4]!=complement(m[5], reverse=True): #complement check
                                    print "Mismatch!", m[0],m[1],m[-2],m[-1], self.strandnums[int(m[0][:2])], self.strandnums[int(m[1][:2])]
                                    rep=reformat(m[4],m[5]) #colour mismatching bases

                        if m[1]=="%02d%s"%(num,alph):
                            #same as above
                            m[5]=strand[beg:end][m[2]:m[3]] #TODO: fix this to use one number
                            if type(m[4]).__name__=='str':
                                if m[4]!=complement(m[5], reverse=True):
                                    print "Mismatch!", m[0],m[1],m[-2],m[-1], self.strandnums[int(m[0][:2])], self.strandnums[int(m[1][:2])]
                                    rep=strand[beg:end][:m[2]]+reformat(m[5],m[4])+strand[beg:end][m[3]:]
                    svg=svg.replace(n,'>%s<'%(rep)) #replace template chunk with DNA sequence
            if end !=len(strand):
                #if last template chunk is not end of DNA sequence, complain!
                print 'Length error: %d-%s, %d vs. %d'%(num, self.strandnums[num], end, len(strand))
        if svg.find('>label<')>0: svg=svg.replace('>label<','>%s<'%label)
        o=open(fname, "w") #write output SVG file
        o.write(svg)
        o.close()
        #here follows svg-> pdf conversion via inkscape, if present. If not, nothing fatal happens.
        if sys.platform=="win32":
            #Inkscape needs absolute file path on win32, Windows usually needs absolute path to inkscape.exe
            fname=os.path.abspath(fname)
            comm=inkscape_path+'inkscape.exe -z -A "'+fname[:-4]+'.pdf" "'+fname+'"'
            print comm
            os.system(comm)
        if sys.platform=="darwin":
            #Mac usually needs absolute path to inkscape as well.
            fname=os.path.abspath(fname)
            #save current directory and cd into directory with inkscape executable
            tp=os.getcwd()
            os.chdir(inkscape_path)
            #run inkscape
            os.system('./inkscape -z -A "'+fname[:-4]+'.pdf" "'+fname+'"')
            #change back to current directory
            os.chdir(tp)
        if sys.platform.startswith('linux'):
            os.system('inkscape -z -A "'+fname[:-4]+'.pdf" "'+fname+'"')
        time.sleep(0.5) #wait for inkscape to finish
        if not os.access(fname[:-4]+'.pdf', os.F_OK):
            print "PDF conversion  not successful. Did you install Inkscape?"

if __name__=='__main__':
    #makes this work as a command-line script filling a SVG template. Try 'python DNAstruct.py -h' from a non-Python terminal session to get instructions.
    myparser=ArgumentParser(description=
    """Fills a given SVG DNA structure template with sequences
    from a CSV sequence directory.  

    The SVG file must contain strand segments as text blocks 01AXXXX 01BXXXX etc.
    (for strand #1), always in 5'->3' direction.

    For complementarity testing, the SVG has to contain a comment  of the form
    "<!--matchmap [['01A','03C',0,1000,0,1],['01A','03C',0,5,0,1],...]-->",
    where the first two items in  each inner bracket are the segment IDs (shorter segment first),
    the second two are the index range of the second segment where it matches the first.
    The program will attempt to highlight mismatches by black/bold font.

    -If the template contains text blocks like 'twist01A', and a twist map is given, they will be
    replaced by the helical twist for the corresponding segment.

    -The lines of the sequence directory should be of the form 'sequenceID[whitespace/comma]sequence.'

    -There also needs to be a correspondence map file relating strand numbers from the SVG to the
    sequence directory with the line structure strand#[whitespace/comma]sequenceID.

    -The program creates an output SVG file; it also attempts PDF conversion if inkscape is found on
    the system. Adjust the path of inkscape's installation directory with -i if necessary. Make sure that
    the monospace font used in the SVG is installed.""",
    epilog=u'Have fun!',formatter_class=RawTextHelpFormatter) #formatter_class necessary for newlines in description string
    #the only positional argument can be omitted because of nargs='?'
    myparser.add_argument('filename', nargs='?', type=str, help='output *.svg file with DNA sequences', default='test.svg')
    #optional arguments. Default values are set at the top of this script.
    myparser.add_argument('-t','--templatefile', type=str, help='SVG template', default=templateDefault)
    myparser.add_argument('-n','--numbermap', type=str, help='maps strand numbers to IDs', default=strandnumsDefault)
    myparser.add_argument('-s','--stranddict', type=str, help='sequence dictionary', default=strandfileDefault)
    myparser.add_argument('-w','--twistfile', type=str, help='next neighbour twist map file', default=twistfile)
    myparser.add_argument('-i','--inkscapedir', type=str, help="inkscape's installation directory", default=inkscape_path)
    myparser.add_argument('-c','--csvcol', type=int, nargs=2, help="selects columns in sequence dictionary file, default first and last, format: 2 space separated integers", default=[0,-1])
    #read command line arguments into namespace.
    myargs = myparser.parse_args()
    print myargs
    print "HELLO!"
    #create DNA structure object.
    mystruct=DNAStructure(temp=myargs.templatefile, sts=myargs.stranddict, nums=myargs.numbermap, stcols=myargs.csvcol)
    #call fillTemplate method
    mystruct.fillTemplate(myargs.filename, twistfile=myargs.twistfile, inkpath=myargs.inkscapedir)
    #b=DNAStructure(temp='SDX.svg')
    #b=DNAStructure(temp='TXAtemplate.svg')
    #b.readTemplate()
    #b.populate()
    #b.writeCSVs()
    #b.fillTemplate('TXAfill.svg')
