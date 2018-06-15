
#------Library for manipulation of DNA Origami sequences 
# Written by Corinna Maass Oct 2013
# Updated by Matan Yah Ben Zion Sept 2014
# Updated by Matan Yah Ben Zion Aug 2015


import sys
sys.path.append('/windows/D/NY/python/modules')
from random import random, shuffle
from itertools import chain
from numpy import arange,array,concatenate
import os,re
from DNAstruct import complement
bs=16
be=1818
ms=896
me=959
es=32
ee=1791

class staple(str):
	"""A staple"""
	def __init__(self, line, delim=" ", typ="belt"):
		# self.line=line#self.__new__(self,line=str)
		dat=line[:-1].split(delim)
		self.label=dat[0]
		self.body=dat[1]
		self.se3=''
		self.se5=''
		self.stem3  = ''
		self.stem5  = ''
		self.bulge3 = ''
		self.bulge5 = ''
		self.seqStaple()

		if typ=="belt":
			self.sb=int(dat[2][1:5])
			self.sh=int(dat[2][6])
			self.eb=int(dat[2][8:12])
			self.eh=int(dat[2][13])
			self.index=dat[2][:14]
			self.num = int(self.label[2:])
			self.batch=''
		if 3<=self.num <=110: #range for belt left wing
			self.cell = int( (self.num-3) / 4)
		elif 119<=self.num <=222: #range for belt right wing
			self.cell = int( (self.num-119)/4) +27
		else: 
			self.cell=''

	def seqStaple(self):
		self.seq = self.se5 + self.stem5 + self.bulge5 + \
		self.body + \
		self.bulge3 + self.stem3 + self.se3 
		return self.seq

	def addSE(self, SE,stems,bulge3='TT', p=3):
		#stplNum = getStapleNum(stdict, cell, stplPos);
		#stpl = stdict[stplNum]	
	
		self.stem5 = stems[0]
		self.bulge3 = bulge3
		self.stem3 = stems[1]	
		self.se3 = SE

		self.seqStaple();
		
		return self
	#        if obj.num <27*32+40: shift=39

	#       else: shift=103
	#      obj.cell=(obj.num -shift)/32.
	#    return obj

	def addLabel(self, prefLabel='', sufLabel=''):
		self.label = prefLabel + self.label+sufLabel
		return self

def loadStaples(fileName='staples.csv'):
	staples = {}
	with open('staples.csv') as f:
		for line in f:
			staples[staple(line).num]=staple(line)
	return staples;

def getStapleFromBelt(stdict, cellnum,stplpos):
	celldict={}
	stplnum= staplePosInCelltoStapleNum(stplpos)
	if num<=26: shift=39
	else: shift=103
	for item in stdict.values():
		if item.index[:7]=='s%04d.%d'%(num*32+shift,stplnum):
      			celldict[item.label]=stdict.pop(item.label)
	return stdict,celldict

def get_cell_0(stdict, num):
  celldict={}
  if num<=26: shift=39
  else: shift=103
  for item in stdict.values():
    if item.index[:7]=='s%04d.0'%(num*32+shift):
      celldict[item.label]=stdict.pop(item.label)
  return stdict,celldict

def get_cell_1(stdict, numcls):
  celldict={}
  if num<=26: shift=39
  else: shift=103
  for item in stdict.values():
    if item.index[:7]=='s%04d.1'%(num*32+shift+17):
      celldict[item.label]=stdict.pop(item.label)
  return stdict,celldict


def get_cell_2(stdict, num):
  celldict={}
  if num<=26: shift=39
  else: shift=103
  for item in stdict.values():
    if item.index[:7]=='s%04d.2'%(num*32+shift):
      celldict[item.label]=stdict.pop(item.label)
  return stdict,celldict

def get_cell_3(stdict, num):
  """returns staple on helix 3 for cell number num"""
  celldict={}
  if num<=26: shift=39
  else: shift=103
  for item in stdict.values():
    if item.index[:7]=='s%04d.3'%(num*32+shift+17):
      celldict[item.label]=stdict.pop(item.label)
  return stdict,celldict


def make_cell_page(stpldict,ncell,fname,title='Origami Belt Cells'):
	beg="""<html>\n<head>\n"""
	beg += """<title>"""+title+"""</title>\n"""
	beg += """</head>\n<body>\n"""
	beg += """<b>"""+title+"""</b>"""
	colors = ['yellow', 'blue']
	stemsColors = {'0-3':0, '0-5': 1, '1-3':1, '1-5':0, '2-3':1, '2-5':0,'3-3':0, '3-5':1}
	#stems are assigned circularily in an A-B fashion. there are two different stems to prevent staple folding on itself.
	#the assignment above is 'helix# - sequence side'
	for cell in range(ncell):
		beg +=	"<p>Cell " + str(cell) + "</p>\n" +\
	      		"""<table style="font-family:monospace">\n"""
		for helixNum in range(4): #iterate through helices
			stapNum = getStapleNum(stpldict,cell,helixNum)
			stap = stpldict[stapNum]
			stem3Color = colors[stemsColors[str(stap.sh)+'-3']]
			stem5Color = colors[stemsColors[str(stap.sh)+'-5']]
			beg += """<tr>\n""" + \
				"""<td>""" + stap.label + "</td>" + \
				"<td>helix# " + str(stap.sh) + " </td>" + \
				"""<td><span style="color:green">""" + stap.se5 + """</span></td>""" + \
				"""<td><span style="color:""" +stem3Color +""" ">""" + stap.stem5 + "</span></td>" + \
				"""<td><span style="color:gray">""" + stap.bulge5 + "</span></td>"+ \
				"""<td><span style="color:black">""" + stap.body + "</span></td>" + \
				"""<td><span style="color:gray">""" + stap.bulge3 + "</span></td>"+ \
				"""<td><span style="color:""" + stem5Color+"""">""" + stap.stem3 + "</span></td>" + \
				"""<td><span style="color:red">""" + stap.se3 + "</span></td>\n"+ \
				"""</tr>\n""";
		beg += "</table>";

#	    try: 
#	    except IndexError: pass
  	beg+="</body>\n</html>"
	with open(fname, "w") as f: f.write(beg)

def addSEtoCell(stdict,cell,SE,stems,bulge3='TT'):
#Adds sticky ends and stems to all 4 staples in a cell in a circular manner so stems work together.

	helicesOrder = [0,2,3,1]
	A = stems[0];
	B = stems[1];
	Ap = complement(A,reverse='True');
	Bp = complement(B, reverse = 'True');

	stdict[getStapleNum(stdict,cell,helicesOrder[0])].addSE(SE,(A,B),bulge3);
	stdict[getStapleNum(stdict,cell,helicesOrder[1])].addSE(SE, (Bp,Ap),bulge3);
	stdict[getStapleNum(stdict,cell,helicesOrder[2])].addSE(SE,(A,B),bulge3);
	stdict[getStapleNum(stdict,cell,helicesOrder[3])].addSE(SE,(Bp,Ap),bulge3);
	

def moveCellsNicks(stdict,cell,shift):
#Adds sticky ends and stems to all 4 staples in a cell in a circular manner so stems work together.
	helicesOrder = [0,2,3,1]

	body0 = stdict[getStapleNum(stdict,cell,helicesOrder[0])].body
	body1 = stdict[getStapleNum(stdict,cell,helicesOrder[1])].body
	body2 = stdict[getStapleNum(stdict,cell,helicesOrder[2])].body
	body3 = stdict[getStapleNum(stdict,cell,helicesOrder[3])].body

	x = shift

	if x!=0:
		stdict[getStapleNum(stdict,cell,helicesOrder[0])].body = body3[-1*x:]+body0[:-1*x];
		stdict[getStapleNum(stdict,cell,helicesOrder[1])].body = body0[-1*x:]+body1[:-1*x];
		stdict[getStapleNum(stdict,cell,helicesOrder[2])].body = body1[-1*x:]+body2[:-1*x];
		stdict[getStapleNum(stdict,cell,helicesOrder[3])].body = body2[-1*x:]+body3[:-1*x];

def addLabeltoCell(stdict,cell,prefLabel='',sufLabel=''):
#Adds prefix and suffix to staple's label	
	helicesOrder = [0,2,3,1]
	
	stdict[getStapleNum(stdict,cell,helicesOrder[0])].addLabel(prefLabel,sufLabel);
	stdict[getStapleNum(stdict,cell,helicesOrder[1])].addLabel(prefLabel,sufLabel);
	stdict[getStapleNum(stdict,cell,helicesOrder[2])].addLabel(prefLabel,sufLabel);
	stdict[getStapleNum(stdict,cell,helicesOrder[3])].addLabel(prefLabel,sufLabel);


def staplePosInCelltoStapleNum(pos):
	stplStartHelix = {'tl':0,'tr':1,'bl':2,'br':3}
	return stplStartHelix[pos]

def getStapleNum(stdict,cell,stplPos): #pos - which of the 4 double helixes (0,1,2,3)
	stplNum = -1
	for stpl in stdict: #find the staple number from its cell number and helix number
		staple = stdict[stpl]		
		if staple.sh == stplPos and staple.cell == cell:
			stplNum = staple.num
	if stplNum ==-1:
		print 'CAN NOT FIND STAPLE at cell=' + str(cell) + ' and helix=' +str(stplPos)
	return stplNum

def getCellsStaples(stdict,cellNum):
	cell = []
	for helixNum in range(4):
		cell.append(stdict[getStapleNum(stdict,cellNum,helixNum)]);

	return cell
	
#	stpl = stdict['EM'+"%3d"%getStapleNum(stdict,cell,stplPos)]
# 	stpl = stdict[stpl]
'''  

  if cell<=26: shift=39
  else: shift=103
  if p==5:
    nicks={0:1,1:3,3:2,2:0}  
  if p==3:  
    nicks={1:0,3:1,2:3,0:2}  
  for item in stdict.values(): #get staple labels (e.g. "EM003")
    if item.index[:7]=='s%04d.%d'%(cell*32+shift+17*(stapnum%2),stapnum):
      labelSE=item.label
    if item.index[:7]=='s%04d.%d'%(cell*32+shift+17*(nicks[stapnum]%2),nicks[stapnum]):
      labelBL=item.label
  if p==5:
    if not stdict[labelBL].se3 in ['', stems[0],stems[1],complement(stems[0],reverse=True),complement(stems[1],reverse=True)]:
      print "Warning: partner strand already has sticky end: "+stdict[labelBL].se3
      return -1
    if stapnum in [0,3]: 
      stemSE=complement(stems[1],reverse=True)
      stemBL=stems[1]
    else:
      stemSE=complement(stems[0],reverse=True)
      stemBL=stems[0]  
    stdict[labelSE].se5=SE+stemSE
    stdict[labelBL].se3=stemBL
  if p==3:
    if not stdict[labelBL].se5 in ['', stems[0],stems[1],complement(stems[0],reverse=True),complement(stems[1],reverse=True)]:
      print "Warning: partner strand already has sticky end: "+stdict[labelBL].se5
      return -1
    if stapnum in [0,3]: 
      stemBL=complement(stems[0],reverse=True)
      stemSE=stems[0]
    else:
      stemBL=complement(stems[1],reverse=True)
      stemSE=stems[1]  
    
    stdict[labelSE].se3=(' TT' if bulge else '')+ ' ' + stemSE+' '+SE
    stdict[labelBL].se5=stemBL + ' '
    
  stdict[labelSE].seq=stdict[labelSE].se5+stdict[labelSE].body+stdict[labelSE].se3

  stdict[labelBL].seq=stdict[labelBL].se5+stdict[labelBL].body+stdict[labelBL].se3
  return 0 '''




def shiftsystem(staplefile, mapfile,sedict,delimiter=";", typ=""):
    oldstaples={}
    newstaples={}
    copystaples={}
    with open(staplefile) as f:
        for line in f:
            oldstaples[line.split(delimiter)[0]]=staple(line, delim=delimiter, typ=typ)
            copystaples[line.split(delimiter)[0]]=staple(line, delim=delimiter, typ=typ)
    with open(mapfile) as f:
        for line in f:   
            p=line.replace("\n","").split('#')[0].split(delimiter)
            print p[0]
            e5=p[1].split("_")
            e3=p[2].split("_")
            newstaples[p[0]]=copystaples[p[0]]
            try:
                n=int(e5[0])
                if n<0: newstaples[p[0]].body=oldstaples[p[0]].body[-n:]
                if n>0: newstaples[p[0]].body=oldstaples[e5[1]].body[-n:]+newstaples[p[0]].body
            except ValueError: pass
            for k in sedict.keys(): e5[2]=e5[2].replace(k,sedict[k])
            newstaples[p[0]].se5=e5[2]
            try:
                n=int(e3[0])
                if n<0: newstaples[p[0]].body=newstaples[p[0]].body[:n]
                if n>0: newstaples[p[0]].body=newstaples[p[0]].body+oldstaples[e3[1]].body[:n]
            except ValueError: pass
            for k in sedict.keys(): e3[2]=e3[2].replace(k,sedict[k])
            newstaples[p[0]].se3=e3[2]
            newstaples[p[0]].seq=newstaples[p[0]].se5+newstaples[p[0]].body+newstaples[p[0]].se3
    return newstaples, oldstaples
        
           
def TMelt(string, mg="12.5e-3", conc="5e-7", comp='', verbose=False):
    """passes a DNA sequence to the external MELTING Java app and returns
    the approximate melting temperature, entropy and enthalpy. Default values for Mg++
    and DNA concentration can be adjusted as well as the location of the melting5.jar
    executable. Note: This function assumes your system's decimal separator to be a dot '.', not a comma ','!"""
    if len(comp) >0: os.system('melting -S%s -C%s -G%s -P%s -Hdnadna -Omelttemp.txt'%(string,comp[::-1],mg,conc))
    else:os.system('melting -S%s -G%s -P%s -Hdnadna -Omelttemp.txt'%(string,mg,conc))
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
    if os.path.exists('melttemp.txt'): os.remove('melttemp.txt')
    return res


def all_permutations(length):
  lib=['G','A','T','C']
  perms=['']
  for i in range(length):
    perms=list(chain.from_iterable([[p+l for l in lib] for p in perms]))
  return list(perms)


def move_nick(move, cell, stapnum, stdict,verbose=False):
  """always assumes shift happens at 5' end of 'stapnum' staple. 
  Shift goes in direction of scaffold strand (i.e., elongates at 5' end of 'stapnum' staple). 
  Sequence is padded from next staple at nick, which is cut accordingly.""" 
  nicks={0:1,1:3,3:2,2:0}
  if cell<=26: shift=39
  else: shift=103
  for item in stdict.values():
    if item.index[:7]=='s%04d.%d'%(cell*32+shift+17*(stapnum%2),stapnum):
      label5=item.label
      body5=item.body
    if item.index[:7]=='s%04d.%d'%(cell*32+shift+17*(nicks[stapnum]%2),nicks[stapnum]):
      label3=item.label
      body3=item.body
  if verbose: 
    print label5, body5
    print label3, body3
  if move<0:
    body3=body3+body5[:move]
    body5=body5[move:]
  else:
    body5=body3[-move:]+body5
    body3=body3[:-move]
  if verbose:
    print label5, body5
    print label3, body3
  stdict[label3].body=body3
  stdict[label3].seq=stdict[label3].se5+stdict[label3].body+stdict[label3].se3
  stdict[label5].body=body5
  stdict[label5].seq=stdict[label5].se5+stdict[label5].body+stdict[label5].se3 
  print label5, stdict[label5].body
  print label3, stdict[label3].body


def make_cell_pageOLD(cdict,ncell,fname,title='Origami Belt Cells'):
	beg="""<html>\n<head>\n"""
	beg += """<title>"""+title+"""</title>\n"""
	beg += """</head>\n<body>\n"""
	beg += """<b>"""+title+"""</b>"""
  	for i in range(ncell):
	    try:
	      cdict,stap0=get_cell_0(cdict,i)
	      stap0=stap0.values()[0]
	      cdict,stap1=get_cell_1(cdict,i)
	      stap1=stap1.values()[0]
	      cdict,stap2=get_cell_2(cdict,i)
	      stap2=stap2.values()[0]
	      cdict,stap3=get_cell_3(cdict,i)
	      stap3=stap3.values()[0]
	      beg=beg+"""
	      <p>Cell %d</p>
	      <table style="font-family:monospace">
	      <tr>
		<td>%s</td>
		<td>helix 0 </td>
		<td><span style="color:red">%s</span></td>
		<td>%s</td>
		<td><span style="color:red">%s</span></td>
		<td><span style="color:magenta">%s</span></td>
		</tr>
	      <tr>
		<td>%s</td><td>helix 1</td>
		<td><span style="color:blue">%s</span></td>
		<td>%s</td>
		<td><span style="color:blue">%s</span></td>
		<td><span style="color:magenta">%s</span></td>
		</tr>
	      <tr>
		<td>%s</td>
		<td>helix 2</td>
		<td><span style="color:green">%s</span></td>
		<td>%s</td>
		<td><span style="color:green">%s</span></td>
		<td><span style="color:magenta">%s</span></td>
		</tr>
	      <tr>
		<td>%s</td>
		<td>helix 3</td>
		<td><span style="color:cyan">%s</span></td>
		<td>%s</td>
		<td><span style="color:cyan">%s</span></td>
		<td><span style="color:magenta">%s</span></td>
		</tr>
		</table>
	      """%(i,stap0.label,stap0.se5,stap0.body,stap0.se3[:-8],stap0.se3[-8:],
		   stap1.label, stap1.se5,stap1.body,stap1.se3[:-8],stap1.se3[-8:],
		   stap2.label, stap2.se5,stap2.body,stap2.se3[:-8], stap2.se3[-8:],
		   stap3.label, stap3.se5,stap3.body,stap3.se3[:-8],stap3.se3[-8:])
	    except IndexError: pass
  	beg+="</body></html>"
	with open(fname, "w") as f: f.write(beg)

def add_SE_OLD(SE,stems,cell,stapnum, stdict, bulge = True, p=3):
  if cell<=26: shift=39
  else: shift=103
  if p==5:
    nicks={0:1,1:3,3:2,2:0}  
  if p==3:  
    nicks={1:0,3:1,2:3,0:2}  
  for item in stdict.values(): #get staple labels (e.g. "EM003")
    if item.index[:7]=='s%04d.%d'%(cell*32+shift+17*(stapnum%2),stapnum):
      labelSE=item.label
    if item.index[:7]=='s%04d.%d'%(cell*32+shift+17*(nicks[stapnum]%2),nicks[stapnum]):
      labelBL=item.label
  if p==5:
    if not stdict[labelBL].se3 in ['', stems[0],stems[1],complement(stems[0],reverse=True),complement(stems[1],reverse=True)]:
      print "Warning: partner strand already has sticky end: "+stdict[labelBL].se3
      return -1
    if stapnum in [0,3]: 
      stemSE=complement(stems[1],reverse=True)
      stemBL=stems[1]
    else:
      stemSE=complement(stems[0],reverse=True)
      stemBL=stems[0]  
    stdict[labelSE].se5=SE+stemSE
    stdict[labelBL].se3=stemBL
  if p==3:
    if not stdict[labelBL].se5 in ['', stems[0],stems[1],complement(stems[0],reverse=True),complement(stems[1],reverse=True)]:
      print "Warning: partner strand already has sticky end: "+stdict[labelBL].se5
      return -1
    if stapnum in [0,3]: 
      stemBL=complement(stems[0],reverse=True)
      stemSE=stems[0]
    else:
      stemBL=complement(stems[1],reverse=True)
      stemSE=stems[1]  
    
    stdict[labelSE].se3=(' TT' if bulge else '')+ ' ' + stemSE+' '+SE
    stdict[labelBL].se5=stemBL + ' '
    
  stdict[labelSE].seq=stdict[labelSE].se5+stdict[labelSE].body+stdict[labelSE].se3

  stdict[labelBL].seq=stdict[labelBL].se5+stdict[labelBL].body+stdict[labelBL].se3
  return 0

