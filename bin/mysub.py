#!/usr/bin/env python
import os,sys

# this is a class which defines default program parameters and allows them to be updated 
# by supplying a list with entries like [a=2] or [-a,2].  
class inputs:
  def __init__(self,input,pre_exe=''):
    self.funcs=['update','get_dict','report','funcs','allv']
    exec(pre_exe)
    self.__input=input.split()
    for i in input.split():
      if not ['int','float','str'].count(i.split("=")[-1]):
        try:exec('self.%s=%s'%tuple(i.split("="))) 
        except: print "bad input %s given in the inputs class"%i; sys.exit()
    self.allv=[s.split("=")[0] for s in input.split()]
  def update(self,list,pre_exe=''):
    exec(pre_exe)
    if type(list)==type(''):
      if not os.path.isfile(list):print "%s is not a file..."%file;sys.exit()
      else: list=open(list).read().replace("\n"," ").split()

    for j,jj in enumerate(list):
      if jj.count('=') or jj[0]=='-':
        if jj.count('='):ii=jj.split('=')
        else:
          if (len(list)-1)>=(j+1):
            if not list[j+1].count('=') and list[j+1][0]!='-':ii=[jj[1:],list[j+1]]
            else: ii=[jj[1:],'']
          else:ii=[jj[1:],''] 
        if dir(self).count(ii[0]):
          try:
            if type(eval('self.'+ii[0]))==type(''):   exec('self.%s="%s"'%tuple(ii))
            elif type(eval('self.'+ii[0]))==type(1):  exec('self.%s=int(%s)'%tuple(ii))
            elif type(eval('self.'+ii[0]))==type(1.): exec('self.%s=float(%s)'%tuple(ii))
          except: print "error executing command line input paramters %s=%s"%tuple(ii); sys.exit()
        elif [s for s in self.__input if s.split("=")[0]==ii[0]]:
          temp=[s.split("=") for s in self.__input if s.split("=")[0]==ii[0]][0]
          try:
            if temp[-1]=='str':     exec('self.%s="%s"'%tuple(ii))
            elif temp[-1]=='int':   exec('self.%s=int(%s)'%tuple(ii))
            elif temp[-1]=='float': exec('self.%s=float(%s)'%tuple(ii))
          except: print "error executing command line input paramters %s=%s"%tuple(ii); sys.exit()
        else: 
          print """input variable "%s" not found \n\noptions are: """%ii[0], 
          print "  ".join([s for s in dir(self) if s[0:2]!="__" and s!='update']);sys.exit()
      else:
        if list[j-1][0]!='-': print 'improper entry: %s not preceded by -<variable>'%jj
  def report(self,inp='inputs.out'):
    if inp:OUT=open(inp,'w')
    else:OUT=sys.stdout
    for i in [s for s in sorted(dir(self)) if s[0]!="_" and not self.funcs.count(s)]:OUT.write("%s=%s "%(i,eval("self.%s"%i)))
    OUT.write('\n');OUT.close()
  def get_dict(self,inp=[]):
    if not inp:
      d={}
      for i in [s for s in dir(self) if s[0]!="_" and not self.funcs.count(s)]:d[i]=eval("self.%s"%i)
    else:
      d={}
      for i in [s for s in dir(self) if s[0]!="_" and not self.funcs.count(s) and inp.count(s)]:d[i]=eval("self.%s"%i)
    return d

# print matrix for viewing...
def pmatrix(input,prec=2,out=''):
  f=" %."+"%s"%prec+"f " ; tmp=''; maxl=len(str(int(max([max(i) for i in input]))))
  if out:OUT=open(out,'w')
  else: OUT=sys.stdout
  for i in input: 
    for j in i : tmp+=(f%j).rjust(4+maxl+prec)
    tmp+="\n"
  OUT.write(tmp)
  if out:OUT.close()

# update dictionary with command line....
def update_dict(dict,comline=sys.argv):
  for i in comline:
    if i.count("=")==1:
      a1=i.split("=")[0]; a2=i.split("=")[1]   # too dangerous .replace("_"," ")
      if dict.has_key(a1):
         if type(dict[a1])==type(1.0): dict[a1]=float(a2)
         if type(dict[a1])==type(1): dict[a1]=int(a2)
         if type(dict[a1])==type("a"): dict[a1]=str(a2)
         if type(dict[a1])==type(True): dict[a1]=a2
      else: 
         print "variable %s not found in dictionary"%a1; tt=dict.keys();tt.sort() 
         for i in tt:print i
         sys.exit()
    elif i.count('=')>1: print "too many = signs in %s"%i; sys.exit()
  return dict

def latex_table(inp,cap='',cstyle=0,rstyle=0):
  table=r"""
\begin{table}[ht]
\centering
\begin{tabular}{|| %s  ||}
\hline\hline
%s
\hline \hline
\end{tabular}
\caption{%s}
\end{table}
"""
  inp=[s.split() for s in inp.strip().split("\n")]
  nf=len(inp[0])
  temp=[s for s in inp if len(s)!=nf]
  if temp:print "table rows must have same length...\n",temp;sys.exit()
  if   cstyle==0:cols=' c '*nf
  elif cstyle==1:cols=' c | '+' c '*(nf-1)
  elif cstyle==2:cols= " | ".join(['c']*nf)
  if rstyle==1:inp[1][0]=r"\hline"+"\n"+inp[1][0]
  if rstyle==2:
    for i in range(1,len(inp)):inp[i][0]=r"\hline"+"\n"+inp[i][0]
  temp="\n".join([" & ".join(s)+r" \\ " for s in inp])
  return (table%(cols,temp,cap)).strip()

def latex_go(inp,name='temp.eps'):
  #from pyx import canvas,epsfile,trafo,bitmap,path,deco,style,color,document
  #can = canvas.canvas()
  #can.text(0,0,inp,[trafo.scale(3)])
  #can.writeEPSfile(name)
  #os.system("gv %s"%(name))
  # ugh! pyx will not work with tables!!!
  temp=r"\documentclass[14pt]{article}\begin{document}%s\end{document}"
  OUT=open("/tmp/tmplatex.tex",'w');OUT.write(temp%inp);OUT.close()
  os.system("latex /tmp/tmplatex.tex > /dev/null")
  os.system("xdvi tmplatex.dvi >& /dev/null")

def lstr(inp,f='',p=0):
  tt=type(inp)
  if tt==type(0.) or tt==type(0) or tt==type(''):return " %s "%inp   #print "\nlstr takes list or tuple not %s"%inp;sys.exit()
  if not f and not p:
    if len([s for s in inp if type(s)==float])==len(inp):p=4
    elif len([s for s in inp if type(s)==int])==len(inp):f='i'
  if p:temp=(" %."+"%sf "%p)*len(inp)
  elif f=='i': temp=" %i "*len(inp)
  else: temp=" %s "*len(inp)
  return temp%tuple(inp)

def gs_orthog(self,input): 
  import numpy as np
  # assumes that input eigenvectors are rows... and returns rows...
  tol=10**-6
  for i,ii in enumerate(input):
    temp=ii.copy()
    for j in input[:i]:temp-=np.dot(np.conjugate(j),ii)*j
    if sqrt(np.dot(temp,temp))>tol:input[i]=temp/(sqrt(np.dot(temp,np.conjugate(temp))))
    else:raise ZeroDivisionError # "gs_orthog failed: vectors are not linearly independent."
  return input

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


#if __name__=='__main__':
#  temp="1 2 3 4\n 5 6 7 8\n 9 10 11 12"
#  temp=latex_table(temp,cstyle=2,rstyle=0)
#  print temp
#  latex_go(temp)
