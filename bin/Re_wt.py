#!/usr/bin/env python3
import os, time, subprocess


#   holder=None
def re_wt(trigger,count_me,bands,trig1,dof,xwt,ywt,zwt,strang,begin,kpt):
   if os.path.exists('OUTCAR'):
      lines=open('OUTCAR','r').readlines()[1:]
      for i, line in enumerate(lines):
         if trigger in lines[i]:
            print('hello darkness my old friend')
            begin=i
            for el in range(i+2,i+2+count_me):
               xwt=float(lines[el].split()[3])
               tmp_nome=str('deltEps_%s.deig' %str(dof))
               cmd="awk '{print $1, $2, $3*(%f)}' " %(xwt)+tmp_nome+" > deltEps_%s_primo.deig" %str(dof)
               tmp_nome=str('deltEps_%s.deig' %str(dof))
               print(os.popen(cmd).read())
               print(xwt)
               dof+=1
               ywt=float(lines[el].split()[4])
               tmp_nome=str('deltEps_%s.deig' %str(dof))
               cmd="awk '{print $1, $2, $3*(%f)}' " %(ywt)+tmp_nome+" > deltEps_%s_primo.deig" %str(dof)
               print(os.popen(cmd).read())
               dof+=1
               print(ywt)
               zwt=float(lines[el].split()[5])
               tmp_nome=str('deltEps_%s.deig' %str(dof))
               cmd="awk '{print $1, $2, $3*(%f)}' " %(zwt)+tmp_nome+" > deltEps_%s_primo.deig" %str(dof)
               print(os.popen(cmd).read())
               dof+=1
               print(zwt)
   for jj in range(1,3*count_me):
      strang+=" deltEps_%s_primo.deig" %str(jj)
#   strang+="deltEps_%s_primo.deig" %str(3*count_me)

   cmd="awk '{a[FNR]=$1;b[FNR]=$2;c[FNR]+=$3} END{for (i=1; i<=FNR; i++) print a[i],b[i],c[i]}'>wannier.deig"+strang
   print(os.popen(cmd).read())


if __name__=='__main__':
   trigger=str(raw_input("Which deg. of fredom do you want from OUTCAR?\n (Please mind your whitespace and if mode is imaginary):"))
   count_me = int(raw_input("How many atoms are there:"))
   bands=int(raw_input("How many Bloch bands did you use:"))
   trig1='band No.'
#   bands= 60
   dof = 1
   xwt=0.
   ywt=0.
   zwt=0.
   strang=""
   begin=10**10
   kpt=1