#!/usr/bin/env python3
import os, time, subprocess


#   holder=None
def oreo(trigger='Degree', trig1='band No.',bands=5,flag='-1',begin=10**10,kpt=1):   
   if os.path.exists('OUTCAR'):
      lines=open('OUTCAR','r').readlines()[1:]
      for i, line in enumerate(lines):
         if trigger in lines[i]:
            print('hello darkness my old friend')
            begin=i
            flag+=1
            kpt=1
            if flag>0: 
               Dells=open('deltEps_%s.deig' %(str(flag)),'w')
         if i>begin:
            if trig1 in lines[i]:
               for j in range(i+1,i+1+bands):
#                 holder=str( lines[j])
                  bandno=int(lines[j].split()[0])
                  celta= float(lines[j].split()[1])
                 
                  Dells.write("%d        %d        %.6f\n" %(bandno,kpt,celta))
               kpt+=1

if __name__=='__main__':
   trigger='Degree'
   trig1='band No.'
   bands=int(raw_input("Enter how many bands?:"))
   flag = -1
   begin=10**10
   kpt=1
   oreo(trigger,trig1,bands,flag,begin,kpt)