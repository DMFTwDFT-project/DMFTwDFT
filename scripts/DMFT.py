#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys, subprocess, os
import numpy as np
import shutil
from shutil import copyfile
import VASP3
import Struct
from INPUT import *
import argparse
from argparse import RawTextHelpFormatter
import pychemia
import re


class Initialize():
	"""DMFTwDFT Initialization and Calculation

	This class contains methods to run the initial DFT and wannier90 calculation
	to generate inputs for the DMFT calculation and performs it. 

	For ionic convergence create a directory DFT_relax under the root directory and have
	DFT input files there.
	
	Run with:
	DMFT.py <options>

	<options>:
	-dft <vasp,siesta>
	-dftexec : Name of the DFT executable. Default: vasp_std
	-relax : This flag turns on DFT relaxation if DFT_relax directory exists
	-dmft : This flag performs dmft calculation
	-hf : This flag performs Hartree-Fock calcualtion
	-force : This flag forces a dmft or hf calculation even if it has been completed
	-kmeshtol : k-mesh tolerance for wannier90


	if the -relax flag is used then the program will check if convergence is reached and rerun
	if necessary. Remember to put an updated version of the DFT inputs for higher convergence 
	inside DFT_relax. This only works with VASP.

	"""

	def __init__(self,args):
		"""
		Contains common functions for all methods.
		This launches the dmft calculation as well.
		"""
		if os.path.exists("para_com.dat"):
			fipa=open('para_com.dat','r')
			self.para_com=str(fipa.readline())[:-1]
			fipa.close()
		else:
			self.para_com=""

		#import the VASP class. This can be used for other DFT codes as well.
		self.DFT = VASP3.VASP_class()

		#dft running directory (current directory)
		self.dir = os.getcwd()

		#name of structure. Required for siesta-> structurename.fdf etc.
		self.structurename = args.structurename

		#kmesh tolerence for wannier mesh
		self.kmeshtol = args.kmeshtol

		#force dmft calculation True of False
		self.force = args.force	

		###################### VASP  ###################################################	
		if args.dft == 'vasp':	
			self.dft = 'vasp'	

			#vasp executable
			self.vasp_exec = 'vasp_std' 
								
			self.gen_win()	
			self.gen_sig()
			if args.relax:
				self.vasp_convergence()

		###################### Siesta  ######################################################		
		if args.dft == 'siesta':		
			self.dft = 'siesta'

			#siesta executable
			self.siesta_exec = 'siesta' 

			self.fdf_to_poscar()
			self.gen_win()	
			self.gen_sig()
		###################################################################################


		#DMFT run					
		if args.dmft:
			self.type = 'DMFT'

		if args.hf:
			self.type = 'HF'

		self.run_dmft()





	def fdf_to_poscar(self):
		"""
		This function converts the siesta .fdf format to POSCAR for further calculations.
		"""
		file = pychemia.code.siesta.SiestaInput(self.structurename+'.fdf') 
		self.st = file.get_structure()
		pychemia.code.vasp.write_poscar(self.st, filepath='POSCAR', newformat=True, direct=True, comment=None) 

	def gen_win(self):
		"""
		This method generates wannier90.win for initial DFT run.
		"""
		
		#generating wannier90.win
		TB=Struct.TBstructure('POSCAR',p['atomnames'],p['orbs'])
		TB.Compute_cor_idx(p['cor_at'],p['cor_orb'])
		print((TB.TB_orbs))
		if list(pV.keys()).count('NBANDS='):
			self.DFT.NBANDS = pV['NBANDS='][0]
		self.DFT.Create_win(TB,p['atomnames'],p['orbs'],p['L_rot'],self.DFT.NBANDS,self.DFT.EFERMI+p['ewin'][0],self.DFT.EFERMI+p['ewin'][1],self.kmeshtol)

		if self.dft =='siesta':

			#Update wannier90.win file then rename it 
			f = open('wannier90.win','a') 			
			f.write('\nbegin unit_cell_cart\n')
			np.savetxt(f,self.st.cell)
			f.write('end unit_cell_cart\n\n')

			#writing the atoms cart block	
			f.write('begin atoms_cart\n')
			aT = (np.array([self.st.symbols]) ).T
			b = self.st.positions
			atoms_cart =np.concatenate((aT,b),axis=1)
			np.savetxt(f,atoms_cart,fmt='%s')
			f.write('end atoms_cart\n\n')

			#writing the mp_grid line
			fi = open(self.structurename+'.fdf')
			data = fi.read()
			fi.close()
			grid = re.findall(r'%block kgrid_Monkhorst_Pack([\s0-9.]*)%endblock kgrid_Monkhorst_Pack',data)[0].split()
			f.write('mp_grid= %s %s %s \n' % (grid[0],grid[5],grid[10]))

			#kpoints
			f.write('\nbegin kpoints\n')
			cmd = 'kmesh.pl '+ grid[0] +' '+grid[5]+' '+grid[10]+' wannier'
			out, err = subprocess.Popen(cmd, shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate() 
			f.write(out.decode('utf-8'))
			if err:
				print(err.decode('utf-8'))
			f.write('end kpoints')
			f.close()
			shutil.copy('wannier90.win',self.structurename+'.win')

	def vasp_run(self,dir):
		"""
		This method runs the inital VASP calculation.
		"""

		#initial VASP run
		print("Running VASP in %s"%dir)
		cmd = 'cd '+dir+ ' && '+ self.para_com+' '+self.vasp_exec #+ " > dft.out 2> dft.error"
		out, err = subprocess.Popen(cmd, shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()
		if err:
			print('DFT calculation failed! Check dft.error for details.\n')
			errdir = dir+os.sep+'dft.error'
			f = open(errdir,'wb')
			f.write(err)
			f.close()
			sys.exit()
		else:
			print('DFT calculation complete.\n')	
			outdir = dir+os.sep+'dft.out'
			f = open(outdir,'wb')
			f.write(out)
			f.close()


	def siesta_run(self,dir):
		"""
		This method runs the initial siesta calculation.
		"""
		#wannier90 pre-processing
		self.run_wan90_pp()	

		#Running siesta
		print('Running Siesta in %s' %dir)	
		cmd = 'cd '+dir+ ' && '+self.siesta_exec+'<'+self.structurename+'.fdf>'+self.structurename+'.out'
		out, err = subprocess.Popen(cmd, shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()


		if os.path.exists(self.structurename+'.out'):

			fi=open(self.structurename+'.out','r')
			done_word=fi.readlines()[-1]
			fi.close()

			if done_word.split()[0] == 'Job':
				print('DFT calculation complete.\n')

			else:
				print('DFT calculation failed!\n')
				sys.exit()

		else:
			print('DFT calculation failed!\n')
			sys.exit()


	def read_outcar(self,outcarpath):
		"""This function reads the OUTCAR file if exists from VASP calculations
		"""
		if os.path.exists(outcarpath+os.sep+'OUTCAR'):
			return pychemia.code.vasp.VaspOutput(outcarpath+os.sep+'OUTCAR')
		else:
			print('OUTCAR not found.')
			return False	

	def vasp_convergence(self):

		"""
		This function checks for convergence inside the DFT_relax directory 
		and copies CONTCAR as POSCAR to root directory. Otherwise it runs vasp for convergence.
		If you want better convergence remember to copy an updated INCAR in the DFT_relax directory.
		"""

		def check_relax(vaspout):
			#Checks for convergence
			ediffg  = abs(pychemia.code.vasp.VaspInput('./DFT_relax/INCAR').EDIFFG) 
			avg_force = vaspout.relaxation_info()['avg_force']
			print('EDIFFG = %f and Average force = %f' % (ediffg,avg_force))
			if avg_force <= ediffg:
				print('Forces are converged.')
				return True
			else:
				print('Forces are not converged.')
				return False

		if os.path.exists('DFT_relax'):
			vaspout = self.read_outcar('DFT_relax')
			if vaspout:
				if vaspout.is_finished and check_relax(vaspout):					
					print('Copying CONTCAR to root directory.\n')
					copyfile('./DFT_relax/CONTCAR','POSCAR')
				else:
					print('Recalculating...')	
					self.vasp_run('./DFT_relax')
					vaspout = self.read_outcar('DFT_relax')							
					if vaspout:
						if vaspout.is_finished and check_relax(vaspout):					
							print('Copying CONTCAR to root directory.\n')
							copyfile('./DFT_relax/CONTCAR','POSCAR')	
						else:
							print('Update convergence parameters. Exiting.')
							sys.exit()	
			else:
				if os.path.isfile('./DFT_relax/INCAR') \
				and os.path.isfile('./DFT_relax/POTCAR') \
				and os.path.isfile('./DFT_relax/POSCAR') \
				and os.path.isfile('./DFT_relax/KPOINTS'):
					print('DFT_relax directory exists. Recalculating...')	
					self.vasp_run('./DFT_relax')
					vaspout = self.read_outcar('DFT_relax')							
					if vaspout:
						if vaspout.is_finished and check_relax(vaspout):					
							print('Copying CONTCAR to root directory.\n')
							copyfile('./DFT_relax/CONTCAR','POSCAR')	
						else:
							print('Update convergence parameters. Exiting.')
							sys.exit()	
				else:
					print('VASP input files missing. Exiting.')
					sys.exit()
		else:
			print('DFT_relax directory does not exist. Convergence failed.')
			sys.exit()

	def update_win(self):
		"""
		This updates the wannier90.win file with the number of bands and fermi energy 
		from the initial DFT calculation.
		"""
		if self.dft == 'vasp':
			#Updating wannier90.win with the number of DFT bands
			self.DFT.Read_NBANDS()
			self.DFT.Read_EFERMI()
			self.DFT.Update_win(self.DFT.NBANDS,self.DFT.EFERMI+p['ewin'][0],self.DFT.EFERMI+p['ewin'][1])

		elif self.dft =='siesta':	
		
			#Reading the Fermi energy and number of bands from siesta output
			fi = open(self.structurename+'.out','r')
			data = fi.read()
			fi.close()

			self.DFT.EFERMI = float(re.findall(r'Fermi\s=[\s0-9+-.]*',data)[0].split()[-1])	
			self.DFT.NBANDS = int(re.findall(r'Siesta2Wannier90.NumberOfBands[\s0-9]*',data)[0].split()[-1])

			self.DFT.Update_win(self.DFT.NBANDS,self.DFT.EFERMI+p['ewin'][0],self.DFT.EFERMI+p['ewin'][1])
			shutil.copy('wannier90.win',self.structurename+'.win')

	def run_wan90_pp(self):
		"""
		This function performs the wannier90 pre-processing required by some DFT codes like siesta. 
		Outputs a .nnkp file which is required for the DFT calculaiton.
		"""
		cmd = 'wannier90.x -pp'+' '+self.structurename
		out, err = subprocess.Popen(cmd, shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()
		if err:
			print(err.decode('utf-8'))
			sys.exit()
		else:
			print(out.decode('utf-8'))	



	def run_wan90(self,filename='wannier90'):
		"""
		Running wannier90.x to generate .chk and .eig files.
		"""

		print('Running wannier90...')
		cmd = "wannier90.x "+filename
		out, err = subprocess.Popen(cmd, shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()
		if err:
			print('wannier90 calculation failed!')
			print(err.decode('utf-8'))
			sys.exit()
		else:
			print('wannier90 calculation complete.')	
			print(out.decode('utf-8'))


	def gen_sig(self):
		"""
		This method generates the initial self energy file sig.inp.
		"""
		cmd = "sigzero.py"
		out, err = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
		if err:
			print(err.decode('utf-8'))
			sys.exit()
		else:
			print('Initial self-energy file generated.\n')


	def copy_files(self):	
		"""
		This creates a directory DMFT or HF in the root directory 
		and copies all the necessary files.
		"""	

		#creating directory for DMFT
		if os.path.exists(self.type):
			if os.path.exists(self.type+"/imp.0/"):
				shutil.rmtree(self.type+"/imp.0/")
				#os.makedirs("DMFT")
		else:	
			os.makedirs(self.type)

		#copying files into DMFT or HF directory
		if self.structurename != None and self.dft != None:
			cmd = "cd "+self.type+" && Copy_input.py ../ "+"-structurename "+self.structurename+" -dft "+self.dft
		else:
			cmd = "cd "+self.type+" && Copy_input.py ../ "
		out, err = subprocess.Popen(cmd, shell=True).communicate()
		if err:
			print('File copy failed!\n')
			print(err)
			sys.exit()
		else:
			print(out)
			print('\n'+self.type+' initialization complete. Ready to run calculation.\n')

	def run_dft(self):
		"""
		This function  calls the dft calculations and the wannier calculations
		"""

		#VASP
		if self.dft =='vasp':
			self.vasp_run(self.dir)
			self.update_win()
			self.run_wan90()
			self.copy_files()

		#Siesta
		elif self.dft =='siesta':
			self.siesta_run(self.dir)	

			#need to rename .eigW to .eig to run wannier90	
			shutil.copy(self.structurename+'.eigW',self.structurename+'.eig')
			self.update_win()
			self.run_wan90(self.structurename)

			#renaming files
			shutil.copy(self.structurename+'.eig','wannier90.eig')
			shutil.copy(self.structurename+'.chk','wannier90.chk')
			shutil.copy(self.structurename+'.win','wannier90.win')
			shutil.copy(self.structurename+'.amn','wannier90.amn')
			self.copy_files()


	def run_dmft(self):
		"""
		This first checks if there is a previous DMFT or HF calculation and runs
		only if that run is incomplete unless forced.
		"""

		#Checking for previous DMFT run in the directory
		pathstr =self.type+os.sep+'INFO_TIME'

		if os.path.exists(pathstr):
			fi=open(pathstr,'r')
			done_word=fi.readlines()[-1]
			fi.close()

			if done_word.split()[0] == 'Calculation':
				print('Existing '+self.type+' calculation is complete.')
				if self.force:
					#forcing DMFT calculation	
					print('-force flag enabled. Restarting '+self.type+'...')
					self.run_dft()
					print('Running '+self.type+' ...')

					if self.dft != None and self.structurename != None:
						if self.type =='HF':
							cmd = 'cd '+self.type+ ' && '+ 'RUNDMFT.py -hf -dft '+self.dft+' -structurename '+self.structurename
						else:
							cmd = 'cd '+self.type+ ' && '+ 'RUNDMFT.py -dft '+self.dft+' -structurename '+self.structurename	

					else:
						if self.type =='HF':
							cmd = 'cd '+self.type+ ' && '+ 'RUNDMFT.py -hf' 
						else:
							cmd = 'cd '+self.type+ ' && '+ 'RUNDMFT.py '

					out, err = subprocess.Popen(cmd, shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()
					if err:
						print(self.type+' calculation failed! Check '+self.type+'.error for details.\n')
						errdir = self.type+os.sep+self.type+'.error'
						f = open(errdir,'wb')
						f.write(err)
						f.close()
						sys.exit()
					else:
						print(self.type+' calculation complete.\n')	
						outdir =self.type+os.sep+self.type+'.out'
						f = open(outdir,'wb')
						f.write(out)
						f.close()

				else:
					#exit when exiting DMFT calculation is complete.
					print('-force flag disabled. Exiting. ')
					sys.exit()


			else:
				#Incomplete DMFT calculation.
				print(self.type+' calculation incomplete.')
				self.run_dft()
				print('Running '+self.type+'...')
				if self.dft != None and self.structurename != None:
					if self.type =='HF':
						cmd = 'cd '+self.type+ ' && '+ 'RUNDMFT.py -hf -dft '+self.dft+' -structurename '+self.structurename
					else:
						cmd = 'cd '+self.type+ ' && '+ 'RUNDMFT.py -dft '+self.dft+' -structurename '+self.structurename	

				else:
					if self.type =='HF':
						cmd = 'cd '+self.type+ ' && '+ 'RUNDMFT.py -hf' 
					else:
						cmd = 'cd '+self.type+ ' && '+ 'RUNDMFT.py '
				out, err = subprocess.Popen(cmd, shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()
				if err:
					print(self.type+' calculation failed! Check '+self.type+'.error for details.\n')
					errdir = self.type+os.sep+self.type+'.error'
					f = open(errdir,'wb')
					f.write(err)
					f.close()
					sys.exit()
				else:
					print(self.type+' calculation complete.\n')	
					outdir =self.type+os.sep+self.type+'.out'
					f = open(outdir,'wb')
					f.write(out)
					f.close()

		else:
			#no DMFT/INFO_TIME found
			self.run_dft()
			print('Running '+self.type+'...')
			if self.dft != None and self.structurename != None:
				if self.type =='HF':
					cmd = 'cd '+self.type+ ' && '+ 'RUNDMFT.py -hf -dft '+self.dft+' -structurename '+self.structurename
				else:
					cmd = 'cd '+self.type+ ' && '+ 'RUNDMFT.py -dft '+self.dft+' -structurename '+self.structurename	

			else:
				if self.type =='HF':
					cmd = 'cd '+self.type+ ' && '+ 'RUNDMFT.py -hf' 
				else:
					cmd = 'cd '+self.type+ ' && '+ 'RUNDMFT.py ' 
			out, err = subprocess.Popen(cmd, shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()
			if err:
				print(self.type+' calculation failed! Check '+self.type+'.error for details.\n')
				errdir = self.type+os.sep+self.type+'.error'
				f = open(errdir,'wb')
				f.write(err)
				f.close()
				sys.exit()
			else:
				print(self.type+' calculation complete.\n')	
				outdir =self.type+os.sep+self.type+'.out'
				f = open(outdir,'wb')
				f.write(out)
				f.close()



if __name__ == "__main__":

	#top level parser
	print('\n---------------------- \n| Welcome to DMFTwDFT |\n----------------------\n')
	des = 'This script performs DFT+DMFT calculations through maximally localized Wannier functions.\n For post-processing, run postDMFT.py.' 
	parser = argparse.ArgumentParser(description=des,formatter_class=RawTextHelpFormatter)
	

	#parser for dft  
	parser.add_argument('-dft',default='vasp', type=str, help='Choice of DFT code for the DMFT calculation.', choices=['vasp','siesta'])
	parser.add_argument('-relax',action='store_true', help='Flag to check for DFT convergence. Program exits if not converged.')
	type_parser=parser.add_mutually_exclusive_group()
	type_parser.add_argument('-dmft',action='store_true',help='Flag to run DMFT. Checks for a previous DMFT calculation and runs only if it is incomplete.')
	type_parser.add_argument('-hf',action='store_true',help='Flag to perform Hartree-Fock calculation to the correlated orbitals.')
	parser.add_argument('-force',action='store_true',help='Flag to force DMFT or HF calculation even if a previous calculation has been completed.')
	parser.add_argument('-structurename', type=str, help='Name of the structure. Not required for VASP. ' )
	parser.add_argument('-kmeshtol',default=0.00001, type=float, help='The tolerance to control if two k-points belong to the same shell in wannier90.')
	args = parser.parse_args()
	Initialize(args)
