#/usr/bin/env python

import re
import os
import sys
import argparse

from scipy.stats import chisquare

def FiltBySite(args):
	File = args.infile[0]
	Format = args.type[0]
	OutFile = args.outfile
	Outtype = args.outtype
	posLeft = args.left
	posLeftout = args.leftout
	Fhead = {'name':'', 'popt':'', 'nloc':0, 'nind':0}
	locMartrix = {}
	locMartrix_new = {}
	if Format == 'onemap':
		headtmp, locMartrix = ParseOM(File)
		Fhead = {'name':'Unknown', 'popt':'Unknown', 'nloc':headtmp[1], 'nind':headtmp[0]}
	elif Format == 'joinmap':
		Fhead, locMartrix = ParseJM(File)
	else:
		pass
	if posLeft:
		leftloc = []
		for each in open(posLeft):
			leftloc.append(each.strip())
		for each in locMartrix:
			if each in leftloc:
				locMartrix_new[each] = locMartrix[each]
				leftloc.remove(each)
		for each in leftloc:
			print "loc '%s' is not in input file" %(each)
	elif posLeftout:
		leftloc = []
		for each in open(posLeftout):
			leftloc.append(each.strip())
		for each in locMartrix:
			if each not in leftloc:
				locMartrix_new[each] = locMartrix[each]
	else:
		pass
	Fhead['nloc'] = len(locMartrix_new)
	if Outtype == 'onemap':
		OutOM([Fhead['nind'],Fhead['nloc']], locMartrix_new, OutFile)
	elif Outtype == 'joinmap':
		OutJM(Fhead, locMartrix_new, OutFile)
		
def CatFile(args): #parameters: Files,Format,Outfile,Outtype
	Files = args.infile
	Format = args.type
	Outfile = args.outfile
	Outtype = args.outtype
	Fhead = {'name':'', 'popt':'', 'nloc':0, 'nind':0}
	locMartrix = {}
	for each in range(len(Files)):
		if Format[each] == 'onemap':
			headtmp, loctmp = ParseOM(Files[each])
			if Fhead['nind'] != 0:
				assert Fhead['nind'] == headtmp[0], "Can't merge loc with different population size."
				Fhead['nloc'] += int(headtmp[1])
			else:
				Fhead['nloc'] = int(headtmp[1])
				Fhead['nind'] = headtmp[0]
			locMartrix = dict(locMartrix, **loctmp)
		elif Format[each] == 'joinmap':
			headtmp, loctmp = ParseJM(Files[each])
			if Fhead['nind'] != 0:
				assert Fhead['nind'] == headtmp['nind'], "Can't merge loc with different population size."
				Fhead['nloc'] += int(headtmp['nloc'])
			else:
				Fhead['nloc'] = int(headtmp['nloc'])
				Fhead['nind'] = headtmp['nind']
			if Fhead['name'] != '':
				assert Fhead['popt'] == headtmp['popt'], "Can't merge loc with different population type."
				if Fhead['name'] != headtmp['name']:
					Fhead['name'] = "Merge"
			else:
				Fhead['name'] = headtmp['name']
				Fhead['popt'] = headtmp['popt']
			locMartrix = dict(locMartrix, **loctmp)
		else:
			pass
	
	if Outtype == 'onemap':
		OMhead = [Fhead['nind'], Fhead['nloc']]
		OutOM(OMhead, locMartrix, Outfile)
	elif Outtype == 'joinmap':
		OutJM(Fhead, locMartrix, Outfile)
	else:
		pass
				

def Filter(args):
	Infile = args.infile[0]
	Outfile = args.outfile
	Outtype = args.outtype
	Format = args.type[0]
	Psd = args.Psd
	Misrate = args.Misrate
	Fhead = {'name':'', 'popt':'', 'nloc':0, 'nind':0}
	locMartrix = {}
	misLoc = []
	sdloc = []
	locMartrix_new = {}
	oldlocnum = 0
	if Format == 'onemap':
		headtmp, locMartrix = ParseOM(Infile)
		Fhead = {'name':'Unknown', 'popt':'Unknown', 'nloc':headtmp[1], 'nind':headtmp[0]}
	elif Format == 'joinmap':
		Fhead, locMartrix = ParseJM(Infile)
	else:
		pass
	oldlocnum = int(Fhead['nloc'])
	for each in locMartrix:
		flag = 1
		locp = CalSD(locMartrix[each][0],locMartrix[each][1:])
		locm = CalMis(locMartrix[each][1:])
		if locp < Psd:
			sdloc.append(each)
			flag = 0
		if locm > Misrate:
			misLoc.append(each)
			flag = 0
		if flag:
			locMartrix_new[each] = locMartrix[each]
	Fhead['nloc'] = len(locMartrix_new)
	if Outtype == 'onemap':
		OutOM([Fhead['nind'],Fhead['nloc']], locMartrix_new, Outfile)
	elif Outtype == 'joinmap':
		OutJM(Fhead, locMartrix_new, Outfile)
	print "Segregation Distortion loc: %s\nMissing loc: %s\nBoth: %s\n" \
		%(len(sdloc), len(misLoc), int(Fhead['nloc'])+len(sdloc)+len(misLoc)-oldlocnum)

def ParseVCF(VCFfile):
	#write later
	pass
	
def OutOM(OMhead, locMartrix, Outfile):
	fout = open(Outfile, 'w')
	fout.write('%s %s\n' %(OMhead[0], OMhead[1]))
	for each in locMartrix:
		fout.write('*%s %s %s\n' %(each, locMartrix[each][0], ','.join(locMartrix[each][1:])))
	fout.close()
	

#parse JoinMap loc file format
def ParseJM(JMfile):
	locMartrix = {}
	JMhead = {}
	segtype = {'abxcd':'A.1', 'efxeg':'A.2', 'hkxhk':'B3.7',
	           'nnxnp':'D2.15', 'lmxll':'D1.10'}
	Mtrans = {'ac':'ac','ad':'ad','bc':'bc','bd':'bd',      #abxcd
			  'ee':'a' ,'eg':'ac','ef':'ba','fg':'bc',      #efxeg
			  'hh':'a' ,'hk':'ab','kk':'b',                 #hkxhk
			  'nn':'a' ,'np':'ab',                          #nnxnp
			  'lm':'ab','ll':'a',                           #lmxll
			  '--':'-'}
	flag = 0
	locn = ''
	for line in open(JMfile):
		line = re.sub(';.*','',line.strip())
		if re.search('^\s*$',line):
			continue
		elif re.search('\=',line):
			line = re.sub('^\s*','',line)
			line = re.sub('\s*$','',line)
			info = re.split('\s*\=\s*',line)
			JMhead[info[0]] = info[1]
		else:
			if JMhead['popt'] == 'CP':
				if not flag:
					locn = re.search('^(\S+)',line).group(1)
					segt = re.search('<(\S+)>', line).group(1)
					locMartrix[locn] = [segtype[segt]]
					flag = 1
				else:
					info = re.split('\s+',line)
					info = [Mtrans[i] for i in info]
					locMartrix[locn] += info
					flag = 0
			#other population type can add here
	nloc = len(locMartrix)
	if nloc != int(JMhead['nloc']):
		print "WARNING: There's %s loc in your file: %s, not equal to %s which in file head info\n" %(nloc, JMfile, JMhead['nloc'])
	return(JMhead, locMartrix)

#parse OneMap loc file format
def ParseOM(OMfile):
	locMartrix = {}
	OMhead = []
	for line in open(OMfile):
		if re.search('^\*',line):
			info = re.search('^\*(\S+)\s+(\S+)\s+(\S+)',line).groups()
			genoinfo = re.split('\s*,\s*',info[2])
			locMartrix[info[0]] = [info[1]] + genoinfo
		elif re.search('^\d',line):
			OMhead = re.split('\s+',line)
		else:
			continue
	nloc = len(locMartrix)
	if nloc != int(OMhead[1]):
		print "WARNING: There's %s loc in your file: %s, not equal to %s which in file head info\n" %(nloc, OMfile, OMhead[1])
	return(OMhead, locMartrix)
	
def OutJM(JMhead, locMartrix, Outfile):
	fout = open(Outfile,'w')
	segtype = {'A.1':'abxcd', 'A.2':'efxeg', 'B3.7':'hkxhk',
	           'D2.15':'nnxnp', 'D1.10':'lmxll'}
	Mtrans = {}
	Mtrans['A.1']   = {'ac':'ac','ad':'ad','bc':'bc','bd':'bd','-':'--'}      #abxcd
	Mtrans['A.2']   = {'a':'ee' ,'ac':'eg','ba':'ef','bc':'fg','-':'--'}      #efxeg
	Mtrans['B3.7']  = {'a':'hh' ,'ab':'hk','b':'kk','-':'--'}                 #hkxhk
	Mtrans['D2.15'] = {'a':'nn' ,'ab':'np','-':'--'}                          #nnxnp
	Mtrans['D1.10'] = {'ab':'lm','a':'ll','-':'--'}                           #lmxll

	fout.write('name = %s\npopt = %s\nnloc = %s\nnind = %s\n\n' \
			%(JMhead['name'],JMhead['popt'],JMhead['nloc'],JMhead['nind']))
	for each in locMartrix:
		fout.write(each + '    ' + '<' + segtype[locMartrix[each][0]] + '>' + '\n')
		geno = ' '.join([Mtrans[locMartrix[each][0]][i] for i in locMartrix[each][1:]])
		fout.write('    ' + geno + '\n')
		
#calculating missing rate
def CalMis(genos):
	totalG = len(genos)
	countMis = 0
	for each in genos:
		if each in ('-','--'):
			countMis += 1
	return float(countMis)/totalG

#calculating P value of Chi-square test for segregation distortion	
def CalSD(genotype,genos):
	genoStat = {}
	Obs = []
	Exp = []
	for each in genos:
		if each in genoStat:
			genoStat[each] += 1
		else:
			genoStat[each] = 1
	if genotype in ('abxcd','A.1'):
		Mac = genoStat['ac']
		Mad = genoStat['ad']
		Mbc = genoStat['bc']
		Mbd = genoStat['bd']
		Obs = [Mac, Mad, Mbc, Mbd]
		tol = sum(Obs)
		Exp = [tol/4.0]*4
	elif genotype in ('aaxab','nnxnp','abxaa','lmxll','D2.15','D1.10'):
		Maa = CalSubset(genoStat, ['aa','nn','a','ll'])
		Mab = CalSubset(genoStat,['ab','np','lm'])
		Obs = [Maa, Mab]
		tol = sum(Obs)
		Exp = [tol/2.0]*2
	elif genotype in ('abxab','hkxhk', 'B3.7'):
		Maa = CalSubset(genoStat, ['aa','hh','a'])
		Mab = CalSubset(genoStat, ['ab','hk'])
		Mbb = CalSubset(genoStat, ['bb','kk','b'])
		Obs = [Maa, Mab, Mbb]
		tol = sum(Obs)
		Exp = [tol/4.0, tol/2.0, tol/4.0]
	elif genotype in ('efxeg', 'abxac', 'A.2'):
		Maa = CalSubset(genoStat, ['aa','ee','a'])
		Mac = CalSubset(genoStat, ['ac','eg'])
		Mba = CalSubset(genoStat, ['ef', 'ba', 'ab'])
		Mbc = CalSubset(genoStat, ['bc', 'fg'])
		Obs = [Maa, Mac, Mba, Mbc]
		tol = sum(Obs)
		Exp = [tol/4.0]*4
	else:
		print "Find a unsupported Marker type %s" %(genotype)

	return Chi_square_test(Obs, Exp)
		
def Chi_square_test(Obs, Exp):
	Pval = chisquare(Obs, Exp)[1]
	return Pval
	
def CalSubset(dic1,set1):
	ins = list(set(set1).intersection(set(dic1.keys())))
	if len(ins) > 1:
		print "Marker is coded in a wrong way!"
	elif len(ins)==1:
		return dic1[ins[0]]
	else:
		return 0

def main():
	parent = argparse.ArgumentParser(add_help=False)
	parent.add_argument('--infile', nargs='+', help='')
	parent.add_argument('--type', nargs='+', help='')
	parent.add_argument('--outfile', required=True, help='')
	parent.add_argument('--outtype', required=True, help='')
	
	parser = argparse.ArgumentParser()
	subparsers = parser.add_subparsers(help="Functions")
	
	parser_CatFile = subparsers.add_parser('CatFile', parents=[parent])
	parser_CatFile.set_defaults(func=CatFile)
	
	parser_Filter = subparsers.add_parser('Filter', parents=[parent])
	parser_Filter.add_argument('--Psd', required=True, type=float, help="")
	parser_Filter.add_argument('--Misrate', required=True, type=float, help="")
	parser_Filter.set_defaults(func=Filter)
	
	parser_FiltBySite = subparsers.add_parser('FiltBySite', parents=[parent])
	parser_FiltBySite.add_argument('--left', required=False, help="")
	parser_FiltBySite.add_argument('--leftout', required=False, help="")
	parser_FiltBySite.set_defaults(func=FiltBySite)
	
	args = parser.parse_args()
	args.func(args)
	
	
if __name__ == '__main__':
	main()
