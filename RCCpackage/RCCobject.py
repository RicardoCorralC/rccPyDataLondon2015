#imprime para cada pdb con su chain en el archivo de argumento el size de la proteina y el number of k-Lambdas maximales
#python programa <archivo con pdbnames y chains>  <K>
#PAra ver si la construccion de k-Lambdas viene de residuos contiguos o no
import networkx as nx
import sys, os
from operator import itemgetter
from heapq import nlargest
from collections import defaultdict, Counter
from networkx.algorithms import isomorphism
import RCCdata as rcd
import RCCutils as rcu
import re

class RCCbase(object):
	def __init__(self,pdb,chain,autochain,chain_segments):
		"""
		This constructor creates only a Residue Interaction Graph of a given protein chain.
		chain_segments can be given as a list of pirs indicating start and end position of a segment in chain.
		If chain_segments is an empty list, all chain is considered.
		"""
		self.pdb = pdb
		self.chain = chain

		if autochain:
			chain = rcu.getFirstChain(pdb)
		#os.system('java -cp /home/rcc/libs/ PDBparser ' + pdb + ' ' + chain  + ' > tmpgraph');
		os.system('python make_RIG.py ' + pdb + ' ' + chain  + ' 5.0 > tmpgraph');
		os.system('rm -rf output*')
		fin = open('tmpgraph',"r")

		self.G = nx.Graph()
		self.HG = nx.Graph()
		for linea in fin:
			if linea.strip():
				a, b = map(str,linea.strip().split())
				a_insegment, b_insegment = True, True

				if chain_segments!=[]:
					a_insegment, b_insegment = False, False
					for segment in chain_segments:
						if  (int(a[3:]) in xrange(segment[0],segment[1]+1)): a_insegment = True
						if  (int(b[3:]) in xrange(segment[0],segment[1]+1)): b_insegment = True

				if a_insegment and b_insegment:
					self.G.add_edge(a,b)
		if len(self.G.nodes()) == 0:
			raise Exception("Wrong graph construction")

class RCC(RCCbase):
	def __init__(self,pdb,chain,autochain=False,chain_segments=[]):
		super(RCC,self).__init__(pdb,chain,autochain,chain_segments)
		self.number_of_classes = sum(map(lambda x:len(x),rcd.setSignatures[3:]))
		self.how_many_signatures = Counter()


		self.osisDictString = defaultdict(set) #aqui cada llave signa en setSignatures tendra un set con los osis del tipo signa
		self.osisDict = defaultdict(set) #aqui cada llave signa en setSignatures tendra un set con los osis del tipo signa
		self.osisDictElements = defaultdict(set) #aqui son todos los residuos pertenecientes a un osi tipo signa
		self.regions = defaultdict(set) #regions[res] = [res] o clase de equivalencia de res
		self.factorSet = set() #aqui el conjunto de regiones, o el factor set V(R)/~
		self.ady = defaultdict(set) #los adyacentes a cada residuo; ady[residuo] = {residuos adyacentes a residuo}
		self.adys = defaultdict(type(self.ady)) #adys[signa][residuo] = {residuos adyacentes a residuo tipo signa}
		self.matches = [] #regiones isomorfas
		self.RCCvector = [0]*self.number_of_classes
		self.RCCvector2 = [0]*self.number_of_classes
		self.metainfo_node = defaultdict(tuple) #rc,c; (set(['ASP139', 'THR135', 'LYS111', 'TYR138']), [1, 1, 2])

		self.createR()
		self.getRegions()

	def getSetSignAASeq(self,cliqueAA,gapped=True,gapchar='_'):  
		"""This is a core RCC method"""
		clique = map(lambda x:re.sub("\D", "", x[3:]),cliqueAA)
		AAname = dict()
		for aa in cliqueAA:
			AAname[int(re.sub("\D","",aa[3:]))] = rcd.AA_three_to_one.get(aa[:3],'X')
		r = sorted(clique, key=lambda item: (int(item.partition(' ')[0])
                                  	if item[0].isdigit() else float('inf'), item))
		list_signature = list()
		how_many_consecutive = 1
		secuencia = str()	
		secuencia += AAname[int(r[0])]
		for i in range(1,len(clique)):
			 if(int(r[i])!=int(r[i-1])+1):
				if gapped: secuencia += gapchar
				list_signature.append(how_many_consecutive)
				how_many_consecutive = 1
			 else:
				how_many_consecutive+=1
			 secuencia += AAname[int(r[i])]
		list_signature.append(how_many_consecutive)
		return dict(list_signature=sorted(list_signature),secuencia=secuencia) 


	def createR(self): 
		clases = set() 
		cliques = 0
		for q in  nx.find_cliques(self.G):
			if (len(q) <3) or (len(q)>6) : continue
			cliques += 1
			tmp_list_sign = self.getSetSignAASeq(q)['list_signature']
			self.how_many_signatures[tuple(tmp_list_sign)] += 1	
			L = ','.join(map(lambda(x):str(x),sorted(tmp_list_sign)))
			self.osisDictString[L].add(','.join(q))
			self.osisDict[L].add(tuple(q))
			map(lambda(i):self.osisDictElements[L].add(i),q)

			rcname =  hash(tuple(q))
			self.metainfo_node[rcname] = (set(q),tmp_list_sign)
			self.HG.add_node(rcname)
			for hn in self.HG.nodes():
				if self.metainfo_node[hn][0] & self.metainfo_node[rcname][0]:
					self.HG.add_edge(hn,rcname)

		classindex = 0
		for K in xrange(3,7):
			for signa in rcd.setSignatures[K]:
				self.RCCvector[classindex] = self.how_many_signatures[tuple(signa)]
				for n in self.HG.nodes():
					if self.metainfo_node[n][1] != signa: continue
					self.RCCvector2[classindex] += self.HG.degree(n)
				classindex += 1

	def setUnion(self,a,b):
		return set(a)|set(b)

	def getRegions(self): #necesita que cretateR haya sido invocado antes
		for k in xrange(3,7):
			for signa in rcd.setSignatures[k]:
				llave = ','.join(map(lambda(x):str(x),tuple(signa)))
				ady = defaultdict(set) 
				for res in self.G.nodes():
					ady[res].add(res) #cada residuo es "adyacente" a si mismo, esto para evitar regiones vacias
					for rcc in self.osisDict[llave]:
						if res in rcc:
							map(lambda(r):ady[res].add(r),rcc)
				self.adys[llave] = ady
		for res in self.G.nodes():
			allady = set()
			for k in xrange(3,7):
				for signa in rcd.setSignatures[k]:
					llave = ','.join(map(lambda(x):str(x),tuple(signa)))
					allady.add(tuple(self.adys[llave][res]))
			region = reduce(self.setUnion,allady)
			self.regions[res] = region
			self.factorSet.add(tuple(region))


	def printRCCs(self,k):
		for signa in rcd.setSignatures[k]:
			llave = ','.join(map(lambda(x):str(x),tuple(signa)))
			print '\n>------------------------------------------------\n>Aqui los osis tipo ' + llave
			for osi in self.osisDictString[llave]:
				print osi
			print '\n>Todos juntos:'
			print ','.join(self.osisDictElements[llave])
			print '\n>Todos juntos solo nums: '
			print ','.join(map(lambda(s):s[3:],self.osisDictElements[llave]))

	def printFS(self):
		for clase in self.factorSet:
			print clase

	def printClasesEquiv(self):
		for res in self.G.nodes():
			print '['+res+'] = ' + ','.join(self.regions[res])

	def colores(self):
		colours = ['actinium','aluminum','americium','antimony','argon','arsenic','astatine','barium','berkelium','beryllium','bismuth','bohrium','boron','bromine','cadmium','calcium','californium','carbon','cerium','cesium','chlorine','chromium','cobalt','copper','curium','deuterium','dubnium','dysprosium','einsteinium','erbium','europium','fermium','fluorine','francium','gadolinium','gallium','germanium','gold','hafnium','hassium','helium','holmium','hydrogen','indium','iodine','iridium','iron','krypton','lanthanum','lawrencium','lead','lithium','lutetium','magnesium','manganese','meitnerium','mendelevium','mercury','molybdenum','neodymium','neon','neptunium','nickel','niobium','nitrogen','nobelium','osmium','oxygen','palladium','phosphorus','platinum','plutonium','polonium','potassium','praseodymium','promethium','protactinium','radium','radon','rhenium','rhodium','rubidium','ruthenium','rutherfordium','samarium','scandium','seaborgium','selenium','silicon','silver','sodium','strontium','sulfur','tantalum','technetium','tellurium','terbium','thallium','thorium','thulium','tin','titanium','tungsten','uranium','vanadium','xenon','ytterbium','yttrium','zinc','zirconium']
		i = 0
		while True:
			yield colours[i]
			i = (i+1) % len(colours)



	def printPyMolRegions(self,disjoint=True): 
		"""Shows PyMol commands to visualize regions"""
		color = self.colores()
		print rcd.PyMol_header(self.chain)
		pintados = set()
		for res in self.G.nodes():
			region = map(lambda (x) : x[3:], self.regions[res])
			if disjoint:
				if pintados & set(region): continue
				pintados |= set(region)
			print 'select '+res+' , (resi ' + ','.join(region) + ' , and chain ' + self.chain + ' )'
			#print self.setToSequece(self.regions[res])
			print 'color ' + color.next() + ' , ' + res




	def printPymolIsomorphicRegions(self,disjoint=True):
		color = self.colores()
		print 'hide all'
		print 'set ray_shadows,0'
		print 'set antialias = 1'
		print 'set ray_trace_fog, 2'
		print 'set fog, on'
		print 'set fog, 2'
		print 'set depth_cue, 2'
		print 'set cartoon_fancy_helices, 1'
		print 'bg_color gray'
		print 'color white, all'

		selsA = set()
		selsB= set()
		colsA = set()
		colsB = set() 		
		pintadosA = set()
		pintadosB = set()
		for par in self.matches:
			if disjoint:
				if pintadosA & set(par[0]): continue
				if pintadosB & set(par[1]): continue
				pintadosA |= set(par[0])
				pintadosB |= set(par[1])
			selection_name = str(abs(hash(','.join(par[0])))) #algun nombre se le tiene que poner a la seleccion
			selsA.add('select '+selection_name+' , (resi ' + ','.join(map(lambda(x):x[3:],par[0])) + ' , and chain ' + self.chain + ' )')
			selsB.add('select '+selection_name+' , (resi ' + ','.join(map(lambda(x):x[3:],par[1])) + ' , and chain ' + self.chain2 + ' )')
			colorsin = color.next()
			colsA.add('color ' + colorsin + ' , ' + selection_name)
			colsB.add('color ' + colorsin + ' , ' + selection_name)			

		print 'For protein ' + self.pdb + ':'
		print 'show cartoon, chain ' + self.chain
		for sel in selsA:
			print sel

		print '-----------------------------------------------'

		print 'For protein ' + self.pdb2 + ':'
		print 'show cartoon, chain ' + self.chain2
		for sel in selsB:
			print sel

		print 'For both proteins, color them with this: '
		for color in colsA:
			print color 

	def printSeqRegions(self):
		for par in self.matches:
			print '-----------------------------'
			print self.setToSequece(par[0])
			print self.setToSequece(par[1])
		print ' '

		for par in self.matches:
			print '-----------------------------'
			print self.setToSequece(par[0],withNum=True,threeLetter=True)
			print self.setToSequece(par[1],withNum=True,threeLetter=True)
		print ' '


