#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
from ete3 import Tree
from ete3 import PhyloTree
import numpy as np 
import glob

def numberOfDuplicationGtoS(tree):
    ntrees, ndups, sptrees = tree.get_speciation_trees()
    return ndups

def numberOfLostGtoS(tree):
    compt = 0;
    for node in tree.traverse("preorder"):
        if(node.name == "L"):
            compt += 1
    return compt

def computeStats():
	nb_ensemb_is_better = 0
	nb_esaie_is_better = 0
	nb_ensemb_is_equal_esaie = 0
	au_esaie_values = []
	au_ensembl_values = []
	nb_leaves_values = []

	au_value_treshold = 0.05
	
	result = open("statistics/statistics_80_20.csv", "w")
	result.write("geneTreeID \t  nbLeaves  \t  NbSharedLeaves \t AU_TreeBest  \t AU_Ensembl \t AU_ReconstructEnsembl \t TreeBestPassTest   \t EnsemblTreePassTest \t ReconstructEnsemblTreePassTest \t TreeBestBetter   \t EnsemblBetter \t ReconstructEnsemblBetter \t REC_TreeBest  \t  REC_Ensembl \t REC_ReconstructEnsembl  \t Treebest_REC_is_best     \t Ensembl_REC_is_best \tEnsembl_REC_reconstruct_is_best \n")
	geneFamilyIds = []
	for root, dirs, files in os.walk("ensemblGeneTree/sequences_80_20"): 
			
			for file in files:
				tmp = file.split("_")
				if file.count("_")>0 and len(tmp) > 0 :
					geneFamilyId = tmp[0]
					if not(geneFamilyId in geneFamilyIds) and os.path.exists("ensemblGeneTree/sequences_80_20/" + geneFamilyId + "_stats.txt") == True:
						geneFamilyIds.append(geneFamilyId)

	for geneFamilyId in geneFamilyIds:
		print(geneFamilyId)
		sharedLeavesFile = open("ensemblGeneTree/sequences_80_20/" + geneFamilyId + "_stats.txt", "r")
		lines = sharedLeavesFile.readlines()
		if len(lines)>0:
			line = lines[0]
			line = line.replace("\n", "")
			parts = line.split("\t")
			percentageSharedLeaves = parts[1]
		else:
			percentageSharedLeaves = -10.0
		tree_lines = []
		newickfile = open("ensemblGeneTree/trees_80_20/" + geneFamilyId + ".nw", "r")
		lines = newickfile.readlines()
		for line in lines:
			line = line.replace("\n", "")
			tree_lines.append(line)
		#tree_lines = lines
		nb_input_trees = len(lines)

		nb_leaves_prec = 0
		for line in lines:
			tree_str = line.replace("\n", "")			
			tree = Tree(tree_str)
			nb_leaves = len(tree.get_leaves())
			if nb_leaves_prec == 0:
				nb_leaves_prec = nb_leaves
			else:
				if nb_leaves != nb_leaves_prec:
					print("Input trees dont have the same number of leaves set")
					continue
					#exit("Input trees dont have the same number of leaves set")
					
		if nb_input_trees!=4:
			print("The number of Input tree should be 4")
			exit("The number of Input tree should be and odd number")

		else:

			orthoGroup_result_file_consel_str = "ensemblGeneTree/sequences_80_20/" + geneFamilyId + "_Aln_ByOrthoGroup.phy._seq_OrthoGroup_consel"
			#orthoGroup_result_file_consel_str =  "ensemblGeneTree/sequences_80_20/" + geneFamilyId + "_Aln_ByEnsemblGeneTree.phy._seq_Ensembl_consel"
			ensembl_result_file_consel_str = "ensemblGeneTree/sequences_80_20/" + geneFamilyId + "_Aln_ByEnsemblGeneTree.phy._seq_Ensembl_consel"
			prot_gene_spe = open("/home/local/USHERBROOKE/kuie2201/Bureau/MarieDegen/SpliceGraph/ensemblGeneTree/transcript_protein_gene/" + geneFamilyId + "_update.txt", "r")
			sptree = PhyloTree("(((Drosophila melanogaster:0.155009,Caenorhabditis elegans strain N2:0.188311)Ecdysozoa:0.0001,((Ciona intestinalis:0.151327,Ciona savignyi:0.154463)Ciona:0.0329096,((Eptatretus burgeri:0.128455,Petromyzon marinus:0.169105)Cyclostomata:0.000769775,((((((((((((Astyanax mexicanus Astyanax_mexicanus-2.0:0.00393522,Astyanax mexicanus Astyanax_mexicanus-1.0.2:0.00769478)Astyanax mexicanus:0.0783682,Pygocentrus nattereri:0.0888718)Characoidei:0.0164521,Electrophorus electricus:0.107985)Characiphysae:0.0045476,Ictalurus punctatus:0.100737)Characiphysae:0.00945255,Danio rerio:0.110624)Otophysi:0.0001,(Clupea harengus:0.102439,Denticeps clupeoides:0.119621)Clupeiformes:0.00740506)Otomorpha:0.00259501,(((((((((((Fundulus heteroclitus:0.0802258,Cyprinodon variegatus:0.0839142)Cyprinodontoidei:0.0001,(((Xiphophorus maculatus:0.00748124,Xiphophorus couchianus:0.0110388)Xiphophorus:0.0241812,Gambusia affinis:0.0308088)Poeciliinae:0.00988622,(((Poecilia latipinna:0.00440636,Poecilia formosa:0.00521364)Poecilia:0.00197579,Poecilia mexicana:0.00794921)Poecilia:0.0164401,Poecilia reticulata:0.0285797)Poecilia:0.0101232)Poeciliinae:0.043166)Cyprinodontoidei:0.0137751,Kryptolebias marmoratus:0.0930248)Cyprinodontiformes:0.012096,(((Oryzias latipes ASM223467v1:0.0148451,Oryzias latipes ASM223471v1:0.0165149)Oryzias latipes:0.00101579,Oryzias latipes ASM223469v1:0.0171142)Oryzias latipes:0.0445773,Oryzias melastigma:0.0589407)Oryzias:0.0512964)Atherinomorphae:0.0101824,Gouania willdenowi:0.107449)Ovalentaria:0.0001,(((((Amphiprion percula:0.00551478,Amphiprion ocellaris:0.00676522)Amphiprion:0.0269473,Acanthochromis polyacanthus:0.0400077)Pomacentridae:0.0212392,Stegastes partitus:0.05178)Pomacentridae:0.0186711,Parambassis ranga:0.0783139)Ovalentaria incertae sedis:0.00775579,((((((Maylandia zebra:0.00223424,Astatotilapia calliptera:0.00236576)Haplochromini:0.000578973,Pundamilia nyererei:0.00881603)Haplochromini:0.000671,Haplochromis burtoni:0.00745101)Haplochromini:0.00816067,Neolamprologus brichardi:0.0156223)Pseudocrenilabrinae:0.0105045,Oreochromis niloticus:0.0230369)Pseudocrenilabrinae:0.0381004,Amphilophus citrinellus:0.0554442)Cichlidae:0.0208118)Ovalentaria:0.0185306)Ovalentaria:0.0001,((((Tetraodon nigroviridis:0.067489,Takifugu rubripes:0.074401)Tetraodontidae:0.0352102,Mola mola:0.0977148)Tetraodontiformes:0.0018189,(((Cottoperca gobio:0.0749839,Gasterosteus aculeatus:0.0930361)Perciformes:0.00120023,Larimichthys crocea:0.0742248)Eupercaria:0.00681344,Labrus bergylta:0.0840581)Eupercaria:0.00662427)Eupercaria:0.0001,(((Cynoglossus semilaevis:0.095666,Scophthalmus maximus:0.100244)Pleuronectoidei:0.0001,((Seriola dumerili:0.0127718,Seriola lalandi dorsalis:0.0127882)Seriola:0.0431147,Lates calcarifer:0.0602003)Carangaria:0.0254459)Carangaria:0.0001,((Mastacembelus armatus:0.0747997,Monopterus albus:0.0860603)Synbranchiformes:0.00518692,(Anabas testudineus:0.0721382,Betta splendens:0.0850318)Anabantoidei:0.00851808)Anabantaria:0.00727071)Percomorphaceae:0.00517341)Percomorphaceae:0.00920235)Percomorphaceae:0.0164438,Hippocampus comes:0.11414)Percomorphaceae:0.00159408,Periophthalmus magnuspinnatus:0.114372)Percomorphaceae:0.00367202,Gadus morhua:0.116188)Acanthomorphata:0.00378108,(Hucho hucho:0.0868644,Esox lucius:0.0989656)Protacanthopterygii:0.032928)Euteleosteomorpha:0.0066721)Clupeocephala:0.0001,(Scleropages formosus:0.116167,Paramormyrops kingsleyae:0.119273)Osteoglossiformes:0.0150268)Osteoglossocephalai:0.00643571,Lepisosteus oculatus:0.156606)Neopterygii:0.00646505,Erpetoichthys calabaricus:0.13214)Actinopterygii:0.0001,(((((((((Erinaceus europaeus:0.101026,Sorex araneus:0.107684)Eulipotyphla:0.000657858,((((Pteropus vampyrus:0.0800211,Myotis lucifugus:0.0826389)Chiroptera:0.00170224,(Equus caballus:0.00638598,Equus asinus asinus:0.00665402)Equus:0.0757903)Laurasiatheria:0.00292536,(((((Ursus maritimus:0.00323691,Ursus americanus:0.00365309)Ursus:0.0138412,Ailuropoda melanoleuca:0.0173938)Ursidae:0.0376385,(Mustela putorius furo:0.0160511,Neovison vison:0.0160589)Mustelinae:0.0387644)Caniformia:0.00639456,((Canis lupus dingo:0.000957111,Canis lupus familiaris:0.00135289)Canis lupus:0.0104553,Vulpes vulpes:0.0120047)Canidae:0.0473409)Caniformia:0.00610117,((Panthera pardus:0.00387187,Panthera tigris altaica:0.00451813)Panthera:0.00754061,Felis catus:0.0116044)Felidae:0.0534155)Carnivora:0.0143547)Laurasiatheria:0.00315901,(((((Ovis aries:0.0110756,Capra hircus:0.0112644)Caprinae:0.0170002,((((Bos indicus x Bos taurus UOA_Brahman_1:0.00247744,Bos indicus x Bos taurus UOA_Angus_1:0.00264256)Bos indicus x Bos taurus:0.0001,Bos taurus:0.00186489)Bos:0.00292802,Bos mutus:0.00540361)Bos:0.0001,Bison bison bison:0.00517296)Bovinae:0.0223216)Bovidae:0.0438906,Tursiops truncatus:0.0658466)Cetartiodactyla:0.0100522,Vicugna pacos:0.0818253)Cetartiodactyla:0.00120462,((((((Sus scrofa strain jinhua:0.00184196,Sus scrofa strain meishan:0.00198804)Sus scrofa:0.000258996,Sus scrofa strain rongchang:0.002011)Sus scrofa:0.000304281,Sus scrofa strain tibetan:0.00244605)Sus scrofa:0.0001,Sus scrofa strain bamei:0.0021991)Sus scrofa:0.000107332,Sus scrofa strain wuzhishan:0.00276613)Sus scrofa:0.000516651,((Sus scrofa strain reference:0.00167198,Sus scrofa strain usmarc:0.00196802)Sus scrofa:0.0001,((Sus scrofa strain largewhite:0.00133982,Sus scrofa strain pietrain:0.00150018)Sus scrofa:3.39977e-05,((Sus scrofa strain hampshire:0.00123984,Sus scrofa strain berkshire:0.00154016)Sus scrofa:0.000120658,Sus scrofa strain landrace:0.00145934)Sus scrofa:4.41161e-05)Sus scrofa:0.000169627)Sus scrofa:0.00096732)Sus scrofa:0.0733376)Cetartiodactyla:0.00883307)Laurasiatheria:0.0155472)Laurasiatheria:0.0001,(((((((((Cavia porcellus:0.00818682,Cavia aperea:0.0125132)Cavia:0.0661817,Chinchilla lanigera:0.0705433)Hystricomorpha:0.00414421,Octodon degus:0.0759569)Hystricomorpha:0.00138677,((Heterocephalus glaber HetGla_female_1.0:0.000503489,Heterocephalus glaber HetGla_1.0:0.000766511)Heterocephalus glaber:0.0562146,Fukomys damarensis:0.0577254)Bathyergidae:0.0193335)Hystricomorpha:0.0207161,(((Urocitellus parryii:0.0121268,Spermophilus dauricus:0.0130632)Marmotini:0.000652288,Ictidomys tridecemlineatus:0.0138277)Marmotini:0.00565301,Marmota marmota marmota:0.0183762)Marmotini:0.0711477)Rodentia:0.00260726,(Dipodomys ordii:0.0915058,Castor canadensis:0.0926542)Castorimorpha:0.00511557)Rodentia:0.0044061,(((((((((Mus spicilegus:0.0085474,Mus spretus strain SPRET/EiJ:0.0089626)Mus:0.00109964,(((((((Mus musculus strain reference CL57BL6:0.0001,Mus musculus strain C57BL/6NJ:0.000756133)Mus musculus:0.000884866,Mus musculus strain NZO/HlLtJ:0.00152513)Mus musculus:0.0001,((((Mus musculus strain A/J:0.000329467,Mus musculus strain BALB/cJ:0.000950533)Mus musculus:6.02928e-05,((Mus musculus strain C3H/HeJ:0.000397156,Mus musculus strain CBA/J:0.000522844)Mus musculus:0.000220714,Mus musculus strain DBA/2J:0.00115929)Mus musculus:0.000316136)Mus musculus:0.00033581,Mus musculus strain AKR/J:0.000688249)Mus musculus:0.0001,(Mus musculus strain FVB/NJ:0.000912933,Mus musculus strain NOD/ShiLtJ:0.00148707)Mus musculus:0.000239914)Mus musculus:0.000368911)Mus musculus:0.0001,(Mus musculus strain 129S1/SvImJ:0.000465156,Mus musculus strain LP/J:0.000604844)Mus musculus:0.000970521)Mus musculus:0.000412647,Mus musculus domesticus strain WSB/EiJ:0.00194964)Mus musculus:0.00280959,Mus musculus musculus strain PWK/PhJ:0.00538439)Mus musculus:0.000400067,Mus musculus castaneus strain CAST/EiJ:0.00513373)Mus musculus:0.00549987)Mus:0.0109037,Mus caroli strain CAROLI_EIJ:0.0204399)Mus:0.0168539,Mus pahari strain PAHARI_EIJ:0.0365829)Mus:0.0231643,Rattus norvegicus:0.0629608)Murinae:0.0162398,Meriones unguiculatus:0.0763181)Muridae:0.000636067,(((((Cricetulus griseus CHOK1GS_HDv1:0.000317,Cricetulus griseus CriGri_1.0:0.000693)Cricetulus griseus:0.000777232,Cricetulus griseus CriGri-PICR:0.000667768)Cricetulus griseus:0.0478828,Mesocricetus auratus:0.0525298)Cricetinae:0.0160348,Peromyscus maniculatus bairdii:0.0637651)Cricetidae:0.00356177,Microtus ochrogaster:0.0703313)Cricetidae:0.0104833)Muroidea:0.0120886,Nannospalax galili:0.0881937)Muroidea:0.00981086,Jaculus jaculus:0.0948615)Myomorpha:0.007332)Rodentia:0.00367649,(Oryctolagus cuniculus:0.0778675,Ochotona princeps:0.0908825)Lagomorpha:0.0251812)Glires:0.0001,((((((((((Pan troglodytes:0.002216,Pan paniscus:0.003324)Pan:0.00429938,Homo sapiens:0.00661063)Homininae:0.00183587,Gorilla gorilla:0.00859768)Homininae:0.00840489,Pongo abelii:0.0171815)Hominidae:0.0027716,Nomascus leucogenys:0.0196173)Hominoidea:0.011061,(((Piliocolobus tephrosceles:0.00970118,Colobus angolensis palliatus:0.0111088)Colobinae:0.00169751,(Rhinopithecus roxellana:0.00209587,Rhinopithecus bieti:0.00301413)Rhinopithecus:0.00966999)Colobinae:0.00484466,((((Cercocebus atys:0.00568422,Mandrillus leucophaeus:0.00676578)Cercopithecinae:0.000535902,(Papio anubis:0.00407742,Theropithecus gelada:0.00410258)Cercopithecinae:0.0024316)Cercopithecinae:0.00109526,((Macaca fascicularis:0.00205989,Macaca mulatta:0.00220011)Macaca:0.000877645,Macaca nemestrina:0.00415235)Macaca:0.00437927)Cercopithecinae:0.00407974,Chlorocebus sabaeus:0.0116629)Cercopithecinae:0.00499439)Cercopithecidae:0.0135309)Catarrhini:0.0177941,(((Cebus capucinus imitator:0.0260843,Saimiri boliviensis boliviensis:0.0262457)Cebidae:0.00269488,Callithrix jacchus:0.0300851)Cebidae:0.0001,Aotus nancymaae:0.0252368)Platyrrhini:0.0237389)Simiiformes:0.0275616,Carlito syrichta:0.0797473)Haplorrhini:0.000100626,(((Prolemur simus:0.0371796,Propithecus coquereli:0.0381704)Lemuriformes:0.00513981,Microcebus murinus:0.0414252)Lemuriformes:0.0275709,Otolemur garnettii:0.0747342)Strepsirrhini:0.00703666)Primates:0.0158798,Tupaia belangeri:0.0929595)Euarchontoglires:0.00575616)Euarchontoglires:0.0001)Boreoeutheria:0.00121202,(((Loxodonta africana:0.0623719,Procavia capensis:0.0761881)Afrotheria:0.0167982,Echinops telfairi:0.0966318)Afrotheria:0.019805,(Dasypus novemcinctus:0.0777143,Choloepus hoffmanni:0.0782357)Xenarthra:0.0226405)Eutheria:0.00463665)Eutheria:0.0132897,((((Vombatus ursinus:0.0353722,Phascolarctos cinereus:0.0375378)Diprotodontia:0.0214841,Notamacropus eugenii:0.0630759)Diprotodontia:0.0122477,Sarcophilus harrisii:0.0663409)Marsupialia:0.00209439,Monodelphis domestica:0.0747086)Marsupialia:0.0424224)Theria:0.00745259,Ornithorhynchus anatinus:0.127283)Mammalia:0.00148067,(((((Pogona vitticeps:0.104868,Anolis carolinensis:0.105202)Iguania:0.0134654,Notechis scutatus:0.111395)Toxicofera:0.0001,Salvator merianae:0.115591)Episquamata:0.0122703,Sphenodon punctatus:0.127363)Lepidosauria:0.000969327,((((((((Gallus gallus:0.0369694,Meleagris gallopavo:0.0430606)Phasianidae:0.00645001,Coturnix japonica:0.043595)Phasianidae:0.00126886,Numida meleagris:0.0489561)Galliformes:0.0297436,(Anas platyrhynchos platyrhynchos:0.0284988,Anser brachyrhynchus:0.0314312)Anatidae:0.0385057)Galloanserae:0.00751628,(((Calidris pugnax:0.016933,Calidris pygmaea:0.018207)Calidris:0.0469876,Melopsittacus undulatus:0.0749724)Neognathae:0.0086076,((((Cyanistes caeruleus:0.014851,Parus major:0.015359)Paridae:0.030815,Ficedula albicollis:0.04559)Passeriformes:0.000514931,(((Lonchura striata domestica:0.0200049,Taeniopygia guttata:0.0208851)Estrildinae:0.0209595,Serinus canaria:0.0412755)Passeroidea:0.0001,(Zonotrichia albicollis:0.0116358,Junco hyemalis:0.0129942)Passerellidae:0.0283156)Passeriformes:0.00860953)Passeriformes:0.0173337,(Lepidothrix coronata:0.0119296,Manacus vitellinus:0.0130504)Pipridae:0.0536023)Passeriformes:0.0105191)Neognathae:0.0112947)Neognathae:0.00150368,((((Apteryx owenii:0.00131082,Apteryx haastii:0.00213918)Apteryx:0.00323517,Apteryx rowi:0.00448483)Apteryx:0.0331671,Dromaius novaehollandiae:0.0391011)Palaeognathae:0.0308303,Nothoprocta perdicaria:0.0671553)Palaeognathae:0.0167124)Aves:0.0215438,Crocodylus porosus:0.112377)Archosauria:0.00203969,(((Chelonoidis abingdonii:0.0225391,Gopherus agassizii:0.0228709)Testudinidae:0.0125796,Chrysemys picta bellii:0.0340404)Testudinoidea:0.0413783,Pelodiscus sinensis:0.0782918)Cryptodira:0.0404071)Archelosauria:0.0172966)Sauria:0.00707042)Amniota:0.0125733,Xenopus tropicalis:0.145985)Tetrapoda:0.0001,Latimeria chalumnae:0.11922)Sarcopterygii:0.017797)Euteleostomi:0.0001,Callorhinchus milii:0.14598)Gnathostomata:0.00454606)Vertebrata:0.0418931)Chordata:0.00968241)Bilateria:0.120573,Saccharomyces cerevisiae strain S288C:0.120573)ENSEMBLTREE:1;", format=1)
			for n in sptree.traverse("postorder"):
				if n.is_leaf():
					n.name = n.name.replace(" ", "")
					n.name = n.name.replace("/", "")
					n.name = n.name.replace("_", "")        
				else:
					n.name = n.name.replace(" ", "")
					n.name = n.name.replace("/", "")
					n.name = n.name.replace("_", "")        
			dictSpe = {}
			dictSpeInv = {}
			compt = 10
			for node in sptree.traverse("preorder"):
				if node.is_leaf():

					var = "A" + str(compt)
					compt = compt + 1
					dictSpe[node.name] = var
					dictSpeInv[var] = node.name
					node.name = var

			lines = prot_gene_spe.readlines()
			mapping_gene_species = {}
			for line in lines:
				line = line.replace("\n", "")
				parts = line.split("\t")
				species_name = parts[2]
				species_name = species_name.replace(" ", "")
				species_name = species_name.replace("/", "")
				species_name = species_name.replace("_", "")        
				#species_name = species_name.replace("-", "")        
				#species_name = species_name + "*"
				#species_name = species_name.upper()				
				mapping_gene_species[parts[1]] = species_name


			if os.path.exists(orthoGroup_result_file_consel_str) : #and os.path.exists(ensembl_result_file_consel_str):
				#ensembl_result_file_consel = open(ensembl_result_file_consel_str, "r")
				ensembl_result_file_consel = open(orthoGroup_result_file_consel_str, "r")
				orthoGroup_result_file_consel = open(orthoGroup_result_file_consel_str, "r")
				lines_ensembl_sequences = ensembl_result_file_consel.readlines()
				lines_orthoGroup_sequences = orthoGroup_result_file_consel.readlines()

				if len(lines_orthoGroup_sequences) - 5 == nb_input_trees:
					nb_result_for_each_tree_building_approach = (nb_input_trees - 1)/2

					au_treeBest = 0.0
					au_nj = 0.0
					au_ensembl = 0.0
					au_reconstruct_ensembl = 0.0

					i_best_treeBest = 4
					i_best_nj = 5
					i_ensembl = 6
					i_reconstruct_ensembl = 7

					line = lines_orthoGroup_sequences[i_best_treeBest]
					parts = line.split()
					au_treeBest = float(parts[4])

					line = lines_orthoGroup_sequences[i_best_nj]
					au_nj = line.split()
					au_treeBest = float(parts[4])
					
					line = lines_orthoGroup_sequences[i_ensembl]
					parts = line.split()
					au_ensembl = float(parts[4])
					
					line = lines_orthoGroup_sequences[i_reconstruct_ensembl]
					parts = line.split()
					au_reconstruct_ensembl = float(parts[4])


					au_treeBest_pass_test = 0
					au_reconstruct_ensembl_pass_test = 0
					au_ensembl_pass_test = 0
					au_orthoGroup_pass_test = 0

					if au_treeBest >= au_value_treshold:
						au_treeBest_pass_test = 1

					#if au_nj >= au_value_treshold:
					#	au_nj_pass_test = 1

					if au_ensembl >= au_value_treshold:
						au_ensembl_pass_test = 1

					if au_reconstruct_ensembl >= au_value_treshold:
						au_reconstruct_ensembl_pass_test = 1

					#if au_treeBest_pass_test ==1 or au_nj_pass_test == 1:
					#	au_orthoGroup_pass_test = 1
					
					au_treebest_is_better = 0
					au_reconstruct_ensembl_is_better = 0
					au_ensembl_is_better = 0
					
					if au_treeBest > au_ensembl and au_treeBest>au_reconstruct_ensembl:
						au_treebest_is_better = 1
					elif au_ensembl > au_treeBest and au_ensembl>au_reconstruct_ensembl:
						au_ensembl_is_better = 1
					elif au_reconstruct_ensembl > au_ensembl and au_reconstruct_ensembl>au_treeBest:
						au_reconstruct_ensembl_is_better = 1	
					

					#i_best_ensembl -= 1
					#print(tree_lines, i_best_treeBest, i_best_nj, i_best_ensembl)
					#print(dictSpe.keys())
					treeFromtreeBest = PhyloTree(tree_lines[0])
					#treeFromNJ = PhyloTree(tree_lines[i_best_nj])
					treeFromEnsembl = PhyloTree(tree_lines[2])
					treeFromReconstructEnsembl = PhyloTree(tree_lines[3])

					rf, max_rf, common_leaves, parts_t1, parts_t2, v5, v6 = treeFromtreeBest.robinson_foulds(treeFromEnsembl)

					
					dictTreeBest = {}
					dictReconstructEnsembl ={}
					dictEnsenblBest = {}
					compt = 10
					for node in treeFromtreeBest.traverse("preorder"):
						if node.is_leaf():
							var = "B" + str(compt)
							compt = compt + 1
							dictTreeBest[var] = node.name
							node.name = var
					
					compt = 10
					for node in treeFromReconstructEnsembl.traverse("preorder"):
						if node.is_leaf():
							var = "C" + str(compt)
							compt = compt + 1
							dictReconstructEnsembl[var] = node.name
							node.name = var
						
					compt = 10
					for node in treeFromEnsembl.traverse("preorder"):
						if node.is_leaf():
							var = "D" + str(compt)
							compt = compt + 1
							dictEnsenblBest[var] = node.name
							node.name = var
					
					if rf == 0:
						print("_____________\t", geneFamilyId)
						continue																
					
					for node in treeFromtreeBest.traverse("postorder"):
						if node.is_leaf():					
							node.name =   dictSpe[mapping_gene_species[dictTreeBest[node.name]]] + "_" + node.name
					#exit()
					#print(sptree.write())
					recon_tree1, events1 = treeFromtreeBest.reconcile(sptree)
					nbDup1 = numberOfDuplicationGtoS(recon_tree1)
					nbLost1 = numberOfLostGtoS(recon_tree1)
					reconciliation_score1 = nbDup1 + nbLost1

					for node in treeFromReconstructEnsembl.traverse("postorder"):
						if node.is_leaf():							
							node.name = dictSpe[mapping_gene_species[dictReconstructEnsembl[node.name]]]  + "_" + node.name
					recon_tree2, events2 = treeFromReconstructEnsembl.reconcile(sptree)
					nbDup2 = numberOfDuplicationGtoS(recon_tree2)
					nbLost2 = numberOfLostGtoS(recon_tree2)
					reconciliation_score2 = nbDup2 + nbLost2

					for node in treeFromEnsembl.traverse("postorder"):
						if node.is_leaf():
							#node.name =  dictSpe[mapping_gene_species[node.name]] + "_" + node.name
							node.name =  dictSpe[mapping_gene_species[dictEnsenblBest[node.name]]] + "_" + node.name
					recon_tree3, events3 = treeFromEnsembl.reconcile(sptree)
					nbDup3 = numberOfDuplicationGtoS(recon_tree3)
					nbLost3 = numberOfLostGtoS(recon_tree3)
					reconciliation_score3 = nbDup3 + nbLost3

					reconciliation_score1_is_better = 0
					reconciliation_score2_is_better = 0
					reconciliation_score3_is_better = 0
					

					if reconciliation_score1 < reconciliation_score2 and reconciliation_score1<reconciliation_score3:
						reconciliation_score1_is_better = 1
					elif reconciliation_score2 < reconciliation_score1 and reconciliation_score2<reconciliation_score3:
						reconciliation_score2_is_better = 1
					elif reconciliation_score3 < reconciliation_score1 and reconciliation_score3<reconciliation_score2:
						reconciliation_score3_is_better = 1


					print("writing result...")
					result.write(geneFamilyId + "\t" + str(nb_leaves) +  "\t" + str(percentageSharedLeaves) +  "\t" + str(au_treeBest) + "\t" + str(au_ensembl) + "\t" + str(au_reconstruct_ensembl) + "\t" 
						+	str(au_treeBest_pass_test)  + "\t" + str(au_ensembl_pass_test) + "\t" + str(au_reconstruct_ensembl_pass_test) + "\t" +  str(au_treebest_is_better) 
						+ "\t" + str(au_ensembl_is_better)  + "\t" + str(au_reconstruct_ensembl_is_better)  + "\t"+ 
						str(reconciliation_score1) + "\t" + str(reconciliation_score2) + "\t" +  str(reconciliation_score3) + "\t" + str(reconciliation_score1_is_better)  + "\t" + str(reconciliation_score2_is_better)  + "\t"+ 	str(reconciliation_score3_is_better) + "\n")					
					
				else:
					print(geneFamilyId, "missing result in consel file")
	result.close()	
	"""
	result.write("\n\n\n")
	result.write("Mean \t " +  str(round(np.mean(nb_leaves_values),3)) + "\t" + str(round(np.mean(au_esaie_values),3)) + "\t" + str(round(np.mean(au_ensembl_values),3)) + "\n")
	result.write("Median \t " +  str(round(np.median(nb_leaves_values),3)) + "\t" + str(round(np.median(au_esaie_values),3)) + "\t" + str(round(np.median(au_ensembl_values),3)) + "\n")
	result.write("std \t " +  str(round(np.std(nb_leaves_values),3)) + "\t" + str(round(np.std(au_esaie_values),3)) + "\t" + str(round(np.std(au_ensembl_values),3)) + "\n")
	result.write("nb_ensemb_is_better \t"+ str(nb_ensemb_is_better) + "\n")
	result.write("nb_esaie_is_better \t"+ str(nb_esaie_is_better) + "\n")
	result.write("nb_ensemb_is_equal_esaie \t"+ str(nb_ensemb_is_equal_esaie) + "\n")
	"""

def sharedLeaves():
	percentageSharedLeaves = []
	for x in glob.glob("ensemblGeneTree/sequences_80_20/*_stats.txt"): 

		sharedLeavesFile = open(x, "r")
		lines = sharedLeavesFile.readlines()

		if len(lines)>0:
			line = lines[0]
			line = line.replace("\n", "")
			parts = line.split("\t")
			if len(parts)>1:
				percentageSharedLeaves.append(float(parts[1]))
	print(percentageSharedLeaves)
	print(np.mean(percentageSharedLeaves), np.std(percentageSharedLeaves))


#computeStats()
sharedLeaves()

"""

if __name__ == "__main__": 
	try:
		computeStats()
	except Exception as e:
		print(e)
		pass
"""
