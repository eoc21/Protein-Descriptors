'''
Created on 4 Dec 2009
Program to extra gene bank files in fasta format
@author: eoc21
'''

from Bio import Entrez
from Bio import SeqIO
import urllib, sys, os
Entrez.email = "eoc210@googlemail.com"
#Read in all files in a cluster and get geneId

pdbFile = open(sys.argv[1],'r')
values = str(pdbFile).split("/")
proteinId = str(values[-1]).split(".")[0]
pdbFileResource = "http://www.rcsb.org/pdb/explore/biologyAndChemistry.do?structureId="+proteinId
dataSource = urllib.urlopen(pdbFileResource)
DS = dataSource.readlines()
ncbiLine = "http://www.ncbi.nlm.nih.gov/entrez/query"
geneIdValue = ""
for i in range(len(DS)):
    if ncbiLine in str(DS[i]):
        splitUpLine = str(DS[i]).split("&nbsp")
        geneId = str(splitUpLine[0]).replace("<td>", "").strip()
geneIdValue =  geneId
#Extract gene sequence
handle = Entrez.efetch(db="nucleotide", rettype="fasta", id=geneIdValue)
seq_record = SeqIO.read(handle, "fasta")
handle.close()
print seq_record.id
print seq_record.seq
#print "%s with %i features" % (seq_record.id, len(seq_record.features))
#v = "MDCLCIVTTKKYRYQDEDTPPLEHSPAHLPNQANSPPVIVNTDTLEAPGYELQVNGTEGEMEYEEITLERGNSGLGFSIAGGTDNPHIGDDPSIFITKIIPGGAAAQDGRLRVNDSILFVNEVDVREVTHSAAVEALKEAGSIVRLYVMRRKPPAEKVMEIKLIKGPKGLGFSIAGGVGNQHIPGDNSIYVTKIIEGGAAHKDGRLQIGDKILAVNSVGLEDVMHEDAVAALKNTYDVVYLKVAKPSNAYLSDSYAPPDITTSYSQHLDNEISHSSYLGTDYPTAMTPTSPRRYSPVAKDLLGEEDIPREPRRIVIHRGSTGLGFNIVGGEDGEGIFISFILAGGPADLSGELRKGDQILSVNGVDLRNASHEQAAIALKNAGQTVTIIAQYKPEEYSRFEAKIHDLREQLMNSSLGSGTASLRSNPKRGFYIRALFDYDKTKDCGFLSQALSFRFGDVLHVIDAGDEEWWQARRVHSDSETDDIGFIPSKRRVERREWSRLKAKDWGSSSGSQGREDSVLSYETVTQMEVHYARPIIILGPTKDRANDDLLSEFPDKFGSCVPHTTRPKREYEIDGRDYHFVSSREKMEKDIQAHKFIEAGQYNSHLYGTSVQSVREVAEQGKHCILDVSANAVRRLQAAHLHPIAIFIRPRSLENVLEINKRITEEQARKAFDRATKLEQEFTECFSAIVEGDSFEEIYHKVKRVIEDLSGPYIWVPARERL"
