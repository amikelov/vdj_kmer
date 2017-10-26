
def inputFileName = args[0]
def Vtable = args[1]
def discont = args[2].toBoolean()
def useFreq = args[3].toBoolean()

assert args.length == 4



//read V genes similarity table

def Vref = [:]
new File(Vtable).splitEachLine(":") {v ->
	Vref.put(v[0],v[1].split(","))
}



def freq
def kmerMap=[:].withDefault{[:].withDefault{0}}
def kmer
int fileLength = 0
def V
def Vmap
def currentKmers = []
new File(inputFileName).splitEachLine("\t") {seq ->
	fileLength++
	if (!(fileLength % 100000)){
		System.err << fileLength.intdiv(100000)
	}
	if(seq[1].length() > 7) {
		seq[1]=seq[1].substring(3,seq[1].length()-3);

		V = seq[2].split(",");
		for (k=3;k<5;k++){
			currentKmers =[]
			for (i=0;i<seq[1].length()-k+1;i++){
				kmer = seq[1].substring(i,i+k);
				if (currentKmers.contains(kmer)) continue
				currentKmers << kmer
				Vmap = kmerMap.get(kmer)
				if (useFreq) {freq = Double.parseDouble(seq[3])
					} else { freq = 1}
				V.each{v->
					Vref.get(v).each {entry ->		
						if (Vmap.containsKey((entry))) {
							kmerMap.get(kmer) << [(entry) : (Vmap.get((entry)) + 1*freq)]
						} else {
							kmerMap.get(kmer)<< [(entry):1*freq]
						}
					}
				}
			}
		}
	}
}

System.err << "\nDone with continious, going on"

if (discont) {
def kmerDiscontMap = [:]
def kmerD
def regex
def nmatch=0
def result = [:].withDefault{0}
def maps =[]
def usedKeys = []
for (km in kmerMap.keySet()){
	for(l=0;l<km.length()-2;l++){
		kmerD = km.substring(0,l+1) + "." + km.substring(l+2,km.length())
		if (!usedKeys.contains(kmerD) && (kmerD.length()>3) ) {
			usedKeys << kmerD
			regex = ~"$kmerD"
			result = [:].withDefault{0}
			maps=[]
			kmerMap.keySet().grep(regex).each{key ->
				maps<<kmerMap.get(key)
			}		 
			maps.collectMany{it.entrySet()}.each{ result[it.key] += it.value }
			kmerDiscontMap<< ["${kmerD}":result]				
		}	
	}
}

kmerMap<<kmerDiscontMap
}

println(fileLength)
BufferedWriter log = new BufferedWriter(new OutputStreamWriter(System.out));
kmerMap.each{k,v ->
	v.each{Vgene,n ->
		log.write("${k}"+"\t"+"${Vgene}" +"\t"+n +"\n")
	
	}
}
log.flush()
