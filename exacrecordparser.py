minfreq=float(raw_input("Enter minimum allele frequency threshold as decimal e.g. 0.05: "))

def parseexacrecord(inputstring):
    recordlist=[]

    i=0
    parsestart=0
    parseend=0
    gotstart = False
    gotend = False 

    while i < len(inputstring):
        if inputstring[i+1:(i+15)]=='{"allele_count':
            parseend=i
            if gotstart == True:
                gotend = True
            
        if inputstring[i+1:(i+12)]=='], "start":':
            parseend=i+1
            #print"********************************************"
            if gotstart == True:
                gotend = True

                    
        if inputstring[i:(i+14)]=='{"allele_count':
            parsestart=i
            gotstart = True

        if i==len(inputstring)-1:
            parseend=i+1
            gotend=True
            
        if gotstart == True and gotend ==True:
            recordlist.append(inputstring[parsestart:parseend])
            parsestart=0
            parseend=0
            gotstart = False
            gotend = False
            #print i, parsestart, parseend, inputstring[parsestart:parseend], gotstart, gotend
        i=i+1   

    return recordlist
  
def parselistandwrite(listofrecords):

    import ast #need this to convert strings in list to dictionary

    i=0
    while i < len(listofrecords):
        allelerecord={}
        allelerecord = eval(listofrecords[i].rstrip(","))
        
        recPos = allelerecord["pos"]
        recChr = allelerecord["chrom"]
        recVarID = allelerecord["variant_id"]
        recRefBase = allelerecord["ref"]
        recAltBase = allelerecord["alt"]
        recAlleleFreq = allelerecord["allele_freq"]
        recChange = allelerecord["major_consequence"] # might want to change this to "category"
        recCoding = allelerecord["HGVSc"]
        recProtein = allelerecord["HGVSp"]
        recTotalAlleles = allelerecord["allele_num"]
        recAltAlleles = allelerecord["allele_count"]
        dictpop_allelecounts = allelerecord["pop_acs"]
        dictpop_alleletotals = allelerecord["pop_ans"]

        dictpop_allelefreq={}
        dictpop_allelefreq["East Asian"]= float(dictpop_allelecounts["East Asian"])/float(dictpop_alleletotals["East Asian"])
        dictpop_allelefreq["Other"]= float(dictpop_allelecounts["Other"])/float(dictpop_alleletotals["Other"])
        dictpop_allelefreq["African"]= float(dictpop_allelecounts["African"])/float(dictpop_alleletotals["African"])
        dictpop_allelefreq["Latino"]= float(dictpop_allelecounts["Latino"])/float(dictpop_alleletotals["Latino"])
        dictpop_allelefreq["South Asian"]= float(dictpop_allelecounts["South Asian"])/float(dictpop_alleletotals["South Asian"])
        dictpop_allelefreq["European (Finnish)"]= float(dictpop_allelecounts["European (Finnish)"])/float(dictpop_alleletotals["European (Finnish)"])
        dictpop_allelefreq["European (Non-Finnish)"]= float(dictpop_allelecounts["European (Non-Finnish)"])/float(dictpop_alleletotals["European (Non-Finnish)"])

        thresholdtest = any( v > minfreq for v in dictpop_allelefreq.itervalues())
        
        if thresholdtest == True:
            newline = bedGene+ "\t"+ bedAmplicon+"\t"+str(recChr)+ "\t"+ str(recPos)+ "\t"+ recRefBase+ "\t"+ recAltBase+ "\t"+ str(recAlleleFreq)+ "\t" + str(dictpop_allelefreq["European (Non-Finnish)"])+ "\t"+ str(dictpop_allelefreq["European (Finnish)"])+ "\t" +str(dictpop_allelefreq["African"])+ "\t"+ str(dictpop_allelefreq["Latino"])+ "\t"+ str(dictpop_allelefreq["East Asian"])+ "\t"+ str(dictpop_allelefreq["South Asian"])+ "\t"+ str(dictpop_allelefreq["Other"])+"\n"
            target.write (newline)
        i=i+1


#start writing the output file

target = open("output.txt", 'w')
header = "Gene"+ "\t"+ "Amplicon"+"\t"+"chr"+ "\t"+ "Position"+ "\t"+ "Ref"+ "\t"+ "Alt"+ "\t"+ "Overall Freq"+ "\t"+ "European (Non-Finnish"+ "\t"+ "European (Finnish)"+ "\t" +"African"+ "\t"+ "Latino"+ "\t"+ "East Asian"+ "\t"+ "South Asian"+ "\t"+ "Other"+"\n"
target.write (header)

import urllib
import csv
with open('CHP2_full.bed', 'r') as f:  #change this file name to be your input file
    reader = csv.reader(f, delimiter='\t')
    for row in reader:
            bedChr=row[0]
            bedStart=row[1]
            bedEnd=row[2]
            bedGenefield = row[5].split('=')
            bedGene=bedGenefield[1]
            bedAmplicon = bedChr+":"+bedStart+"-"+bedEnd
            target_url="http://exac.hms.harvard.edu/rest/region/"+bedChr+"-"+bedStart+"-"+bedEnd
            print bedGene, bedAmplicon, target_url
            #import urllib
            sample = urllib.urlopen(target_url).read()
            listofrecords = parseexacrecord(sample)
            parselistandwrite(listofrecords)

target.close()

