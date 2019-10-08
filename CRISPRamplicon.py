#!/usr/bin/python

import sys
import getopt
import optparse
from optparse import OptionParser
from collections import defaultdict
import subprocess
import os
import multiprocessing
import itertools

def hamming2(str1, str2):
    diffs = 0
    for ch1, ch2 in zip(str1, str2):
        if ch1 != ch2:
           diffs += 1
    return diffs


def fasta2dic(fastafile):
    with open(fastafile, 'r') as fasta:
	fastadic=defaultdict(str)
	for line in fasta:
	    if line.startswith(">"):
		header=line.strip()[1:].split(' ')[0]
		continue
	    fastadic[header]+=line.strip()
    return fastadic

def makeRC(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a':'t','c':'g','g':'c','t':'a'}
    reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))
    return reverse_complement

def callblastn(db, query):
    outfile=os.path.basename(query)+'_vs_'+os.path.basename(db)
##    subprocess.call(['blastn', '-db', db , '-task', 'blastn-short', '-query', query, '-num_threads','6', '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore btop qlen slen', '-out', outfile, '-max_target_seqs', '100000000', '-evalue', '0.001'])
    return outfile 

def makeblastdb(db):
    subprocess.call(['makeblastdb', '-in', db , '-dbtype', 'nucl', '-parse_seqids'])
    
def parsespacers(blastfile):
    with open(blastfile, "r") as br:
        query=''
        subject2startend={}
        qsubdic={}
        for line in br:
            linearray=line.split('\t')
            if linearray[0] != query:
                subject2startend={}
                query=linearray[0]
            subject=linearray[1]
            sstart,send=int(linearray[8]),int(linearray[9])
            subject2startend.setdefault(subject,[]).append([sstart,send])
            qsubdic[query]=subject2startend
        return qsubdic
                
def sortdiccoords(qsubdic,db,m):
    batchlist=[]
    m=int(m)
    for repeats in qsubdic:
        fname=repeats+'.'+os.path.basename(db)+'.entry'
        batchlist.append(fname)
        f=open(fname, "w")
        for reads in qsubdic[repeats]:
            qsubdic[repeats][reads].sort(key=lambda x:x[0])
            if len(qsubdic[repeats][reads])>1:
                #print repeats, reads, qsubdic[repeats][reads]
                locs=qsubdic[repeats][reads]
                x=0
                while len(locs)-1>x:
                    s1,s2,e1,e2=locs[x][0],locs[x+1][0],locs[x][1],locs[x+1][1]
                    if s1<e1 and s2-e1<m+2 and s2-e1>m-2:
                        line=' '.join([reads, str(e1+1)+'-'+str(s2-1),'plus','\n'])
                        f.write(line)
                    if s1>e1 and e2-s1<m+2 and e2-s1>m-2:
                        line=' '.join([reads, str(s1+1)+'-'+str(e2-1),'minus','\n'])
                        f.write(line)
                    x=x+1
        f.close()
        if os.stat(fname).st_size==0:
            os.remove(fname)
            batchlist.remove(fname)
    return batchlist

def blastdbcmd(entrybatch, db):
    spacerfasta=entrybatch+'.fasta'
##    subprocess.call(['blastdbcmd', '-entry_batch', entrybatch, '-db', db, '-out', spacerfasta])
    return spacerfasta
    
def fasta2readspacer(fastadict):
    fasta2read={}
    for header in fastadict:
	fasta2read.setdefault(fastadict[header],[]).append(header)
    return fasta2read

def comparefasta(fasta2r):
    comparisonlist=list(itertools.combinations(fasta2r,2))
    print str(len(comparisonlist))+' pairwise comparisons'
    for keys in comparisonlist:
        if hamming2(keys[0], keys[1]) <=3 or hamming2(keys[1],keys[0])<=3 or hamming2(keys[1][1:],keys[0])<=3 or hamming2(keys[1],keys[0][1:])<=3:  
          if keys[1] in fasta2r and keys[0] in fasta2r:
            if len(fasta2r[keys[0]])> len(fasta2r[keys[1]]):
		fasta2r[keys[0]]=fasta2r[keys[0]]+fasta2r[keys[1]]
		fasta2r.pop(keys[1], None)
	    else:
		fasta2r[keys[1]]=fasta2r[keys[1]]+fasta2r[keys[0]]
		fasta2r.pop(keys[0], None)
    return fasta2r

def compareloci(loci2amplicon):
    comparisonlist=list(itertools.combinations(loci2amplicon,2))
    print str(len(comparisonlist))+' pairwise comparisons'
    for keys in comparisonlist:
        if hamming2(keys[0], keys[1]) <1 and keys[1] in loci2amplicon and keys[0] in loci2amplicon:
            #if len(loci2amplicon[keys[0]])> len(loci2amplicon[keys[1]]):
            if len(keys[0].split(',')) > len(keys[1].split(',')):
		loci2amplicon[keys[0]]=loci2amplicon[keys[0]]+loci2amplicon[keys[1]]
		loci2amplicon.pop(keys[1], None)
	    else:
		loci2amplicon[keys[1]]=loci2amplicon[keys[1]]+loci2amplicon[keys[0]]
		loci2amplicon.pop(keys[0], None)
    return loci2amplicon

def splitdictionary(fasta2read,size):
    keys=fasta2read.keys() 
    splitkeys=[keys[i:i+int(size)] for i in xrange(0,len(keys),int(size))] 
    dictionarylist=[]
    for keyslist in splitkeys:
        d={key:fasta2read[key] for key in keyslist}
	dictionarylist.append(d)
    return dictionarylist
	
def combinedictionary(listofd):
    combinedD={}
    for d in listofd:
	for key in d:
	    if key not in combinedD:
	        combinedD[key]=list(set(d[key]))
	    else:
		combinedD[key]=list(set(combinedD[key]+d[key]))
    return combinedD

def combinedictionary1(listofd):
    combineD=listofd[0]
    i=1
    while i < len(listofd):
        #dicpair=(combinedD.keys(), listofd[i].keys())
        combinations=list(itertools.product(combineD.keys(), listofd[i])) 
        print len(combineD), 'combined'
        combineD=combinedictionary([combineD, listofd[i]])
        for keys in combinations:
            if hamming2(keys[0], keys[1]) <=3 or hamming2(keys[1][1:],keys[0])<=3 or hamming2(keys[0][1:],keys[1])<=3: 
              if keys[1] in combineD and keys[0] in combineD:
                if len(combineD[keys[0]])> len(combineD[keys[1]]):
		    combineD[keys[0]]=combineD[keys[0]]+combineD[keys[1]]
		    combineD.pop(keys[1], None)
	        else:
		    combineD[keys[1]]=combineD[keys[1]]+combineD[keys[0]]
		    combineD.pop(keys[0], None)
	#combineD=comparefasta(fasta2r)
        i=i+1
    return combineD
def combinedictionary2(listofd):
    combineD=listofd[0]
    i=1
    while i < len(listofd):
        #dicpair=(combinedD.keys(), listofd[i].keys())
        combinations=list(itertools.product(combineD.keys(), listofd[i])) 
        print len(combineD), 'combined'
        combineD=combinedictionary([combineD, listofd[i]])
        for keys in combinations:
           if hamming2(keys[0], keys[1]) <1:
             if keys[1] in combineD and keys[0] in combineD:
                if len(keys[0].split(','))> len(keys[1].split(',')):
		    combineD[keys[0]]=combineD[keys[0]]+combineD[keys[1]]
		    combineD.pop(keys[1], None)
	        else:
		    combineD[keys[1]]=combineD[keys[1]]+combineD[keys[0]]
		    combineD.pop(keys[0], None)
	#combineD=comparefasta(fasta2r)
        i=i+1
    return combineD

def renamedict(combineddictionary):
    spacercluster=1
    spacerclustername='SP'+str(1)
    newdictionary={}
    newdictionary2={}
    f=open('spacer2cluster.csv', 'w')
    for fasta in combineddictionary:
        for spacers in combineddictionary[fasta]:
            #print spacers+'\t'+spacerclustername
            newdictionary[spacers.split('|')[1]]=spacerclustername
            newdictionary2.setdefault(spacerclustername,[]).append(spacers.split('|')[1])
            f.write(spacers.split('|')[1]+'\t'+spacerclustername+'\n')
        spacercluster+=1
        spacerclustername='SP'+str(spacercluster)
    #write the spacer cluster name file with lengths
    f.close()
    f1=open('cluster2spacer.csv', 'w')
    for cluster in newdictionary2:
	f1.write(cluster+'\t'+','.join(newdictionary2[cluster])+'\t'+str(len(newdictionary2[cluster]))+'\n')
    #write loci to spacer one, within each loci find the most common amplicon
    f1.close()
    return newdictionary, newdictionary2

def entrybatchprocessing(entrybatch,spacer2cluster):
    amplicon2spacer={}
    amplicon2cluster={}
    with open(entrybatch) as coordinates:         
        for line in coordinates:
            #print  line
            amplicon=line.split(' ')[0]
            if 'plus' in line:
	        spacer=amplicon+':'+line.split(' ')[1]
	    elif 'minus' in line:
	        spacer=amplicon+':c'+line.split(' ')[1].split('-')[1]+'-'+line.split(' ')[1].split('-')[0]
            amplicon2spacer.setdefault(amplicon,[]).append(spacer)
            amplicon2cluster.setdefault(amplicon,[]).append(spacer2cluster[spacer])
    return amplicon2spacer,amplicon2cluster
#def amplicon2loci(

def main():
    parser=OptionParser(usage="spacer2blast.py -t taxon -q repeat.fasta -d blastdb -o outfile") 
    parser.add_option("-q", "--query", dest="q",  help="queryfile") 
    parser.add_option("-d", "--dbfile", action="store", dest="db", default=False,help="dbfile") 
    parser.add_option("-o", "--outputfile", action="store", dest='out', default='outfile',help="outputfile") 
    parser.add_option("-m", "--spacerlength", action="store", dest='m', default=32,help="average spacer length") 
    (options,args)=parser.parse_args() 
    if len(args) !=0: 
        parser.error("need input and output") 
    print "This initial blast may take a while.... please wait"
    #makeblastdb(options.db)
    blastreport=callblastn(options.db, options.q)
    blastreport= os.path.basename(options.q)+'_vs_'+os.path.basename(options.db)    
    print "BLAST done"
    qsubdic=parsespacers(blastreport)    
    print len(qsubdic)
    print "parsed spacers"
    batchlist=sortdiccoords(qsubdic,options.db, options.m)
    #print batchlist
    #print "retrieve entry batch"
    for entrybatch in batchlist:
       spacerfasta=blastdbcmd(entrybatch,options.db)
       print "calling blastdbcmd for "+entrybatch
    print batchlist
    numthreads=12
    dic=fasta2dic(spacerfasta)
    print "grabbing spacers"
    print "total spacers:" +str(len(dic))
    fasta2read=fasta2readspacer(dic)
    print "total unique spacers:" + str(len(fasta2read))
    print "splitting spacer dictionary: by 5000"
    dlistsplit=splitdictionary(fasta2read,5000)
    print str(len(dlistsplit)) + ' dictionaries'
    print "using 8 threads"
    print "this may take a few minutes depending on how many threads"
    pool=multiprocessing.Pool(processes=numthreads)
    result=pool.map(comparefasta, (dlistsplit))
    combineD=combinedictionary(result)
    print len(combineD), "spacers condensed"
    dlistsplit=splitdictionary(combineD,10000)
    print str(len(dlistsplit)) + ' dictionaries'
    print "using 8 threads"
    print "this may take a few minutes depending on how many threads"
    pool=multiprocessing.Pool(processes=numthreads)
    result=pool.map(comparefasta, (dlistsplit))
    combineD=combinedictionary(result)
    print len(combineD), "spacers condensed"
    dlistsplit=splitdictionary(combineD,10000)
    print str(len(dlistsplit)) + ' dictionaries'
    print "using 8 threads"
    print "this may take a few minutes depending on how many threads"
    pool=multiprocessing.Pool(processes=numthreads)
    result=pool.map(comparefasta, (dlistsplit))
    combineD=combinedictionary(result)
    print len(combineD), "spacers condensed"
    dlistsplit=splitdictionary(combineD,10000)
    print str(len(dlistsplit)) + ' dictionaries'
    print "using 8 threads"
    print "this may take a few minutes depending on how many threads"
    pool=multiprocessing.Pool(processes=numthreads)
    result=pool.map(comparefasta, (dlistsplit))
    combineD=combinedictionary(result)
    print len(combineD), "spacers condensed"
    dlistsplit=splitdictionary(combineD,10000)
    print str(len(dlistsplit)) + ' dictionaries'
    print "using 8 threads"
    print "this may take a few minutes depending on how many threads"
    pool=multiprocessing.Pool(processes=numthreads)
    result=pool.map(comparefasta, (dlistsplit))
    combineD=combinedictionary(result)
    print len(combineD), "spacers condensed"
    dlistsplit=splitdictionary(combineD,15000)
    print str(len(dlistsplit)) + ' dictionaries'
    print "using 8 threads"
    print "this may take a few minutes depending on how many threads"
    pool=multiprocessing.Pool(processes=numthreads)
    result=pool.map(comparefasta, (dlistsplit))
    combineD=combinedictionary(result)
    print len(combineD), "spacers condensed"
    dlistsplit=splitdictionary(combineD,100)
    result=pool.map(comparefasta, (dlistsplit))
    #combineD=combinedictionary(result)
    #print len(combineD)
    combineD=combinedictionary1(result)
    spacer2cluster,cluster2spacer=renamedict(combineD)
    print "loci assignment"
    amp2spacer, amp2cluster=entrybatchprocessing(batchlist[0],spacer2cluster)
    #print amp2spacer
    print 'total loci:' +str(len(amp2spacer))
    print "loci reduction"
    amp2clusterstring={k:','.join(v) for k,v in amp2cluster.items()}
    #print "total loci: " +str(len(amp2clusterstring))
    #fasta2readspacer(loci #do this part and reduce the same way as spacers)
    clusterstring2amp=fasta2readspacer(amp2clusterstring)
    print "total unique loci: " + str(len(clusterstring2amp)) 
    dlistsplit1=splitdictionary(clusterstring2amp, 2000)
    #print len(dlistsplit1) #remake dlistplit so there cannot be a a dict with just 1
    #for each loci pick longest or most common amplicon
    result1=pool.map(compareloci, (dlistsplit1)) # remake comparefasta to compare loci depends on length of loci
    combineD1=combinedictionary(result1)
    print "total unique loci: " + str(len(clusterstring2amp)) 
    dlistsplit1=splitdictionary(clusterstring2amp, 5000)
    result1=pool.map(compareloci, (dlistsplit1)) # remake comparefasta to compare loci depends on length of loci
    combineD1=combinedictionary(result1)
    print "total unique loci: " + str(len(clusterstring2amp)) 
    dlistsplit1=splitdictionary(clusterstring2amp, 5000)
    result1=pool.map(compareloci, (dlistsplit1)) # remake comparefasta to compare loci depends on length of loci
    combineD1=combinedictionary(result1)
    print "total unique loci: " + str(len(clusterstring2amp)) 
    dlistsplit1=splitdictionary(clusterstring2amp, 5000)
    result1=pool.map(compareloci, (dlistsplit1)) # remake comparefasta to compare loci depends on length of loci
    combineD1=combinedictionary(result1)
    print "total reduced loci: " +str(len(combineD1))
    dlistsplit1=splitdictionary(clusterstring2amp, 500)
    result1=pool.map(compareloci, (dlistsplit1)) # remake comparefasta to compare loci depends on length of loci
    combineD1=combinedictionary2(result1)
    print "total reduced loci: " +str(len(combineD1))
   # print combineD1.keys()
    locicount=1
   # renamed={item.split('_')[1]:combineD1[item] for item in combineD1}
    #for item in combineD1:
    #    print 'L'+str(locicount)+'\t'+item+'\t'+','.join(combineD1[item])+'\t'+str(len(combineD1[item]))
    f=open('locistats.csv', 'w')
    f1=open('longestamplicon.csv', 'w')
    #for item in renamed:
    #    f.write('L'+str(locicount)+'\t'+item+'\t'+','.join(renamed[item])+'\t'+str(len(renamed[item])))
    #    for amplicon in item[renamed]:
    for item in combineD1:
        f.write('L'+str(locicount)+'\t'+item+'\t'+','.join(combineD1[item])+'\t'+str(len(combineD1[item]))+'\n')
        longest=combineD1[item][0]
        for read in combineD1[item]:
            #print amp2spacer[read], len(amp2spacer[longest]), locicount
            if len(amp2spacer[read]) > len(amp2spacer[longest]):
                #print len(amp2spacer[read]), len(amp2spacer[longest])
                longest=read
        f1.write('L'+str(locicount)+'\t'+longest+'\t'+str(len(amp2spacer[longest]))+'\t'+str(len(combineD1[item]))+'\n')
	locicount+=1
    f.close()
    f1.close() 
    #if len(args) !=1: 
    #    parser.error("Need a repeat fasta file") 
    #print options
    #print args


if __name__== '__main__':
    main()
    #callblastn(db,options.out)

