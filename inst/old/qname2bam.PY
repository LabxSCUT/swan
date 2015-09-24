#!/usr/bin/env python
import argparse, pysam, re
def main():

  parser = argparse.ArgumentParser()
  parser.add_argument("bamFile", metavar="bamFile", \
            help="the input bam file")
  parser.add_argument("qnameFile", metavar="qnameFile", \
            help="the input qname file")
  parser.add_argument("outBamFile", metavar="outBamFile", \
            help="the output bam file")
  parser.add_argument("-c", "--chrName", dest="chrName", default="11", \
            help="chromosome name")

  arg_namespace = parser.parse_args()
  bamFile = vars(arg_namespace)['bamFile']
  qnameFile = vars(arg_namespace)['qnameFile']
  outBamFile = vars(arg_namespace)['outBamFile']
  chrName = vars(arg_namespace)['chrName']

  #print set((open(qnameFile, "rU").read()).split("\n"))
  qnameSet=set((open(qnameFile, "rU").read()).split("\n"))
  #print "\n".join(qnameSet)

  samInput = pysam.Samfile( bamFile, "rb" )
  samOutput = pysam.Samfile( outBamFile, "wb", template=samInput )
  for read in samInput.fetch():
    #print read.qname
    if read.qname in qnameSet:
      #print read
      samOutput.write(read)

  samOutput.close()

  samSoft = pysam.Samfile( outBamFile, "rb" )

  def cigar_search(pattern, string):
    try:
      return re.search(pattern,string).group('size')
    except AttributeError:
      return '0'

  #Pysam use 0-based coordinates
  #create a list of clipping position
  al_pattern=re.compile("(?P<size>[0-9]+)[MX=IN]")
  head_pattern=re.compile("^(?P<size>[0-9]+)S")
  tail_pattern=re.compile("(?P<size>[0-9]+)S$")
  allSoft=samSoft.fetch()
  allCigar=[read.cigar for read in allSoft]
  #allCigar=["30S30M40S","51M49S","31S69M","100M"] # "100M" will throw an exception to be catched and assign 0
  #allCigar=["30S30M40S","51M49S","31S69M"]
  #re.search(al_pattern,cigar).group('size')
  allHead=[cigar_search(head_pattern,cigar) for cigar in allCigar]
  allTail=[cigar_search(tail_pattern,cigar) for cigar in allCigar]
  allAlign=[cigar_search(al_pattern,cigar) for cigar in allCigar]
  #calculate the soft position: if head: lpos=pos-head; else if tail: rpos=pos+align
  #pile up soft at each position
  #pile up cvg at each position

  samInput.close()

if __name__=="__main__":
  main()

#chr  pos strand  soft  cvg
