import os, sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord 
import itertools
from Bio.SeqIO.QualityIO import _get_sanger_quality_str
from optparse import OptionParser

use = "Usage: %prog [options]"

parser = OptionParser(usage = use)

parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=False, help="set mode to verbose.")
parser.add_option("-t", "--treshold", dest="threshold",action="store_true", default=False, help="indicate quality threshold")
parser.add_option("-T", "--temp", dest="temp",action="store_true", default=False, help="indicate the location of the temp file")
parser.add_option("--pe", dest="pairend", action="store_true", default=False, help="treat files as pair-end")

options, args = parser.parse_args()

if options.verbose:
	print "Mode is set to verbose"


#print options.pairend

#print len(args)
#print options.temp

if len(args) == 4 and options.pairend and options.temp: 
  if args[1] == args[2]:
    print 'Error: You have provided the two files that have the same name'
    exit()
  print 'Threshold', args[0]
  print 'Filtering file', args[1]
  print 'Filtering file', args[2]
  print 'Temp file',args[3]
elif len(args) == 3 and options.pairend and not options.temp:
  if args[1] == args[2]:
    print 'Error: You have provided the two files that have the same name'
    exit()
  print 'Threshold', args[0]
  print 'Filtering file', args[1]
  print 'Filtering file', args[2]
elif len(args) == 3 and not options.pairend and options.temp: 
  print 'Threshold', args[0]
  print 'Filtering file', args[1]
  print 'Temp file', args[2]
elif len(args) == 2 and not options.pairend and not options.temp:
  print 'Threshold', args[0]
  print 'Filtering file', args[1]


def filter_reads_pe(iter1,iter2, qual_th):


     for (first, second) in itertools.izip(iter1,iter2):
         qual1 = first.letter_annotations["phred_quality"]
         qual2 = second.letter_annotations["phred_quality"]
         qual1_l = len(qual1)          
         qual2_l = len(qual2)          
         seq1_l =  len(first.seq)
         seq2_l =  len(second.seq)         
         num1 = len([qual1 for qual1 in qual1 if qual1 > int(qual_th)]) 
         num2 = len([qual2 for qual2 in qual2 if qual2 > int(qual_th)]) 
         qual_num2 = qual2_l  - num2  
         qual_num1 = qual1_l  - num1          
         # check if qual and seq are of same length 
         if qual1_l == seq1_l and qual2_l == seq2_l:
           if qual_num2 <= 3 and qual_num1 <= 3:
               yield first
               yield second




def filter_reads_se(records, qual_th):

     for record in records:
         qual = record.letter_annotations["phred_quality"]
         qual_l = len(qual)
         seq_l = len(record.seq)
         num = len([qual for qual in qual if qual > int(qual_th)])
         # check if qual and seq are of same length 
         if qual_l == seq_l:
           if qual_l - num <= 3:
             yield record



def isodd(num):
        return num & 1 and True or False



if (len(args) == 4 and options.pairend and options.temp) or (len(args) == 3 and options.pairend and not options.temp):

  print 'Filtering reads for pair-end data'

  file_in1 = args[1]
  file_in2 = args[2]

  file_out1 = file_in1 + ".fq1"
  file_out2 = file_in2 + ".fq1"


  records_1 = SeqIO.parse(open(file_in1,"rU"), "fastq")
  records_2 = SeqIO.parse(open(file_in2,"rU"), "fastq")


  # transform to fq1 

  handle1 = open(file_out1, "w")
  handle2 = open(file_out2, "w")

  num = 0 

  for rec in filter_reads_pe(records_1, records_2, args[0]):
           num = num + 1
           if isodd(num):        
              handle1.write("@%s\t%s\t+%s\t%s\n" % (
              rec.description,
              rec.seq,
              rec.description,  
              _get_sanger_quality_str(rec)))
           else:
              handle2.write("@%s\t%s\t+%s\t%s\n" % (   
              rec.description,
              rec.seq,
              rec.description,  
              _get_sanger_quality_str(rec)))
  handle1.close()
  handle2.close()


  # nmb of reads filtered 
  filtered = num/2


  print "%i records filtered through" % (filtered) 
  #print "%i records filtered %s" % (filtered, file_out1)  
  print "%i records filtered through" % (filtered)


  # concatinate, sort and remove duplicates

  sorted_reads1 = file_out1 + ".uniq.fq1"
  sorted_reads2 = file_out2 + ".uniq.fq1"


  if options.temp:
    #print 'temp files will be stored in', str(args[3])
    cmd1 = 'cat ' + file_out1 + '| sort -T' + str(args[3]) + ' ' + '| uniq > ' + sorted_reads1
    cmd2 = 'cat ' + file_out2 + '| sort -T' + str(args[3]) + ' ' + '| uniq > ' + sorted_reads2
    #print cmd1
    #print cmd2
  else:
    cmd1 = 'cat ' + file_out1 + '| sort | uniq > ' + sorted_reads1
    cmd2 = 'cat ' + file_out2 + '| sort | uniq > ' + sorted_reads2 

  os.system(cmd1)
  os.system(cmd2)


  # transform back to fastq files, clean fq1 files

  in_fq1 = open(sorted_reads1, "Ur")

  fastq_final = file_in1 + ".Quality.fq"

  out_fq = open(fastq_final , "w")
     
  for line in in_fq1:
           rec = line.split('\t') 
           out_fq.write("%s\n%s\n%s\n%s" % (
           rec[0],
           rec[1],
           rec[2],  
           rec[3]))
         
  in_fq1.close
  out_fq.close



  in_fq1 = open(sorted_reads2, "Ur")

  fastq_final = file_in2 + ".Quality.fq"

  out_fq = open(fastq_final , "w")
     
  for line in in_fq1:
           rec = line.split('\t') 
           out_fq.write("%s\n%s\n%s\n%s" % (
           rec[0],
           rec[1],
           rec[2],  
           rec[3]))
         
  in_fq1.close
  out_fq.close

  cmd = 'rm ' + sorted_reads1 + ' ' + sorted_reads2 + ' ' + file_out1 + ' ' + file_out2

  os.system(cmd)

elif (not options.pairend and len(args) == 3 and options.temp) or (len(args) == 2 and not options.pairend and not options.temp):
  
    print 'Filtering reads for single-end data'
 
    file_in1 = args[1]
    file_out1 = file_in1 + ".fq1"
    records_1 = SeqIO.parse(open(file_in1,"rU"), "fastq")
    handle1 = open(file_out1, "w")

    num = 0 

    for rec in filter_reads_se(records_1, args[0]):
              num = num + 1  
              handle1.write("@%s\t%s\t+%s\t%s\n" % (
              rec.description,
              rec.seq,
              rec.description,  
              _get_sanger_quality_str(rec)))
          
    handle1.close()
    filtered = num 
    print "%i records filtered through" % (filtered)  
    #print "%i records filtered %s" % (filtered, file_out1)  
 
    
    # concatinate, sort and remove duplicates
    sorted_reads1 = file_out1 + ".uniq.fq1"

    if options.temp:
       #print 'temp files will be stored in', str(args[2])
       cmd1 = 'cat ' + file_out1 + '| sort -T ' + str(args[2]) + ' ' + '| uniq > ' + sorted_reads1
       #print cmd1 
    else:
       cmd1 = 'cat ' + file_out1 + '| sort | uniq > ' + sorted_reads1

    os.system(cmd1) 

    # transform back to fastq files, clean fq1 files

    in_fq1 = open(sorted_reads1, "Ur")

    fastq_final = file_in1 + ".Quality.fq"

    out_fq = open(fastq_final , "w")
     
    for line in in_fq1:
           rec = line.split('\t') 
           out_fq.write("%s\n%s\n%s\n%s" % (
           rec[0],
           rec[1],
           rec[2],  
           rec[3]))
         
    in_fq1.close
    out_fq.close

    cmd = 'rm ' + sorted_reads1 + ' ' + file_out1

    os.system(cmd)

elif (not options.pairend and len(args) >= 3 and not options.temp) or (not options.pairend and len(args) >= 4 and options.temp):
    print 'Error: You need to provide one file'
    exit()
elif (options.pairend and len(args) > 3 and not options.temp) or (options.pairend and len(args) >= 3 and options.temp):
    print 'Error: You need to provide two files'
    exit()
elif options.pairend and len(args) < 3:
    print 'Error: You need to provide two files'
    exit()
