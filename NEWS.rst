News
====

1.6.2
-----
*Release date: 16 July 2014*

* enforce the use of pymongo version==2.8


1.6.1
-----

*Release date: 16 September 2014*

* added taxon tables for viruses
* added gra (genome relative abundance) calculations
* metrics are now calculated for i) reads aligning to pathogen only (excluding those aligning simultaneously to the human ref (i.e digital subtraction)) and ii) reads aligning simultaneously to both pathogen and human refs
* paired reads from the human ref pair-end BAM file are considered valid only if aligned as proper pairs
* changed the filtering approach from collapsing of FASTQ and filtering to filtering
* fixed a bug to do with the RAM usage when reading a xeno file during subtraction


1.4.3
-----

*Release date: 2 October 2012*

* Fixes Error message when running subtraction and alignment is not found
* More work on fixing long running queries

1.4.2
-----

*Release date: 2 September 2012*

* Removes MongoKit dependency

1.4.1
-----

*Release date: 4 July 2012*

* Fixes cursor timeout for long running statistics queries

1.4.0
-----

*Release date: 1 June 2012*

* Adds statistics while they are being run rather than all at the end
* Saves unique id for genes as uid

1.3.0
-----

*Release date: 28 Feb 2012*

* Adds a new field in mapped for AS tag
* Links mapped genes to NCBI page
* Uses max sample coverage for project max coverage
* Adds timer to wait for file collapsing to finish when intersecting
* Using a path/to/file input no longer breaks intersecting

1.2.6
-----

*Release date: 6 Feb 2012*

* Subtraction filters out unmapped when building mapped reads

1.2
---

*Release date: 20 Jan 2012*

* Utility to create FastQ files from unmapped reads
* Utility to return intersection of FastQ files
* Utility to return filter FastQ files
* Added mapq threshold flag to subtraction
* Gbloader can now load bacteria and fungal genomes
* Saves genome sequences to database using GridFS
* Reduced memory footprint when calculating statistics
* use subprocesses instead of os.system
* mapq scores filter needs to include 0 <= mapq <=3 due to multiple mapping reads

1.1
---

*Release date: 25 Nov 2011*

* Adds support for pair-end reads

1.0.1
-----

*Release date: 19 Oct 2011*

* Added support for MongoDB authentication

1.0
---

*Release date: 12 Oct 2011*

* Initial release
