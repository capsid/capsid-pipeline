#!/usr/bin/perl

=head1 NAME

readscan.pl - Program to identify pathogenic/contaminant reads in a host genome.

=head1 SYNOPSIS

B<readscan.pl> [B<subcommand>] [OPTION]... [FILE]..

B<subcommand> can be one of the following

=over 4

=item * stats

=back

=head1 DESCRIPTION

The program finds the presence of known pathogenic/contaminant sequences in the genome of the host.

=head1 readscan.pl stats 

Arguments

[B<SAM FILE>] 
sorted sam file.

=head1 PREREQUISITES

Readscan depends on 

=over 4

=item * perl

=item * Platform LSF L<http://en.wikipedia.org/wiki/Platform_LSF> 

=item * Unix utilities gzip,sort,cat etc. 

=back

=head1 AUTHOR

Written by Raeece Naeem

=head1 REPORTING BUGS

Report bugs to <raeecenaeem.mohamedghazzali@kaust.edu.sa>.

=head1 COPYRIGHT

Copyright (C) 2012, Raeece Naeem.
The file is licensed under the terms of the
GNU General Public License 3.
See <http://www.gnu.org/licenses/>.

=cut

	
use Getopt::Long ;
use File::Spec;
use strict;
use FindBin;
use Cwd 'abs_path';
use POSIX;

# name of the executing script;
my $program_name="readscan.pl";

our @arguments;
our $help;
our $sub_command;
our $usage = "
Usage : $program_name stats
\n";
GetOptions('help'=>\$help,'<>'=>,\&subcommand);
if(defined $help or $sub_command eq ""){
	print $usage;
}

sub subcommand{
	$sub_command=shift;

	if($sub_command eq "stats"){
		readscan_stats();
	}
	else{
		print "$usage\n";
	}
}

sub parse_cigar{
	my $cigar_str=shift;
	my @cigar_arr=();
	while($cigar_str=~ m/[0-9]+[MIDPSHN]/g){
		my $op=substr($&,length($&)-1,1);
		my $count=substr($&,0,length($&)-1);
		push @cigar_arr,[$op,$count];
	}
	return \@cigar_arr;
}

sub readscan_stats{
	my $i_sam_file;
	my $refstat_file;
	my $taxon_dir;
	my $cluster_file;
	my $help;
	my $verbose;
	my $data;
	my $gi_taxid_dmp;
	my $nodes_dmp;
	my $names_dmp;
	my $usage = "
		Usage : $program_name stats  -R <ref stat file> <input sam file>
		Required:-
		\t-R <file> : Refstat file
		\t-T <taxon directory> : Taxon directory where the ncbi taxon files gi_taxid_nucl.dmp, nodes.dmp and names.dmp are present
		Optional:-
		\t-C <file> : Cluster file
		\t-d : --data : data mode (tab separated file)
		\t-v : --verbose : verbose mode
		\n";
	GetOptions('Refstat=s'=>\$refstat_file,'Taxon=s'=>\$taxon_dir,'Cluster=s'=>\$cluster_file,'data' => \$data,'verbose' => \$verbose,'help'=>\$help,'<>'=>\&readscan_arguments);
	defined $help and  die $usage;
	$i_sam_file=@arguments[0];
# now validate all the parameters 
	(!defined $i_sam_file or !defined $refstat_file or !defined $taxon_dir)  and die $usage;
	if(! -e "$i_sam_file" || ! -e "$refstat_file"){
		die "please check if $i_sam_file and $refstat_file files are present\n";
	}
	if(! -d "$taxon_dir"){
		die "$taxon_dir does not exist\n";
	}
	else{
		my $abs_path= File::Spec->rel2abs($taxon_dir) ;
		$gi_taxid_dmp=$abs_path."/"."gi_taxid_nucl.dmp";
		$nodes_dmp=$abs_path."/"."nodes.dmp";
		$names_dmp=$abs_path."/"."names.dmp";
		if(! -e "$gi_taxid_dmp" || ! -e "$nodes_dmp" || ! -e "$names_dmp"){
			die "please check if $gi_taxid_dmp $nodes_dmp and $names_dmp files are present\n";
		}
	}
	my $sorted_sam="$i_sam_file"."_sorted";
	system("gunzip -c $i_sam_file | sort -T . -k1,1 -k3,3 > $sorted_sam");
	my $gra_ref=gra("$sorted_sam",$refstat_file);	
	system("rm $sorted_sam");
	my %gra=%{$gra_ref};
	open (ISAMFILE, "gunzip -c $i_sam_file |");
	my $read_count=0;
	my $sum_of_weights=0;
	my $contig_name="";
	my %uniq_read_names=();
	my $contig_start;
	my $contig_end;
	my $newreference=1;
	my $newcontig=1;
	my %ref_coverage=();
	my $sum_of_products=();
	my %ref_weighted_mean_contig_length=();
	my %ref_hit=();
	my $clength=0;
	while (<ISAMFILE>) {
		my $samline=$_;
		chomp;
		my ($read_name, $flag, $reference,$a,$b,$cigar,$c,$d,$e,$read,$quality,$aln_tag) = split("\t");
		$uniq_read_names{$read_name}=1;
		my ($x,$y,$aln_score)=split(":",$aln_tag);
		my $len=length($read);
		my $read_start=$a;
		my $read_end=$a+$len;
		my $last_softclipcount=0;
		my $first_softclipcount=0;
		my $cigar_arr=parse_cigar($cigar);
		my $counter=0;
		my $ins_count=0;
		my $del_count=0;
		for my $bits (@$cigar_arr){
			my $op=@$bits[0];
			my $count=@$bits[1];
			if($op eq "S"){
				if($first_softclipcount==0){
					$first_softclipcount=$count;
				}
				else{
					$last_softclipcount=$count;
				}
			}
			elsif($op eq "I"){
				$ins_count+=$count;
			} elsif($op eq "D"){
				$del_count+=$count;
			}
		}
# update the read_start and read_end with the first and last  occurrence of a match             
#print "read start = $read_start read end = $read_end\n";
		$read_end=$read_end-($first_softclipcount+$last_softclipcount+$ins_count)+$del_count;
#print "read start = $read_start read end = $read_end\n";
		$ref_hit{$reference}=$ref_hit{$reference}+1;
		$read_count++;
# check if the reference has changed
		if($reference eq $contig_name){
			$newreference=0;
		}
		else{
			$newreference=1;
		}
# if the current read intersects with previous contig extend contig and continue 
		if(!($newreference)){
			my @result=range_intersect($contig_start,$contig_end,$read_start,$read_end);
			if(@result[0]>0){
				$contig_start=@result[0];
				$contig_end=@result[1];
#print "chain $contig_start $contig_end  $read_start $read_end \n";
				next;
			}
			else{
				$newcontig=1;
			}
		}
# initialize for first read of a contig
		if($newreference || $newcontig){
# print previously  ended contig
			if($contig_name){
				$clength=$contig_end-$contig_start;
				$ref_coverage{$contig_name}=$ref_coverage{$contig_name}+($contig_end-$contig_start);
				$sum_of_weights+=$read_count;
				$sum_of_products=$sum_of_products+(($contig_end-$contig_start)*$read_count);

#print "contig ->\t$contig_name\t$contig_start\t$contig_end\tlen=$clength\treads=$read_count\n";
#print "contig -> $contig_name $contig_start $contig_end len=$clength reads=$read_count $sum_of_products $sum_of_weights\n";
				if($newreference){
					$ref_weighted_mean_contig_length{$contig_name}=$sum_of_products/$sum_of_weights;
#print "REF -> $contig_name $contig_start $contig_end len=$clength reads=$read_count $sum_of_products $sum_of_weights\n";
					$sum_of_products=0;
					$sum_of_weights=0;
				}

				$read_count=1;
			}
# initialize a new contig
			$contig_name=$reference;
			$contig_start=$read_start;
			$contig_end=$read_end;
			$newcontig=0;
		}
	}
	$clength=$contig_end-$contig_start;
	#print "contig  -> $contig_name $contig_start $contig_end len=$clength reads=$read_count\n";
	$ref_coverage{$contig_name}=$ref_coverage{$contig_name}+($contig_end-$contig_start);
	$sum_of_products=$sum_of_products+(($contig_end-$contig_start)*$read_count);
	$sum_of_weights+=$read_count;
	if($sum_of_weights==0){
		$ref_weighted_mean_contig_length{$contig_name}=0;
	}
	else{
		$ref_weighted_mean_contig_length{$contig_name}=$sum_of_products/$sum_of_weights;
	}
	close(ISAMFILE);
#for my $contig (keys %ref_coverage){
#print $ref_hit{$contig}."\t".$ref_coverage{$contig}."\t $contig \n";
#}

	if($refstat_file){
		open (REFFILE, $refstat_file);
# first line print the uniq number of reads
		my $uniq_number_of_reads=scalar(keys %uniq_read_names);
		#print "UNIQUE NUMBER OF READS = $uniq_number_of_reads\n";
		#printf "%5s\t%20s\t%20s\t%20s\t%5s\t%20s\t%s\n",'GRA','NO_OF_ALIGNS','BASES_COVERED','REF_LENGTH','PERC_COVERAGE','MEAN_CONTIG_LENGTH','REF_NAME';
		my %reffile_refname_length=();
		my %reffile_refname_reffullname=();
		while (<REFFILE>) {
			chomp;
			my @ref=split(/>/,$_);
			my @refname_samtoolslike=split(/ /,$ref[1]);
			$reffile_refname_length{$refname_samtoolslike[0]}=$ref[0];
			$reffile_refname_reffullname{$refname_samtoolslike[0]}=">".$ref[1];
		}
# FOREACH LOOP
		my %genome_perc_coverage=();
		#my %product_score=();
		my %ref_name=();
		foreach my $value (keys %ref_coverage){
			if(length($value)>=2){
				my $glen=$reffile_refname_length{$value};
				my $perc_coverage=($ref_coverage{$value}/$glen)*100;
				$genome_perc_coverage{$value}=$perc_coverage;
#$weighted_coverage{$value}=((($ref_coverage{$value}+($ref_sum_of_squares{$value}/$ref_coverage{$value}))/2)/$reffile_refname_length{$value})*100.0;
				#$product_score{$value}=($ref_coverage{$value}*$ref_hit{$value});
#print "$perc_coverage $ref_coverage{$value} $ref_hit{$value} @refname \n";
			}
		}

# ok everthing per sequence computed now group them in to taxon
my %gi_statsline=();
my %gi_gra=();


		if($cluster_file){
			my %cluster_score=();
			my %cluster_hash=();
			my $cluster_id=0;
			open(CLUSTERFILE,$cluster_file);
			while (<CLUSTERFILE>) {
				chomp;
				my @ref=split("\t",$_);
				my $max_coverage=0;
				my @sorted_ref=reverse (sort {$ref_weighted_mean_contig_length{$a} <=> $ref_weighted_mean_contig_length{$b}} @ref);
				$cluster_hash{$cluster_id}=\@sorted_ref;
				$cluster_score{$cluster_id}=$ref_weighted_mean_contig_length{$sorted_ref[0]};
				$cluster_id++;
			}		
			close(CLUSTERFILE);	
			foreach my $cluster_id (reverse (sort {$cluster_score{$a} <=> $cluster_score{$b} } keys %cluster_score)){
				my @ref=@{$cluster_hash{$cluster_id}};
				for my $value (@ref){
					printf "%20d\t%20d\t%20d\t%5.2f\t%10.1f\t%s\n",$ref_hit{$value},$ref_coverage{$value},$reffile_refname_length{$value},$genome_perc_coverage{$value},$ref_weighted_mean_contig_length{$value},$reffile_refname_reffullname{$value};
				}
				print "-----------------------------------------\n";
				print "-----------------------------------------\n";
				print "-----------------------------------------\n";
				print "-----------------------------------------\n\n";
			}
		}
		else{
			foreach my $value (reverse (sort {$gra{$a} <=> $gra{$b} } keys %gra)){
			my $stats_line=	sprintf "%20d\t%20d\t%20d\t%5.2f\t%10.1f\t%s\n",$ref_hit{$value},$ref_coverage{$value},$reffile_refname_length{$value},$genome_perc_coverage{$value},$ref_weighted_mean_contig_length{$value},$reffile_refname_reffullname{$value};
			my ($a,$gi,@dontbother)=split('\|',$value);
			$gi_statsline{$gi}=$stats_line;
			$gi_gra{$gi}=$gra{$value};
			}
		}
		close(REFFILE);
		# load the nodes file in to lookup
		my %nodes=();
		open(NODES,$nodes_dmp);
		while (<NODES>){
			my $nodes_line=$_;
			chomp;
			my ($node,$parent_node,$rank,@dont_bother)=split('\|');
		# remove leading / trailing spaces
			$node =~ s/^\s+|\s+$//g;
			$rank =~ s/^\s+|\s+$//g;
			$parent_node =~ s/^\s+|\s+$//g;
			$nodes{$node}="$parent_node"."-".$rank;
		}
		close(NODES);
		# load the node names into lookup
		my %node_names=();
		open(NODENAMES,$names_dmp);
		while (<NODENAMES>){
			my $nodes_names_line=$_;
			chomp;
			my ($node,$name,$a,$name_type,@dont_bother)=split('\|');
			$name_type =~ s/^\s+|\s+$//g;
			if($name_type eq "scientific name"){
				$node =~ s/^\s+|\s+$//g;
				$name =~ s/^\s+|\s+$//g;
				$node_names{$node}=$name;
			}
		}
		close(NODENAMES);
		# load gi - taxon into lookup
		# open (GI_TAXON,$gi_taxid_dmp);
		# my %taxon_list=();
		# # load all the taxon ids of the GIs into a list
		# while (<GI_TAXON>){
		# 	my $gi_taxon_line=$_;
		# 	chomp;
		# 	my ($gi,$taxon)=split(" ");
		# 	if(exists $gi_gra{$gi}){
		# 	$taxon_list{$gi}=$taxon;
		# 	}
		# }
		# close(GI_TAXON);	
		my %gi_group_taxon=();
		my %taxon_parent=();
		my @taxon_types=("species","genus","family","order","class","phylum");
		for my $gi (keys %gi_gra){
			my $taxon=get_taxon($gi, \%gi_gra, $taxon_dir);
			next if (! $taxon);
			my $rank="";
			my $parent_node="";
			my $leaf_taxon=$taxon;
			my %group_taxon=();
			my $exitsearch=1;
# lookup recursively until the root node is reached and collect all the group taxons
			do{
				($parent_node,$rank)=split("-",$nodes{$taxon});
# if the rank matches one of species,genus,family,order,class or phylum 
				if(grep (/$rank/,@taxon_types)){
					$group_taxon{$rank}=$taxon;
				}
				$taxon=$parent_node;
			}while($taxon!=1);

			# if(!exists($group_taxon{"species"})){
			# 	$group_taxon{"species"}=$leaf_taxon;
			# }
			# if(!exists($group_taxon{"genus"})){
			# 	$group_taxon{"genus"}=$group_taxon{"species"};
			# }
			# if(!exists($group_taxon{"family"})){
			# 	$group_taxon{"family"}=$group_taxon{"genus"};
			# }
			# if(!exists($group_taxon{"order"})){
			# 	$group_taxon{"order"}=$group_taxon{"family"};
			# }
			# if(!exists($group_taxon{"class"})){
			# 	$group_taxon{"class"}=$group_taxon{"order"};
			# }
			# if(!exists($group_taxon{"phylum"})){
			# 	$group_taxon{"phylum"}=$group_taxon{"class"};   
			# }

			my $rank = "sequence";
			my $current = $leaf_taxon;
			foreach my $parent_rank (@taxon_types) {
				if (exists($group_taxon{$parent_rank})) {
					$taxon_parent{$current."-".$rank} = $group_taxon{$parent_rank};
					$current = $group_taxon{$parent_rank};
				}
				$rank = $parent_rank;
			}

			# $taxon_parent{$leaf_taxon."-sequence"}=$group_taxon{"species"};
			# $taxon_parent{$group_taxon{"species"}."-species"}=$group_taxon{"genus"};
			# $taxon_parent{$group_taxon{"genus"}."-genus"}=$group_taxon{"family"};
			# $taxon_parent{$group_taxon{"family"}."-family"}=$group_taxon{"order"};
			# $taxon_parent{$group_taxon{"order"}."-order"}=$group_taxon{"class"};
			# $taxon_parent{$group_taxon{"class"}."-class"}=$group_taxon{"phylum"};

			for my $group (keys %group_taxon){
				$gi_group_taxon{$gi."-".$group}=$group."-".$group_taxon{$group}."-".$leaf_taxon;        
			}
		}
		# now summarise the results

		# but first, print a header in data mode
		if ($data) {
			print "level\tindex\tparent\tid\tname\tscore\n";
		}

		for my $taxon_type(@taxon_types){
			print uc($taxon_type)."\n" unless ($data);
			my $prefix = ($data) ? uc($taxon_type)."\t" : "";
			my $previous_group_taxon_id="";
			my $total_gra=0;
			my %group_total_gra=();
			for my $gi_group (sort {$gi_group_taxon{$a} cmp $gi_group_taxon{$b} } grep(/$taxon_type/,keys %gi_group_taxon)){
				my($gi,$group)=split("-",$gi_group);
				my($group,$group_taxon_id,$leaf_taxon_id)=split("-",$gi_group_taxon{$gi_group});
				if($previous_group_taxon_id eq "" or ($previous_group_taxon_id eq $group_taxon_id)){
					$total_gra+=$gi_gra{$gi};
				}
				else{
		# print the total_gra for the previous chunk    
					$group_total_gra{$previous_group_taxon_id}=$total_gra;
					$total_gra=$gi_gra{$gi};
				}
				$previous_group_taxon_id=$group_taxon_id;
			}
		# last chunk
			$group_total_gra{$previous_group_taxon_id}=$total_gra if ($previous_group_taxon_id ne "");

		# print descending order of total_gra
			my $number=1;
			for my $taxon(sort {$group_total_gra{$b} <=> $group_total_gra{$a}} keys %group_total_gra){
				my $taxon_parent_key=$taxon."-".$taxon_type;
				print $prefix.$number."\t".($taxon_parent{$taxon_parent_key} ? "ti:".$taxon_parent{$taxon_parent_key} : "")."\t"."ti:".$taxon."\t".$node_names{$taxon}."\t".$group_total_gra{$taxon}."\n";
				$number++;
			}
		}
		# first print the sequences sorted by gra
		print "SEQUENCE\n" unless ($data);
		my $prefix = ($data) ? "SEQUENCE\t" : "";
		my $number=1;
		for my $gi (sort {$gi_gra{$b} <=> $gi_gra{$a}} keys %gi_gra){
			my $taxon = get_taxon($gi, \%gi_gra, $taxon_dir);
			my $taxon_parent_key=$taxon."-sequence";
			if(! $taxon){
				print $prefix.$number."\t"."UNKNOWN TAXON  "."\t".$gi_gra{$gi}."\t".$gi_statsline{$gi};
			}
			else{
				print $prefix.$number."\t"."ti:".$taxon_parent{$taxon_parent_key}."\t"."gi:".$gi."\t".$node_names{$taxon}."\t".$gi_gra{$gi}."\t".$gi_statsline{$gi};
			}
			$number++;
		}

	}

}

sub readscan_arguments{
	my $argument=shift;
	push (@arguments,$argument);
}

use File::SortedSeek ':all';

my $taxon_index;

sub munge_taxon_line {
	my ($line) = @_;
	my ($gi) = split(/\t/, $line);
	return 0 + $gi;
}
 
sub get_taxon {
	my ($gi, $gi_gra, $taxon_dir) = @_;

	if (! defined($taxon_index)) {
		my $abs_path= File::Spec->rel2abs($taxon_dir);
		my $taxon_data_file = File::Spec->catfile($taxon_dir, "gi_taxid_nucl.dmp");
		open my $taxon_index_fh, "<", $taxon_data_file or die("Can't open $taxon_data_file: $!");
		$taxon_index = $taxon_index_fh;
	}

	if (exists $gi_gra->{$gi}){
		my $tell = File::SortedSeek::numeric($taxon_index, $gi, \&munge_taxon_line);
		my $line = <$taxon_index>;
		chomp($line);
		my ($found_gi, $found_taxon) = split(/\t/, $line);
		if ($gi == $found_gi) {
			return $found_taxon;
		}
	}
}

sub range_intersect{
	my ($s1,$e1,$s2,$e2)=($_[0],$_[1],$_[2],$_[3]);
	my @result;
	if (!($s1>$e2 || $s2>$e1)){
		my $min_val= ($s1<$s2) ? $s1 : $s2;
		my $max_val= ($e1>$e2) ? $e1 : $e2;
		@result=($min_val,$max_val);
	}
	else{
		@result=(-1);
	}
	return @result;
}

sub gra{

	my $i_sam_file=$_[0];
	my $refstat_file=$_[1];

	(!defined $i_sam_file or !defined $refstat_file)  and die "something wrong with $i_sam_file or $refstat_file";
	if(! -e "$i_sam_file" || ! -e "$refstat_file"){
		die "please check if $i_sam_file and $refstat_file files are present\n";
	}

	my $diff=0.01;
	my %glength=();
	open (REFFILE, $refstat_file);
	while (<REFFILE>) {
		chomp;
		my @ref=split(/>/,$_);
		my @refname_samtoolslike=split(/ /,$ref[1]);
		$glength{$refname_samtoolslike[0]}=$ref[0];
	}
	close(REFFILE);
	my %RGHit=();
	my %GHit=();
	my %rname=();
	my $total_hits=0;
	my $total_no_of_reads=0;
	my $previous_read_name="";
	my $previous_ref_name="";

# init phi
	my %phi=();
	open (ISAMFILE, $i_sam_file);
#print "init phi..\n";
#my $start = new Benchmark;
	while (<ISAMFILE>) {
		my $samline=$_;
		chomp;
		my ($read_name, $flag, $reference,$a,$b,$cigar,$c,$d,$e,$read,$quality,$aln_tag) = split("\t");
		if(!($previous_read_name eq $read_name)){
			$total_no_of_reads++;
		}
		if(!($previous_read_name eq $read_name && $previous_ref_name eq $reference)){
			$GHit{$reference}++;
		}
		$previous_read_name=$read_name;
		$previous_ref_name=$reference;
	}
	for my $genome (keys %GHit){
		$phi{$genome}=$GHit{$genome}/$total_no_of_reads;
	}
#my $end = new Benchmark;
#my $tdiff = timediff($end, $start);
#print "init phi complete..\n";
#print "Time taken was ", timestr($tdiff, 'all'), " seconds\n";
# compute phi_dash
	my %phi_dash=();
	my $i=1;
	my $flag=1;
	my $denominator_sum=0;
	do{
# re initialise phi_dash for the next iteration
#$start = new Benchmark;
		%RGHit=();
		%phi_dash=();
		$flag=1;
#print "iteration $i\n";
		$previous_read_name="";
# seek to first line 
		seek (ISAMFILE,0,0);
		while (<ISAMFILE>) {
			my $samline=$_;
			chomp;
			my ($read_name, $flag, $reference,$a,$b,$cigar,$c,$d,$e,$read,$quality,$aln_tag) = split("\t");
			if($previous_read_name eq $read_name){
				$RGHit{$reference}=(1/$glength{$reference});
			}
			else{
# process the previous chunk of of r->{G1,G2,G3,..}
				$denominator_sum=0;
				for my $ref (keys %RGHit){
					$denominator_sum+=($RGHit{$ref}*$phi{$ref});
				}
				for my $ref (keys %RGHit){
					my $prob_value=(($RGHit{$ref}*$phi{$ref})/$denominator_sum);
					$phi_dash{$ref}=$phi_dash{$ref}+=$prob_value;
				}
# re initialise RGHit for the next read 
				%RGHit=();
				if (! defined $glength{$reference}) {
					die("Internal error: can't find reference data: $reference");
				}
				$RGHit{$reference}+=(1/$glength{$reference});
			}
			$previous_read_name=$read_name;
		}
#process the last chunk
		$denominator_sum=0;
		for my $ref (keys %RGHit){
			$denominator_sum+=($RGHit{$ref}*$phi{$ref});
		}
		for my $ref (keys %RGHit){
			my $prob_value=(($RGHit{$ref}*$phi{$ref})/$denominator_sum);
			$phi_dash{$ref}=$phi_dash{$ref}+=$prob_value;
		}
# normalise phi_dash
		for my $genome (keys %phi_dash){
			$phi_dash{$genome}=$phi_dash{$genome}/$total_no_of_reads;
		}
# check if phi == phi_dash
		for my $genome (keys %phi){
			my $delta ||= $diff;  # default value of delta
				if(abs($phi{$genome}-$phi_dash{$genome})>$delta) {
					$flag=0;
					last;
				}
		}
# set phi = phi_dash
		for my $genome (keys %phi){
#print "$genome \t ".$phi{$genome}."\t".$phi_dash{$genome}."\n";
			$phi{$genome}=$phi_dash{$genome};
		}
		$i++;
#my $end = new Benchmark;
#my $tdiff = timediff($end, $start);
#print "Time taken for this iteration was ", timestr($tdiff, 'all'), " seconds\n";
	}while($flag==0);
	close(ISAMFILE);

	my %gra=();
	my $total_normalised_phi=0;
	for my $genome (keys %phi){
		$total_normalised_phi+=($phi{$genome}/$glength{$genome});
	}
	for my $genome (keys %phi){
		$gra{$genome}=($phi{$genome}/$glength{$genome})/($total_normalised_phi);
	}
	return \%gra;	
}

1;