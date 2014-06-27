#!/usr/bin/env perl -w

use common::sense;

use Getopt::Long;
use MongoDB;
use JSON;

my $database_name = "capsidstaging";
my $database_host = "fervmongo1.res.oicr.on.ca";
my $database_port = 27017;
my $database_username = "capsid";
my $database_password = "mko0MKO)";

my $result = GetOptions(
	"database" => \$database_name,
	"host" => \$database_host,
	"port" => \$database_port,
	"username" => \$database_username,
	"password" => \$database_password
);

sub open_database {
    my (@args) = @_;
    
    my @options = (@args, host => "$database_host:$database_port", db_name => $database_name);
    push @options, username => $database_username if ($database_username);
    push @options, password => $database_password if ($database_password);
    push @options, query_timeout => -1;
    push @options, auto_reconnect => 1;

    my $conn = MongoDB::MongoClient->new(@options);
    $conn->connect();
    my $database = $conn->get_database($database_name);
    return $database;
}

sub close_database {
  my ($database) = @_;
}

sub process {
	my $database = open_database();

	my $collection = $database->get_collection("feature");
	my $query = $collection->find({type => "gene"});

	my %features = ();
	my @remove_queue = ();
	my $json = JSON->new()->convert_blessed(1);

	while (my $object = $query->next()) {
		my $key = "$object->{genome},$object->{start},$object->{end}";

		if (exists($features{$key})) {
			# We are about to remove something. So let's log what we are removing first
			say $json->encode($object);

			# And now, let's queue it for removal.
			push @remove_queue, $object->{_id};
		}

        $features{$key}++;
    }

    # Now for the nasty bit, try to remove.
    foreach my $id (@remove_queue) {

    	# We're dry-running here.
    	my $record = $collection->remove({_id => $id});
    	if (! defined($record)) {
    		die("Can't find record: $record")
    	}
    }

	close_database($database);
}

process();

1;