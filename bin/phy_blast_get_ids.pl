#!/usr/bin/env perl

use strict;
use warnings;

use Carp ();
use FindBin;
use IO::File ();
use File::Spec ();
use XML::Simple ();
use Parallel::ForkManager ();
use LWP::UserAgent ();
use Data::Dumper;
use  Time::HiRes qw(sleep);

sub config {
    my ($file) = @_;
    local $@; # keep it to ourselves!
    my $data = do $file;
    if ($@) {
        Carp::carp "ERROR loading configuration from '$file': $@\n";
        return;
    } elsif (!defined $data) {
        Carp::carp "ERROR loading configuration from '$file': $!\n";
        return;
    }
    $data;
}

# return { webenv => a, querykey => b, count => c}
sub query_entrez {
    my ($cf) = @_;
    my $xml = get_entrez_data($cf, op => 'search');
	#print $xml, $/;
    my $xs = XML::Simple->new();
    my $ref = $xs->XMLin($xml);
	return $ref;
}

sub get_entrez_data {
    my ($cf, %params) = @_;

	my $ua = LWP::UserAgent->new;
	$ua->from('dnalcadmin@cshl.edu');

    my $SORTBY = 'MDAT'; # or ACCN

    my $query;
    if (exists $params{op} && $params{op} eq 'search') {
        # XXX user sort=ACCN
        my @terms = sort @{$cf->{TERMS}};
		$query = "esearch.fcgi?db=nuccore&usehistory=y&sort=$SORTBY"
			. '&term=(300:20000[SLEN])+AND+(' . join('+OR+', map {qq{"$_"}} @terms) . ')';
        if ($cf->{ONLY_VERTEBRATA}) {
            $query .= '+AND+Vertebrata[porgn:__txid7742]';
        }
        elsif ($cf->{ONLY_METAZOA}) {
            $query .= '+AND+Metazoa[porgn:__txid33208]';
        }
        if ($cf->{EXCLUDE_TERMS} && @{$cf->{EXCLUDE_TERMS}}) {
            $query .= '+NOT+' . join('+NOT+', map {$_} @{$cf->{EXCLUDE_TERMS}});
        }
    }
    elsif(exists $params{webenv} && exists $params{from}) {
		my $batch_size = $cf->{BATCH_SIZE} || 500;

        # rettype for nuccore: acc, fasta, seqid, native(xml)
        $query = 'efetch.fcgi?db=nuccore'
		    . '&rettype=acc&WebEnv=' . $params{webenv}
			. '&query_key=' . $params{qkey}
			. '&retstart=' . $params{from}
			. '&retmax=' . $batch_size;

		#$ua->default_header( 'Accept-Encoding' => 'gzip, deflate');
		#print STDERR  "+ batch_size: $batch_size\tfrom: $params{from}", 
		#	"\t", $params{file}, "\n";
    }

    my $url = $cf->{SERVICE} . '/' . $query;
    #print STDERR ' ~~ ', $query, $/;
    #exit 0;
	my $resp = (defined $params{file} && $params{file} ne '' && $params{from} ne '')
		? $ua->get($url, ':content_file' => $params{file})
		: $ua->get($url);
	#print STDERR  $resp, $/;
    my $content = '';
	if ($resp->is_success && !$params{file}) {
		$content = $resp->content;
	}
	else {
		#print STDERR  "+ fasta: ", $params{file}, $/;
		if (!$resp->is_success) {
			print STDERR ' ***: ', $url, $/;
			print STDERR ' ***: ', $resp->status_line , $/;
		}
		else {
			my $fsize = -s $params{file} || 0;
			if ($fsize < 1024) {
				print STDERR  ' *** file size: ', $fsize, $/;
			}
		}
		$content = '';
	}

    return ($resp->code, $content);
}

# main

my $cf = config("$FindBin::Bin/UTIL_PHY_BLAST");

# perform search
# http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term=mitochondrion[Filter]&usehistory=y
#   
# esearch.fcgi?db=<database>&term=<query1>&usehistory=y
# # esearch produces WebEnv value ($web1) and QueryKey value ($key1) 
#
# esearch.fcgi?db=<database>&term=<query2>&usehistory=y&WebEnv=$web1
# # esearch produces WebEnv value ($web2) that contains the results 
# #of both searches ($key1 and $key2)

mkdir $cf->{TMP};

#2. parse result, get WebEnv & QueryKey & Count
#   TODO store these values into a .query file so we can attempt 
#   to fetch files that didn't get downloaded (I see 500 and 502 errors quite often)

my $ref = {};
my $queryinfopath = File::Spec->catfile($cf->{TMP}, '.queryinfo');

if (-f $queryinfopath) {
	print STDERR  '*** Using stored query info: [', $queryinfopath, ']', $/;
	if (open(my $infh, '<', $queryinfopath)) {
		my $queryinfostring = '';
		local $_;
		while (<$infh>) {
			$queryinfostring .= $_;
		}
		#print STDERR "\n[", $queryinfostring, "]\n";
		eval $queryinfostring;
		close $infh;
	}
}
else {

	$ref = query_entrez($cf);
	if (open(my $outfh, '>', $queryinfopath)) {
		local $Data::Dumper::Purity = 1;
		#print STDERR Data::Dumper->Dump( [$ref], ['$ref']), "\n";
		print $outfh Data::Dumper->Dump( [$ref], ['$ref']), "\n";
		close $outfh;
	}
}

#my ($webenv, $count, $qkey) = ('NCID_1_123480816_130.14.18.34_9001_1443466331_23493941_0MetA0_S_MegaStore_F_1', '5316651', '1');

my ($webenv, $count, $qkey) = map {$ref->{$_}} qw/WebEnv Count QueryKey/;
print STDERR  "\nWebEnv:\t$webenv\nCount:\t$count\nQueryKey:\t$qkey\n\n";

#exit 0;

#3. Increment retstart by retmax (BATCH_SIZE) value each iteration till count is reached
# http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&WebEnv=NCID_1_243765038_130.14.22.215_9001_1395256759_1567063032&+query_key=1&rettype=fasta&retstart=0&retmax=500

# 4. save all the fasta data into a file

my $batch_size = $cf->{BATCH_SIZE} || 500;

# setting it to 0 with not fork (good for debugging)
my $pm = new Parallel::ForkManager(4);
$pm->run_on_start(
    sub { my ($pid, $args ) = @_;
        my ( $from, $count, $file_name ) = @$args;
        #print " -- NEXT $from/$count\tfile: $file_name\n";
    }
);
$pm->run_on_finish(
    sub { my ($pid, $exit_code, $args) = @_;
        my ( $from, $count, $file_name ) = @$args;
        my $size = int((-s $file_name) / 1024);
        print " ++ DONE $from:\t$file_name $size KB\n";
     }
);


#my ($from, $file_num) = (10_000_000, 20_000);
my ($from, $file_num) = (- $batch_size, 0);

while ($count >= $from) {
    #last if $file_num > 20_000;
    #

	my $file_name = File::Spec->catfile( $cf->{TMP}, sprintf("%d-%05d.txt", $batch_size, $file_num) );

	if ( -e $file_name ) {
		$file_num++;
		$from += $batch_size;
		next;
	}

	$file_num++;
	$from += $batch_size;

    $pm->start([$from, $count, $file_name]) and next; # do the fork

	my ($code, $fasta, $tries) = (0, '', 3);
	while ( $tries-- && $code != 200 ) {
            ($code, $fasta) = get_entrez_data(
                $cf, 
                webenv => $webenv, 
                qkey => $qkey, 
                from => $from,
                file => $file_name,
            );
			if ($code != 200) {
				print STDERR  " >>>>> HTTP_RESPONSE_CODE: ", $code, ' >>>>> tries left: ', $tries, $/;
				sleep 5;
			}
    }

    # no longer needed, GZIPPED output
	if ($fasta =~ /<ERROR>(.*)<\/ERROR>/) {
		print $1, "\n";
		sleep  3;
        $pm->finish; # do the exit in the child process
		#next;
	}
	#last if $from > 4000;
	#print STDERR  "\n";

	sleep 2;
    $pm->finish; # do the exit in the child process
} # END main while

$pm->wait_all_children;


