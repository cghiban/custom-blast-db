#!/usr/bin/env perl

use strict;
use warnings;

use Carp ();
use FindBin;
use IO::File ();
use File::Spec ();
use XML::Simple ();
#use LWP::Simple qw(get);
use LWP::UserAgent ();
use Data::Dumper;

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
    return unless $xml;
	#print $xml, $/;
    my $xs = XML::Simple->new();
    my $ref = $xs->XMLin($xml);
	return $ref;
}

sub get_entrez_data {
    my ($cf, %params) = @_;

	my $ua = LWP::UserAgent->new;
	$ua->from('dnalcadmin@cshl.edu');

    my $query;
    if (exists $params{op} && $params{op} eq 'search') {
        # XXX user sort=ACCN
		$query = 'esearch.fcgi?db=nuccore&usehistory=y'
			. '&term=(300:20000[SLEN])+AND+(' . join('+OR+', map {qq{"$_"}} @{$cf->{TERMS}}) . ')';
        if ($cf->{ONLY_VERTEBRATA}) {
            $query .= '+AND+Vertebrata[porgn:__txid7742]';
            #$query .= '+NOT+"Canis lupus"[porgn]';
            #$query .= '+NOT+"Homo sapiens"[porgn]';
        }
        if ($cf->{EXCLUDE_TERMS}) {
            $query .= '+NOT+' . join('+NOT+', map {$_} @{$cf->{EXCLUDE_TERMS}});
        }
    }
    elsif(exists $params{webenv} && exists $params{from}) {
		my $batch_size = $cf->{BATCH_SIZE} || 500;

        $query = 'efetch.fcgi?db=nuccore'
		    . '&rettype=fasta&WebEnv=' . $params{webenv}
			. '&query_key=' . $params{qkey}
			. '&retstart=' . $params{from}
			. '&retmax=' . $batch_size;

		#$ua->default_header( 'Accept-Encoding' => 'gzip, deflate');
		print STDERR  "+ batch_size: $batch_size\tfrom: $params{from}", 
			"\t", $params{file}, "\n";
    }

    my $url = $cf->{SERVICE} . '/' . $query;
    print STDERR  $url, $/;
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
	if ($ref && open(my $outfh, '>', $queryinfopath)) {
		local $Data::Dumper::Purity = 1;
		#print STDERR Data::Dumper->Dump( [$ref], ['$ref']), "\n";
		print $outfh Data::Dumper->Dump( [$ref], ['$ref']), "\n";
		close $outfh;
	}
}

#my ($webenv, $count, $qkey) = ('NCID_1_123480816_130.14.18.34_9001_1443466331_23493941_0MetA0_S_MegaStore_F_1', '5316651', '1');

my ($webenv, $count, $qkey) = map {$ref->{$_}} qw/WebEnv Count QueryKey/;
print STDERR  "\nWebEnv:\t$webenv\nCount:\t$count\nQueryKey:\t$qkey\n\n";


#3. Increment retstart by retmax (BATCH_SIZE) value each iteration till count is reached
# http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&WebEnv=NCID_1_243765038_130.14.22.215_9001_1395256759_1567063032&+query_key=1&rettype=fasta&retstart=0&retmax=500

# 4. save all the fasta data into a file

my $batch_size = $cf->{BATCH_SIZE} || 500;

my ($from, $file_num) = (0, 0);
#my ($from, $file_num) = (1906000, 8975);

while ($count > $from) {
	my $file_name = File::Spec->catfile( $cf->{TMP}, sprintf("%d-%05d.fa", $batch_size, $file_num) );
	if (-e $file_name) {
		$file_num++;
		$from += $batch_size;
		next;
	}
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
				sleep 10;
			}
    }

    # no longer needed, GZIPPED output
	if ($fasta =~ /<ERROR>(.*)<\/ERROR>/) {
		print $1, "\n";
		sleep 3;
		next;
	}
	$file_num++;

	$from += $batch_size;
	#last if $from > 4000;
	print STDERR  "\n";
	sleep 2;
}

#5. get the taxonomy data (if not downloaded already)
# ftp.ncbi.nlm.nih.gov/pub/taxonomy
# ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/gi_taxid.readme
# ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/gi_taxid_nucl.zip
# ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/gi_taxid_prot.zip

#6. Format the file into a BLAST database using makeblastdb, you can encode taxid 
# using -taxid_map <File_In> option
#

#7. Organism name can be retrieved using blastdbcmd if you have taxdb files installed: ftp.ncbi.nlm.nih.gov/blast/db/
#-          More information on eutils is at: http://www.ncbi.nlm.nih.gov/books/nbk25501


__DATA__
<?xml version="1.0" ?>
<!DOCTYPE eSearchResult PUBLIC "-//NLM//DTD esearch 20060628//EN" "http://eutils.ncbi.nlm.nih.gov/eutils/dtd/20060628/esearch.dtd">
<eSearchResult><Count>2175191</Count><RetMax>20</RetMax><RetStart>0</RetStart><QueryKey>1</QueryKey><WebEnv>NCID_1_204219298_130.14.18.34_9001_1395325163_741127391</WebEnv><IdList>
<Id>594551835</Id>
<Id>594551833</Id>
<Id>594551356</Id>
<Id>594551354</Id>
<Id>594551352</Id>
<Id>594551350</Id>
<Id>594551348</Id>
<Id>594543235</Id>
<Id>594543221</Id>
<Id>594543161</Id>
<Id>594543147</Id>
<Id>594543139</Id>
<Id>594543138</Id>
<Id>594543137</Id>
<Id>594543136</Id>
<Id>594543135</Id>
<Id>594543134</Id>
<Id>594543133</Id>
<Id>594543132</Id>
<Id>594543131</Id>
</IdList><TranslationSet/><TranslationStack>   <TermSet>    <Term>barcode[Filter]</Term>    <Field>Filter</Field>    <Count>458621</Count>    <Explode>N</Explode>   </TermSet>   <TermSet>    <Term>mitochondrion[Filter]</Term>    <Field>Filter</Field>    <Count>2149635</Count>    <Explode>N</Explode>   </TermSet>   <OP>OR</OP>  </TranslationStack><QueryTranslation>barcode[Filter] OR mitochondrion[Filter]</QueryTranslation></eSearchResult>
