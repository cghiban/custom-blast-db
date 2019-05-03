#!/usr/bin/env perl

use strict;
use warnings;
use IO::File ();
use Bio::SeqIO ();
use Data::Dumper;
use Digest::SHA qw/sha256_hex/;
use Readonly;

Readonly my $SEQ_LEN_MIN => 300;
Readonly my $SEQ_LEN_MAX => 10_000;

# this is to filterout seqs with len outside of [300, 10_000]
# this is to filterout seqs that have the same DNA seq over 10 times
# (allows for seq to be repeated only 10 times)

my $dir = shift;

$dir =~ s/\/$// if defined $dir;

unless ($dir && -d $dir) {
	print STDERR "\nError: Dir not specified!\n",
		"Usage: $0 path/to/fasta/files\n\n";
	exit 1;
}


my %ids = ();
my %sha = ();
my $counter = 0;
my $okish = 0;
my $out = Bio::SeqIO->new(-fh => \*STDOUT , -format => 'Fasta');

#for my $file ("$dir/500-00000.fasta") {
for my $file (sort <$dir/*.fasta>) {
	print STDERR ' + ', $file, $/;
	#next;
	my $in = Bio::SeqIO->new(-file => $file , -format => 'Fasta');
	while ( my $seq = $in->next_seq ) {
		my $seq_id = $seq->display_id;

		$counter ++;

        if (0 && $seq_id =~ /^gi\|.*\|(?:gb|dbj|emb|ref|pdb)\|([^\s]+)\|/) {
            #print STDERR '    -- ', $1, "\t", ($ids{$1} || 1), "\t", $seq_id, "\n";
            $seq_id = $1;

            # if we've seen this seq id, discard it.
            next if ($ids{ $seq_id }++ > 1);

            # update the seq_id for this sequence
            $seq->display_id($seq_id);
        }

		#$ids{ $seq_id } ++;
		next if (++$ids{ $seq_id } > 1);


		my $slen = $seq->length;
		if ($slen < $SEQ_LEN_MIN || $slen > $SEQ_LEN_MAX) {
			#print STDERR '   filtered out: ', "\t", $slen, "\t", $seq_id, "\n";
			next;
		}

		my $desc = $seq->primary_seq->desc;
		next if $desc =~ /UNVERIFIED:|PREDICTED:/;

        if ($desc =~ / >/) {
            my @descs = split / >/, $desc;
            $desc = $descs[0] if @descs > 1;
            #print STDERR $seq_id, " ~~ @descs\n";
            #print STDERR $seq_id, " == $desc\n";
            $seq->desc($desc);
        }

		my $short_desc = join('', (split /\s+/, $desc)[0..3]);
		my $hkey = join('-', $short_desc, sha256_hex($seq->seq));
		if ($sha{ $hkey } ++ > 4) {
			print STDERR  '   filtered out: ', $hkey, ' ', sprintf("%3d", $sha{$hkey}), $/;
			next;
		}

		#if ( $ids{ $seq_id } > 1 ) {
		#	$seq->display_id( $seq_id . '|' . $ids{ $seq_id } );
		#	print 'alter: ', $seq->display_id, $/;
		#}
        print STDERR $seq->display_id, "\t", $seq->id, "\t", $short_desc, "\n";
		$out->write_seq($seq);
		$okish ++;
	}
}

print STDERR  'passed filtering: ', $okish, $/;
print STDERR  'total: ', $counter, $/;

