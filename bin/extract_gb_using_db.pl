#!perl -w

use strict;
use Bio::SeqIO;
#use Bio::Seq::RichSeq;
use Mojo::SQLite;
use Data::Dumper;
use Cwd 'cwd';

# set if you're only interested in one gene/products
# undef if you want everything
#my $FILTER_BY_PRODUCT = 'cytochrome oxidase subunit I';
my $FILTER_BY_PRODUCT = '16S ribosomal RNA';


sub get_db {
    my ($file) = @_;
    my $db = Mojo::SQLite->new('sqlite:' . $file);
    $db->on(connection => sub {
        my ($db, $dbh) = @_;
        $dbh->do('pragma foreign_keys=ON');
        $dbh->do('pragma journal_mode=WAL');
        #$dbh->do('pragma journal_size_limit=1000000');
    });

    #print $db->db->query('select sqlite_version() as version')->hash->{version}, "\n";
    return $db->db;
}

sub search {
    my ($db, @accns) = @_;
    my $data = {};
    my @rez = $db->select('entries', [qw(file version)], {
            version => { -in => \@accns}
        })->hashes->each;
    for (@rez) {
        #print Dumper($_), "\n";
        if (exists $data->{ $_->{file}}) {
            push @{ $data->{ $_->{file}} }, $_->{version};
        }
        else {
            $data->{ $_->{file} } = [ $_->{version} ];
        }
    }

    $data;
}

#
#-- main
#

my $dbfile = $ARGV[0];

# this is the file with a list of accesstions
my $input_file = $ARGV[1]; 
print STDERR "\nWorking on $input_file\n";
#

my $db = get_db($dbfile);
my $cwd = cwd();

# XXX - better name for this scalar
# here we keep the .gb files and the accessions for which we need to get the fasta data
my %zzz = ();

my @accns = ();
open(my $fh, '<', $input_file) or die $!;
while (<$fh>) {
    chomp;
    next if /^$/;
    push @accns, $_;
}
close $fh;

#print STDERR "@accns\n";
my $dbdata = search($db, @accns);
#print STDERR Dumper($data), $/;
for my $gb (keys %$dbdata) {
    my %accns = map { $_ => 1} @{$dbdata->{$gb}};
    my $short_gb = $gb; $short_gb =~ s/^$cwd//;
    print STDERR "\t", $short_gb, "\n";
    #print STDERR " -- ", join (" ", @accns), "\n";
    #next;

    my $in = new Bio::SeqIO(-file => $gb, -format=>"genbank");
    while (my $seq = $in->next_seq()) {
        my $seqver = $seq->version ?  $seq->id . '.' . $seq->version : $seq->id;
        next unless (defined $accns{$seqver});
        #print STDERR ' -- ', $seqver, "\t", $seq->accession, "\n";
        #next;
        #print $seq->id, "\t", $seq->accession, "\t", $seq->description, $/;
        my %data;
        for my $ft ($seq->get_SeqFeatures) {
            #print $ft->get_all_tags, $/;
            #print "  #len=", $ft->length, "\t", $ft->start, "..", $ft->end, $/;
            for my $tag (qw(db_xref organelle organism product)) {
                if ($ft->has_tag($tag)) {
                    my @values = $ft->get_tag_values($tag);
                    my $value = join '', @values;
                    if ($tag eq 'product') {
                        if ($value eq $FILTER_BY_PRODUCT) {
                            $data{$tag} = $value;
                            $data{seq} = $ft->seq->seq;
                            $data{desc} = $ft->seq->desc;
                        }
                        else {
                            next;
                        }
                    }
                    else {
                        $data{$tag} = $value;
                    }
                }
                #print Dumper($ft), $/;
            }
        }
        my $strseq = $data{seq} || $seq->seq || '';
        my $seq_len = length ($strseq);
        #next if ($strseq eq '');
        next if ( $seq_len < 150 || $seq_len > 9000);

        my $taxonid = $data{db_xref}; $taxonid =~ s/^.*://;
        my $s = join "; ", (
            ">" . $seqver, 
            "organism=$data{organism}",
            "taxon_id=$taxonid",
            "classification=".join(",", $seq->species->classification),
            $data{desc} || $seq->desc
            );
        print $s, "\n", $strseq, "\n\n";
        
        #next;
        #last;
    } # end while
    #last;
} #end for
