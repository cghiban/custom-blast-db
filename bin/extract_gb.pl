
use strict;
use Bio::SeqIO;
#use Bio::Seq::RichSeq;
use Data::Dumper;

# set if you're only interested in one gene/products
# undef if you want everything
my $FILTER_BY_PRODUCT = 'cytochrome oxidase subunit I';

my $input_file = $ARGV[0];
#print 'Working on ', $input_file, $/;

my $in = new Bio::SeqIO(-file => $input_file, -format=>"genbank");
while (my $seq = $in->next_seq()) {
    #print STDERR $seq->version, "\t", $seq->accession, "\n";
    #print $seq->id, "\t", $seq->accession, "\t", $seq->description, $/;
    my %data;
    for my $ft ($seq->get_SeqFeatures) {
        #print $ft->get_all_tags, $/;
        #print "  #len=", $ft->length, "\t", $ft->start, "..", $ft->end, $/;
        for my $tag (qw(db_xref organelle organism product)) {
            if ($ft->has_tag($tag)) {
                my @values = $ft->get_tag_values($tag);
                my $value = join '', @values;
                #print STDERR "  *", $tag, " ", $value, "\n";
                #next;
                if ($tag eq 'product') {
                    if ($value eq $FILTER_BY_PRODUCT) {
                        $data{$tag} = $value;
                        $data{seq} = $ft->seq->seq;
                        $data{desc} = $ft->seq->desc;
                        #print $ft->seq, $/;
                    }
                    else {
                        next;
                    }
                }
                else {
                    $data{$tag} = $value;
                }
                #print "\t$tag: ", $value, $/;
            }
            #print Dumper($ft), $/;
        }
    }
    #print Dumper(\%data), $/;
    my $strseq = $data{seq} || $seq->seq || '';
    my $seq_len = length ($strseq);
    #next if ($strseq eq '');
    next if ( $seq_len < 100 || $seq_len > 9000);
    #next if ( length($strseq) <= 2000 && $data{organism} !~ /Esox flaviae|Saccoglossus kowalevskii/ );

    my $taxonid = $data{db_xref}; $taxonid =~ s/^.*://;
    my $s = join "; ", (
        ">" . $seq->id . ($seq->version ? '.' . $seq->version :  ''), 
        "organism=$data{organism}",
        "taxon_id=$taxonid",
        "classification=".join(",", $seq->species->classification),
        $data{desc} || $seq->desc
        );
    #print STDERR length($strseq), "\t", $s, "\n" if length($strseq) > 2000;
    print $s, "\n", $strseq, "\n\n";
    
    #next;
    #last;
}
