#!perl -w

use strict;

use Mojo::Pg;
use Data::Dumper;
#use Time::HiRes qw( gettimeofday tv_interval );
use File::Basename;
use IO::Dir;
use IO::File;

my $pg = Mojo::Pg->new('postgresql://cornel@localhost/biosql2');
my $db = $pg->db;

sub get_sql_lineage {
    my ($tid) = @_;

    my $sql = q{
        WITH RECURSIVE recursetree(parent_id, rank, taxon_name, level) AS (
          SELECT taxon.parent_taxon_id, taxon.node_rank, taxon_name.name as lineage, 0 as level
          FROM taxon
          JOIN taxon_name ON (taxon.taxon_id = taxon_name.taxon_id AND taxon_name.name_class='scientific name')
          WHERE taxon.ncbi_taxon_id = ?
          
          UNION ALL
          
          SELECT t.parent_taxon_id, t.node_rank, tn.name as lineage, (rt.level + 1) as level
          FROM taxon t
          JOIN taxon_name tn ON (t.taxon_id = tn.taxon_id AND tn.name_class='scientific name')
          JOIN recursetree rt on t.taxon_id = rt.parent_id
          WHERE rt.taxon_name != 'root'
        )
        SELECT * FROM recursetree
        WHERE rank in ('superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')
        ORDER BY level desc
    };
    
    my @ranks = (qw( k p c o f g s));
    my %data = map {$_ => ''} @ranks;
    my $output = '';
    my $results = $db->query($sql, $tid);
    while ( my $r = $results->hash) {
        my $rank = $r->{rank}; $rank = 'kingdom' if $rank eq 'superkingdom';
        #print $$r{rank}, "\t", $rank;
        my $name = $r->{taxon_name};

        $data{ substr($rank, 0, 1) } = $$r{taxon_name};
    }
    $data{s} =~ s/^$data{g}\s+//;
    #print STDERR $tid, "\t", Dumper(\%data), "\n" if (scalar(keys %data) < 4);

    return join '; ', map { join '__', $_, $data{$_} } @ranks;
}

#print get_sql_lineage(309913);
#print get_sql_lineage(51367);
#print get_sql_lineage(1208065);
#print get_sql_lineage(9606);
#print get_sql_lineage(1425170);


# open fasta file
# read line by line
# if ^>
#   - get a new ID
#   - get it's taxon_id and get the lineage
#   - write lineage in taxonomy.txt
#   - write new seqId and it's sequence in a new files

sub process_input {

    # ofh - otus fh
    # tfh - taxonomy file handler
    #
    my ($fname, $cache, $ofh, $tfh, $counter, $skipped) = @_;
    my $infh = IO::File->new($fname);

    my ($seq_id, $org, $taxonid, $fulltx, $desc);

    my $passed_check;
    while (my $line = <$infh>) {
        if ($line =~ /^>/) {
            #if (defined $seq_len && $seq_len < $min_seq_len) {
            #    $min_seq_len = $seq_len;
            #}
            #if (defined $seq_len && $seq_len < 50) {
            #    print "\t", $seq_id, "\t", $taxonid, "\n";
            #}
            $passed_check = 0;
            #$seq_len = 0;

            $$counter ++;
            #($seq_id, $org, $taxonid, $fulltx) = $line =~ m/^>(\S+); organism=(.+?); taxon_id=(\d+);(?:\s+classification=(.*);)/;
            ($seq_id, $org, $taxonid, $fulltx, $desc) = split /; /, $line;
            (undef, $org) = split '=', $org;
            (undef, $taxonid) = split '=', $taxonid;
            #print STDERR  $line, $/ unless ($org);
            #print STDERR $taxonid, $/;
            $passed_check = 1 if $taxonid =~ /^\d+$/;
        }
        elsif ($passed_check) {
            #print STDERR $seq_id, "\t", $taxonid, $/;
            chomp $line;
            my $sha = Digest::SHA->new(256);
            $sha->add($_) for ($org, $taxonid, $line);
            my $hexsha = $sha->hexdigest;
            if (!exists $cache->{"$org-$hexsha"} ) {
                my $lineage = get_sql_lineage($taxonid);
                $$skipped ++ and next if ($lineage =~ /.__; .__; .__; /);
                $cache->{"$org-$hexsha"} ++;

                print $ofh ">", $$counter, "\n";
                print $ofh $line, "\n\n";
                print $tfh $$counter, "\t", $lineage, "\n";
            }
            else {
                $$skipped ++;
            }

            $passed_check = 0;
            #last if $$counter > 20_000;
        }
    }

    #if ($seq_len < $min_seq_len) {
    #    $min_seq_len = $seq_len;
    #}
    #print '$len=', $min_seq_len, "\n";
}

#-----------------------------------------------------------
# main

my $dirname = shift;

my $ofname= File::Spec->catfile('output', 'fs_otus.txt');
my $otname= File::Spec->catfile('output', 'fs_otu_taxonomy.txt');

my $ofh = IO::File->new($ofname, 'w');
my $tfh = IO::File->new($otname, 'w');

my $skipped = 0;
my $counter = 0;
my $cache = {};
my $dir = IO::Dir->new($dirname);
if (defined $dir) {
    while (defined($_ = $dir->read)) {
        next unless /\.fasta$/;
        my $f = File::Spec->catfile($dirname, $_);
        print $f, "\n";
        process_input($f, $cache, $ofh, $tfh, \$counter, \$skipped);
    }
}

print 'processed ', $counter, ' sequences', "\n";
print 'skipped ', $skipped, ' sequences', "\n";

my $max = 0;
for (sort {$cache->{$b} <=> $cache->{$a} } keys %$cache) {
    print " ++ hexsha: ", $_, "\t", $cache->{$_}, "\n" if $cache->{$_} > 1;
    last if $max++ >= 20;
}

__END__
my $sp = join(' ', @ARGV) || 'Homo sapiens';
my $t0 = [gettimeofday];

#my $db = Bio::DB::Taxonomy->new(-source => 'entrez');
my $db = Bio::DB::Taxonomy->new(
                -source => 'flatfile',
                -nodesfile => 'taxdata/nodes.dmp',
                -namesfile => 'taxdata/names.dmp'
        );

print '* ', tv_interval($t0, [gettimeofday]);

sub get_lineage_by_taxon_id {
    my ($tid) = @_;
    my $taxon = $db->get_taxon(-taxonid => $tid);
    print STDERR $taxon, $/;
    return;
    my $sp = '';
    # use NCBI Entrez over HTTP
    #my $taxonid = $db->get_taxonid('Homo sapiens');

    my $t1 = [gettimeofday];

    # get a taxon
    #my $taxon = $db->get_taxon(-taxonid => $taxonid);
    my $tree = $db->get_tree($sp);
    my $t2 = [gettimeofday];
    print '* ', tv_interval($t1, $t2);


    my $level = 0;
    my @descendents = $tree->get_root_node->get_all_Descendents;
    for (@descendents) {
        next if $_->rank eq 'no rank';
        print '  ' x $level++, $_->id, '  ', $_->rank, "\t", $_->scientific_name;
        #if ($_->division) {
        #    print "division is ", $_->division;
        #}
    }
    my $t3 = [gettimeofday];
    print '* ', tv_interval($t1, $t3);
}

sub get_lineage {
    my ($sp) = @_;
    # use NCBI Entrez over HTTP
    #my $taxonid = $db->get_taxonid('Homo sapiens');

    my $t1 = [gettimeofday];

    # get a taxon
    #my $taxon = $db->get_taxon(-taxonid => $taxonid);
    my $tree = $db->get_tree($sp);
    my $t2 = [gettimeofday];
    print '* ', tv_interval($t1, $t2);


    my $level = 0;
    my @descendents = $tree->get_root_node->get_all_Descendents;
    for (@descendents) {
        next if $_->rank eq 'no rank';
        print '  ' x $level++, $_->id, '  ', $_->rank, "\t", $_->scientific_name;
        if ($_->division) {
            print "division is ", $_->division;
        }
    }
    my $t3 = [gettimeofday];
    print '* ', tv_interval($t1, $t3);
}

get_lineage_by_taxon_id(9606);

__END__

get_lineage('Homo sapiens');
print "#----------------------\n";
get_lineage('Thermofilum librum');
print "#----------------------\n";
get_lineage('Zea mays');
print "#----------------------\n";

__DATA__



WITH RECURSIVE recursetree(parent_id, rank, taxon_name, level) AS (
    SELECT taxon.parent_taxon_id, 
	taxon.node_rank, taxon_name.name as lineage,
	0 as level
    FROM taxon
    JOIN taxon_name ON (taxon.taxon_id = taxon_name.taxon_id AND taxon_name.name_class='scientific name')
    WHERE taxon.ncbi_taxon_id=1246958
  UNION ALL
    SELECT t.parent_taxon_id, 
        t.node_rank, tn.name as lineage,
        (rt.level + 1) as level
    FROM taxon t
    JOIN taxon_name tn ON (t.taxon_id = tn.taxon_id AND tn.name_class='scientific name')
    JOIN recursetree rt on t.taxon_id = rt.parent_id
    WHERE rt.taxon_name != 'root'
  )
SELECT * FROM recursetree
--WHERE rank in ('kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species');

