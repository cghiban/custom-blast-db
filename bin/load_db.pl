#!perl -w

use strict;
use Bio::SeqIO;
use Mojo::SQLite;
use File::Spec;
use Cwd;
use Time::HiRes qw( usleep gettimeofday tv_interval);
use Data::Dumper;

sub get_db {
    my ($file) = @_;
    my $db = Mojo::SQLite->new('sqlite:' . $file);
    $db->on(connection => sub {
        my ($db, $dbh) = @_;
        $dbh->do('pragma foreign_keys = ON');
        $dbh->do('pragma journal_mode=WAL');
        #$dbh->do('pragma journal_size_limit=1000000');
    });

    $db->db->query(
        'create table if not exists entries (
            id rowid,
            accession text not null,
            version text not null,
            organism text not null,
            modified_on date not null,
            length integer not null,
            file text not null,

            unique (version)
        )'
    );

    #print $db->db->query('select sqlite_version() as version')->hash->{version}, "\n";
    return $db->db;
}

#-- main

my $dbfile = shift;
my $gbfile = shift;

my $db = get_db($dbfile);

die  "Missing input .gb file...\n" unless ($gbfile && -f $gbfile);

my $basedir = Cwd::cwd();

my @new_entries = ();

my $in = new Bio::SeqIO(-file => $gbfile, -format=>"genbank");
while (my $seq = $in->next_seq) {

    my $version = $seq->id . ($seq->version ? '.' . $seq->version :  '');

    my $strseq = $seq->seq || '';
    next if (length ($strseq) < 100);

    #print length($strseq), "\t", $version, "\t", 
    #    $seq->species->binomial, $/;
    my @dates = $seq->get_dates;

    #print $seq->as_string, $/;
    #last;
    #
    push @new_entries, {
                accession => $seq->id,
                version => $version,
                organism => $seq->species->binomial,
                modified_on => $dates[0],
                length => length($strseq),
                file => File::Spec->rel2abs($gbfile, $basedir),
            };
}

#print ' - about to add ', scalar( @new_entries ), " new entries!\n";


my %entries = map { $_->{version} => 1} 
    $db->select('entries', ['version'], {
        version => { -in => [map {$_->{version}} @new_entries ] }
    })->hashes->each;
#print Dumper(\%entries), $/;
#__END__

srand((localtime)[0]);
usleep(1 + rand(3));

if (@new_entries) {
    my $_t0 = [gettimeofday];
    my $cnt = 0;

    my $tx = $db->begin;
    for my $data (@new_entries) {
        next if (exists $entries{$data->{version} });
        #print ' new new: ', $data->{version}, "\n";
        my $e = $db->insert('entries', $data);
        $cnt ++;
    }
    $tx->commit;
    my $_telapsed = tv_interval ( $_t0 );

    print ' - added ', $cnt, " entries in ", $_telapsed, "s\n";
    print ' - ', $db->query('select count(*) as num from entries')->hash->{num}, " entries!\n";
}

#print ' - ', $db->query('select count(*) as num from entries')->hash->{num}, " entries!\n";

