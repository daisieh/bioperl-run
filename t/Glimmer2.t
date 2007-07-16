#-*-Perl-*-
# ## Bioperl Test Harness Script for Modules
# #
use strict;
use vars qw($NTESTS);

BEGIN {
    $NTESTS = 217;
    
    eval { require Test::More; };
    if( $@ ) {
        use lib 't/lib';
    }
    use Test::More;
    
    plan tests => $NTESTS;
    use_ok('Bio::Tools::Run::Glimmer');
    use_ok('Bio::Root::IO');
    use_ok('Bio::Seq');
}

my $verbose = 1 if $ENV{'BIOPERLDEBUG'};

my $fasta_file = Bio::Root::IO->catfile('t','data','H_pylori_J99.fasta');
my $icm_file   = Bio::Root::IO->catfile('t','data','H_pylori_J99.glimmer2.icm');

my $factory    = Bio::Tools::Run::Glimmer->new(-program => 'glimmer2',
                                               -model => $icm_file);
isa_ok $factory, 'Bio::Tools::Run::Glimmer';

my $seqstream = Bio::SeqIO->new(-file => $fasta_file, -format => 'fasta');
my $seq = $seqstream->next_seq();

SKIP: {
    my $glimmer_present = $factory->executable();
    
    unless ($glimmer_present) {
        skip("glimmer2 program not found. Skipping tests 5 to $NTESTS.", ($NTESTS - 4));
    }
    
    my $glimmer2 = $factory->run($seq);
    isa_ok $glimmer2, 'Bio::Tools::Glimmer';
    
    my $first_gene = $glimmer2->next_prediction();
    isa_ok $first_gene, 'Bio::SeqFeature::Generic';
    is $first_gene->seq_id(), 'gi|15611071|ref|NC_000921.1|';
    
    while (my $gene = $glimmer2->next_prediction()) {
        isa_ok $gene, 'Bio::SeqFeature::Generic';
    }
}

1; 
