# This is -*-Perl-*- code
## Bioperl Test Harness Script for Modules
##
# $Id$

use strict;

BEGIN {
    use Bio::Root::Test;
    test_begin(-tests => 10, -requires_module =>'IO::String');


	use_ok('Bio::Tools::Run::Phylo::Hyphy::SLAC');
	use_ok('Bio::Tools::Run::Phylo::Hyphy::FEL');
	use_ok('Bio::Tools::Run::Phylo::Hyphy::REL');
	use_ok('Bio::Tools::Run::Phylo::Hyphy::Modeltest');
	use_ok('Bio::Tools::Run::Phylo::Hyphy::BatchFile');
}


SKIP: {
	my $rel = Bio::Tools::Run::Phylo::Hyphy::REL->new();

	test_skip(-requires_executable => $rel, -tests => 4);

	my $alignio = Bio::AlignIO->new(-format => 'fasta',
						 -file   => 't/data/hyphy1.fasta');

	my $treeio = Bio::TreeIO->new(-format => 'newick',
						 -file   => 't/data/hyphy1.tree');

	my $aln = $alignio->next_aln;
	my $tree = $treeio->next_tree;
	my $debug = test_debug();
	$debug = 1;
	my ($rc,$results);

#
#     print "rel:\n";
# 	my $rel = Bio::Tools::Run::Phylo::Hyphy::REL->new();
# 	$rel->alignment($aln);
# 	$rel->tree($tree);
# 	($rc,$results) = $rel->run();
# 	print ">>>" . $results . "\n";
# 	if (($rc == 0) && ($debug == 1)){
# 		warn("ERROR in REL module $rc:" . $rel->error_string() . "\n");
# 	}
# 	ok ($rc != 0, "REL module");
#
#     print "fel:\n";
# 	my $fel = Bio::Tools::Run::Phylo::Hyphy::FEL->new();
# 	$fel->alignment($aln);
# 	$fel->tree($tree);
#     $fel->save_tempfiles(1);
# 	($rc,$results) = $fel->run();
# 	if (($rc == 0) && ($debug == 1)){
# 		warn("ERROR in FEL module $rc:" . $fel->error_string() . "\n");
# 	}
# 	ok ($rc != 0, "FEL module");
#
#
#     print "modeltest:\n";
# 	my $modeltest = Bio::Tools::Run::Phylo::Hyphy::Modeltest->new();
# 	$modeltest->alignment($aln);
# 	$modeltest->tree($tree);
# #	$modeltest->save_tempfiles(1);
# 	($rc,$results) = $modeltest->run();
# #	print ">>>" . $results . "\n";
# 	if (($rc == 0) && ($debug == 1)){
# 		warn("ERROR in Modeltest module $rc:" . $modeltest->error_string() . "\n");
# 	}
# 	ok ($rc != 0, "Modeltest module");
#
    print "batchfile:\n";
	my $bf_exec = Bio::Tools::Run::Phylo::Hyphy::BatchFile->new(-params => {'bf' => "ModelTest.bf", 'order' => [$aln, $tree, '4', 'AIC Test', ""]});
	$bf_exec->alignment($aln);
	$bf_exec->tree($tree);
 	$bf_exec->save_tempfiles(1);
 	($rc,$results) = $bf_exec->run();
	if (($rc == 0) && ($debug == 1)){
		warn("ERROR in Batchfile module $rc:" . $bf_exec->error_string() . "\n");
	}
 	ok ($rc != 0, "Batchfile module");

    print "slac:\n";
	my $slac = Bio::Tools::Run::Phylo::Hyphy::SLAC->new();
	$slac->alignment($aln);
	$slac->tree($tree);
 	$slac->save_tempfiles(1);
 	print $slac->tempdir . "\n";
#	$slac->outfile_name("/Users/daisie/Desktop/slac.out");
	($rc,$results) = $slac->run();
# 	foreach my $i (keys (%$results)) {
#         print ">>>". $i . ": " . $results->{$i} . "\n";
#     }
	if (($rc == 0) && ($debug == 1)){
		warn("ERROR in SLAC module $rc:" . $slac->error_string() . "\n");
	}
	ok ($rc != 0, "SLAC module");
}
