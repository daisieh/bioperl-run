=head1 NAME

Bio::Tools::Run::Phylo::Hyphy::BatchFile - Wrapper around the Hyphy BatchFile analysis

=head1 SYNOPSIS

  use Bio::Tools::Run::Phylo::Hyphy::BatchFile;
  use Bio::AlignIO;
  use Bio::TreeIO;

  my $alignio = Bio::AlignIO->new(-format => 'fasta',
  			         -file   => 't/data/hyphy1.fasta');

  my $aln = $alignio->next_aln;

  my $treeio = Bio::TreeIO->new(
      -format => 'newick', -file => 't/data/hyphy1.tree');

  my $BatchFile = Bio::Tools::Run::Phylo::Hyphy::BatchFile->new();
  $BatchFile->alignment($aln);
  $BatchFile->tree($tree);
  my ($rc,$results) = $BatchFile->run();

=head1 DESCRIPTION

This is a wrapper around the BatchFile analysis of HyPhy ([Hy]pothesis
Testing Using [Phy]logenies) package of Sergei Kosakowsky Pond,
Spencer V. Muse, Simon D.W. Frost and Art Poon.  See
http://www.hyphy.org for more information.

This module will generate the correct list of options for interfacing
with TemplateBatchFiles/Ghostrides/Wrapper.bf.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Support

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and
reponsive experts will be able look at the problem and quickly
address it. Please include a thorough description of the problem
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via the
web:

  http://redmine.open-bio.org/projects/bioperl/

=head1 AUTHOR - Daisie Huang

Email daisieh@zoology.ubc.ca

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tools::Run::Phylo::Hyphy::BatchFile;
use vars qw(@ISA @VALIDVALUES $PROGRAMNAME $PROGRAM);
use strict;
use Bio::Root::Root;
use Bio::AlignIO;
use Bio::TreeIO;
use Bio::Tools::Run::Phylo::Hyphy::Base;
use Bio::Tools::Run::WrapperBase;

@ISA = qw(Bio::Root::Root Bio::Tools::Run::Phylo::Hyphy::Base);

=head2 Default Values

Valid and default values for BatchFile are listed below.  The default
values are always the first one listed.  These descriptions are
essentially lifted from the python wrapper or provided by the author.

INCOMPLETE DOCUMENTATION OF ALL METHODS

=cut

BEGIN {
    @VALIDVALUES =
        (
         {'geneticCode' => [ "Universal","VertebratemtDNA","YeastmtDNA","Mold/ProtozoanmtDNA",
                             "InvertebratemtDNA","CiliateNuclear","EchinodermmtDNA","EuplotidNuclear",
                             "Alt.YeastNuclear","AscidianmtDNA","FlatwormmtDNA","BlepharismaNuclear"]},
         {'tempalnfile' => undef }, # aln file goes here
         {'temptreefile' => undef }, # tree file goes here
        );
}

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Tools::Run::Phylo::Hyphy::BatchFile->new();
 Function: Builds a new Bio::Tools::Run::Phylo::Hyphy::BatchFile object
 Returns : Bio::Tools::Run::Phylo::Hyphy::BatchFile
 Args    : -alignment => the Bio::Align::AlignI object
           -save_tempfiles => boolean to save the generated tempfiles and
                              NOT cleanup after onesself (default FALSE)
           -tree => the Bio::Tree::TreeI object
           -params => a hashref of parameters (all passed to set_parameter)
           -executable => where the hyphy executable resides

See also: L<Bio::Tree::TreeI>, L<Bio::Align::AlignI>

=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  my ($aln, $tree, $st, $params, $exe,
      $ubl) = $self->_rearrange([qw(ALIGNMENT TREE SAVE_TEMPFILES
				    PARAMS EXECUTABLE)],
				    @args);
  defined $aln && $self->alignment($aln);
  defined $tree && $self->tree($tree);
  defined $st  && $self->save_tempfiles($st);
  defined $exe && $self->executable($exe);

  $self->set_default_parameters();
  if( defined $params ) {
      if( ref($params) !~ /HASH/i ) {
	  $self->warn("Must provide a valid hash ref for parameter -FLAGS");
      } else {
	  map { $self->set_parameter($_, $$params{$_}) } keys %$params;
      }
  }
  return $self;
}

=head2 update_ordered_parameters

 Title   : update_ordered_parameters
 Usage   : $BatchFile->update_ordered_parameters();
 Function: updates all of the parameters needed for the ordered input redirect in HBL.
 Returns : nothing
 Args    : none


=cut

sub update_ordered_parameters {
    my ($self) = @_;
	unless (defined ($self->{'_params'}{'order'})) {
		$self->throw("No ordered parameters for HYPHY were defined.");
	}
	for (my $i=0; $i< scalar @{$self->{'_params'}{'order'}}; $i++) {
		my $item = @{$self->{'_params'}{'order'}}[$i];
		if (ref ($item) =~ m/Bio::SimpleAlign/) {
			$item = $self->{'_params'}{'tempalnfile'};
		} elsif (ref ($item) =~ m/Bio::Tree::Tree/) {
			$item = $self->{'_params'}{'temptreefile'};
		}
		$self->{'_orderedparams'}[$i] = {$i, $item};
	}
	$self->SUPER::update_ordered_parameters();
}

=head2 run

 Title   : run
 Usage   : my ($rc,$results) = $BatchFile->run($aln);
 Function: run the BatchFile analysis using the default or updated parameters
           the alignment parameter must have been set
 Returns : Return code, Hash
 Args    : L<Bio::Align::AlignI> object,
	   L<Bio::Tree::TreeI> object [optional]


=cut

sub run {
	my ($self,$aln,$tree) = @_;

	unless (defined($self->{'_prepared'})) {
		$self->prepare($aln,$tree);
	}
	unless (defined($self->{'_params'}{'bf'})) {
		$self->throw("no custom batch file specified!");
	}
	my $rc = 1;
	my $results = "";
	my $commandstring;
	my $error_msg;
	my $output_str = "";
	my $BatchFileexe = $self->executable();
	my $outfile = $self->outfile_name;
	unless ($BatchFileexe && -e $BatchFileexe && -x _) {
		$self->throw("unable to find or run executable for 'HYPHY'");
	}

	my $kid_pid = open (KID, "-|"); # forks two processes; the child is for handling just the HYPHY exec.
	unless ($kid_pid) { # child
		#runs the HYPHY command
		$commandstring = $BatchFileexe . " BASEPATH=" . $self->program_dir . " " . $self->{'_wrapper'};
		my $pid = open(RUN, "$commandstring |") or $self->throw("Cannot open exe $BatchFileexe");

		my $waiting = waitpid $pid,0;
		# waitpid will leave a nonzero error in $? if the HYPHY command crashes, so we should bail gracefully.
		my $error = $? & 127;
		if ($error != 0) {
			print "Error: " . $self->program_name . " ($waiting) quit unexpectedly with signal $error";
			exit 0; # exit with 0 so that the parent knows HYPHY aborted.
		}
		#otherwise, return the results and exit with 1 so that the parent knows we were successful.
		while (my $line = <RUN>) {
			$results .= "$line";
		}
		close(RUN);
		print $results;
		exit 1;
	} else { # parent
		# read in the results from the child process
		while (my $line = <KID>) {
			$results .= "$line";
		}
		close KID;

		# process the errors from $? and set the error values.
		$rc = $? >> 8;
		if (($results =~ m/err/i) || ($rc == 0)) { # either the child process had an error, or HYPHY put one in the output.
			$rc = 0;
			$self->warn($self->program_name . " reported an error - see error_string for the program output");
			$results =~ m/(err.+)/is;
			$self->error_string($1);
		}
		open(OUTFILE, ">$outfile") or $self->throw("cannot open $outfile for writing");
		print OUTFILE $results;
		close(OUTFILE);

		unless ( $self->save_tempfiles ) {
			unlink($self->{'_wrapper'});
			$self->cleanup();
		}
	}

	return ($rc,$results);
}


=head2 create_wrapper

 Title   : create_wrapper
 Usage   : $self->create_wrapper
 Function: It will create the wrapper file that interfaces with the analysis bf file
 Example :
 Returns :
 Args    :


=cut

sub create_wrapper {
   my $self = shift;
   my $batchfile = '"' . $self->{'_params'}{'bf'} . '"';
   $self->SUPER::create_wrapper($batchfile);
}

=head2 set_parameter

 Title   : 	set_parameter
 Usage   : 	$hyphy->set_parameter($param,$val);
 Function: 	Sets the named parameter $param to $val if it is a non-numeric parameter
 			If $param is a number, sets the corresponding value of the ordered redirect array.
 Returns : 	boolean if set was success, if verbose is set to -1
           	then no warning will be reported
 Args    : 	$param => name of the parameter
           	$value => value to set the parameter to
 See also: 	L<no_param_checks()>

=cut


sub set_parameter {
	my ($self,$param,$value) = @_;
	if ($param =~ /\d+/) {
		$self->{'_params'}{'order'}[$param] = $value;
	} else {
		$self->{'_params'}{$param} = $value;
	}
	return 1;
}


=head2 set_default_parameters

 Title   : set_default_parameters
 Usage   : $BatchFile->set_default_parameters(0);
 Function: (Re)set the default parameters from the defaults
           (the first value in each array in the
	    @VALIDVALUES class variable)
 Returns : none
 Args    : boolean: keep existing parameter values


=cut

sub set_default_parameters {
	my ($self,$keepold) = @_;
	unless (defined $keepold) {
		$keepold = 0;
	}
	for (my $i=0; $i< scalar (@VALIDVALUES); $i++) {
		my $elem = @VALIDVALUES[$i];
		keys %$elem; #reset iterator
		my ($param,$val) = each %$elem;
		# skip if we want to keep old values and it is already set
		if (ref($val)=~/ARRAY/i ) {
			$self->{'_orderedparams'}[$i] = {$param, $val->[0]};
		} else {
			$self->{'_orderedparams'}[$i] = {$param, $val};
		}
	}
}


=head1 Bio::Tools::Run::WrapperBase methods

=cut

=head2 no_param_checks

 Title   : no_param_checks
 Usage   : $obj->no_param_checks($newval)
 Function: Boolean flag as to whether or not we should
           trust the sanity checks for parameter values
 Returns : value of no_param_checks
 Args    : newvalue (optional)


=cut

=head2 save_tempfiles

 Title   : save_tempfiles
 Usage   : $obj->save_tempfiles($newval)
 Function:
 Returns : value of save_tempfiles
 Args    : newvalue (optional)


=cut

=head2 tempdir

 Title   : tempdir
 Usage   : my $tmpdir = $self->tempdir();
 Function: Retrieve a temporary directory name (which is created)
 Returns : string which is the name of the temporary directory
 Args    : none


=cut

=head2 cleanup

 Title   : cleanup
 Usage   : $BatchFile->cleanup();
 Function: Will cleanup the tempdir directory after a run
 Returns : none
 Args    : none


=cut

=head2 io

 Title   : io
 Usage   : $obj->io($newval)
 Function:  Gets a L<Bio::Root::IO> object
 Returns : L<Bio::Root::IO>
 Args    : none


=cut

sub DESTROY {
    my $self= shift;
    unless ( $self->save_tempfiles ) {
	$self->cleanup();
    }
    $self->SUPER::DESTROY();
}

1;
