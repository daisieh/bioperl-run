# $Id $
#
# BioPerl module for Bio::Tools::Run::Phylo::Phylip::Base
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Run::Phylo::Phylip::Base - DESCRIPTION of Object

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
email or the web:

  bioperl-bugs@bioperl.org
  http://bugzilla.bioperl.org/bioperl-bugs/

=head1 AUTHOR - Jason Stajich

Email jason@bioperl.org

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tools::Run::Phylo::Phylip::Base;
use vars qw(@ISA %DEFAULT);
use strict;

# Object preamble - inherits from Bio::Root::Root

use Bio::Root::Root;
use Bio::Tools::Run::WrapperBase;
@ISA = qw(Bio::Root::Root  Bio::Tools::Run::WrapperBase);

BEGIN {
    %DEFAULT = ( 'OUTFILE'   => 'outfile',
		 'TREEFILE'  => 'treefile',
		 );
}

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Tools::Run::Phylo::Phylip::Base();
 Function: Builds a new Bio::Tools::Run::Phylo::Phylip::Base object 
 Returns : an instance of Bio::Tools::Run::Phylo::Phylip::Base
 Args    :


=cut

=head2 outfile

 Title   : outfile
 Usage   : $obj->outfile($newval)
 Function: Get/Set default PHYLIP outfile name ('outfile' usually)
 Returns : value of outfile
 Args    : newvalue (optional)


=cut

sub outfile{
   my $self = shift;
   $self->{'_outfile'} = shift if @_;
   return $self->{'_outfile'} || $DEFAULT{'OUTFILE'}; 
}


=head2 treefile

 Title   : treefile
 Usage   : $obj->treefile($newval)
 Function: Get/Set the default PHYLIP treefile name ('treefile' usually)
 Returns : value of treefile
 Args    : newvalue (optional)


=cut

sub treefile{
   my $self = shift;
   $self->{'_treefile'} = shift if @_;
   return $self->{'_treefile'} || $DEFAULT{'TREEFILE'};
}



1;