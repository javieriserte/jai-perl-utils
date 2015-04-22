#!perl
################################################################################
## This script ....
## Author: ....
## Usage:
##   use JAI::Cov::AddContacts
################################################################################

################################################################################
## Module name
package JAI::Cov::AddContacts;
################################################################################

################################################################################
## Define Modules Needed
use strict;
use warnings;
use version;
use IO::File;
use Getopt::Long;
use Bio::SeqIO;
use Readonly;
use autodie;
use IPC::Open3;
use Data::Dumper::Simple;
use Exporter 'import';
use JAI::PDB::Contacts;
use JAI::Cov::XmlRepr;
use JAI::Lists::Reduced qw(concat);
################################################################################

################################################################################
## Exports
our @EXPORT_OK = qw(add_contacts);
our %EXPORT_TAGS = ( DEFAULT => [qw(&add_contacts)],
  Both => [qw(&add_contacts)] );
################################################################################

################################################################################
## Constants
Readonly my $FALSE => 0;
Readonly my $TRUE  => 1;
################################################################################

################################################################################
## Version
our $VERSION = 1.0000;
################################################################################

################################################################################
## add_contacts
sub add_contacts {
  my $repr       = shift;
  my $executable = shift;
  my $chains = JAI::Cov::XmlRepr::get_tag_data( $repr, 'first_chain' ) . q(,) .
    JAI::Cov::XmlRepr::get_tag_data( $repr, 'second_chain' );
  my $pdb_file = 'tmp_pdb.pdb';
  JAI::Cov::XmlRepr::tag_data_to_file( $repr, 'pdb', $TRUE , $pdb_file);
  my $contact_mapper = JAI::PDB::Contacts->new($executable);
  $contact_mapper->set_criteria( 'closest,6+chains,' . $chains );
  $contact_mapper->show_chains();
  my @result        = $contact_mapper->calculate_contact_map($pdb_file);
  my $result_string = concat(@result);
  JAI::Cov::XmlRepr::add_tag( $repr, 'pdb_contacts', $result_string, $TRUE );
  return;
}
################################################################################
