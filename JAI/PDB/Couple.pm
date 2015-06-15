#!perl
################################################################################
## This script ....
## Author: ....
## Usage:
##   use JAI::PDB::Couple
################################################################################

################################################################################
## Module name
package JAI::PDB::Couple;
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
use Symbol;
################################################################################

################################################################################
## Exports
our @EXPORT_OK = qw(couple_models);
our %EXPORT_TAGS = ( DEFAULT => [qw(&couple_models)],
  Both => [qw(&couple_models)] );
################################################################################

################################################################################
## Constants
Readonly my $CONSTANT => 0;
################################################################################

################################################################################
## Version
our $VERSION = 1.0000;
################################################################################

################################################################################
## couple_models
sub couple_models {
  ##############################################################################
  ## Read function parameters
  my $pdb_file = shift;
  my $out_file = shift;
  ##############################################################################

  ##############################################################################
  ## Create a PDB data stream
  my $stream = Bio::Structure::IO->new( -file => $pdb_file,
    -format => 'PDB' );
  ##############################################################################

  ##############################################################################
  ## Iterate over all structures in the stream, usully just one.
  while ( my $structure = $stream->next_structure ) {
    my @models    = $structure->get_models();
    my $entry     = Bio::Structure::Entry->new( -id => 'default' );
    my $out_model = Bio::Structure::Model->new( -id => 'default' );
    $entry->add_model($out_model);
    for my $model (@models) {
      my @chains = $structure->get_chains($model);
      for my $chain (@chains) {
        my $out_chain = Bio::Structure::Chain->new( -id => $chain->id() );
        $entry->add_chain( $out_model, $out_chain );
        my @residues = $structure->get_residues($chain);
        for my $residue (@residues) {
          my $out_residue = Bio::Structure::Residue->new( -id => $residue->id() );
          $entry->add_residue( $out_chain, $out_residue );
          my @atoms = $structure->get_atoms($residue);
          for my $atom (@atoms) {
            $entry->add_atom( $out_residue, $atom );
          }
        }
      }
    }
    my $out_stream = Bio::Structure::IO->new(
      -file   => '>' . $out_file,
      -format => 'PDB' );
    $out_stream->write_structure($entry);
  }
  ##############################################################################
  return;
}
################################################################################
