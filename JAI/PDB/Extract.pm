#!perl
################################################################################
## This script ....
## Author: ....
## Usage:
##   use JAI::PDB::Extract
################################################################################

################################################################################
## Module name
package JAI::PDB::Extract;
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
use Bio::Structure::IO::pdb;
use Bio::Structure::IO;
use Bio::Structure::Entry;
use English '-no_match_vars';
use Symbol;
################################################################################

################################################################################
## Exports
our @EXPORT_OK = qw(extract_models extract_chains);
our %EXPORT_TAGS = ( DEFAULT => [qw(&extract_models &extract_chains)],
  Both => [qw(&extract_models &extract_chains)] );
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
## split_models
##   Retrieve from a pdb the models that correspond to the ids given.
##   Each model goes to a different file.
sub split_models {
  my $pdb_file   = shift;
  my $model_ids  = shift;
  my $out_prefix = shift;
  my $stream     = Bio::Structure::IO->new( -file => $pdb_file,
    -format => 'PDB' );
  while ( my $structure = $stream->next_structure ) {
    my @models = $structure->get_models();
    my %model_map = map { $ARG->id() => $ARG } @models;
    for my $model ( @{$model_ids} ) {
      my $current_model = $model_map{$model};
      my $entry = Bio::Structure::Entry->new( -id => 'default' );
      my $out_model = Bio::Structure::Model->new( -id => $current_model->id() );
      $entry->add_model($out_model);
      my @chains = $structure->get_chains($current_model);
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
      my $out_stream = Bio::Structure::IO->new(
        -file   => '>' . $out_prefix . q(_) . $model . '.pdb',
        -format => 'PDB' );
      $out_stream->write_structure($entry);
    }
  }
  return;
}
################################################################################

################################################################################
## extract_models
##   Retrieve from a pdb the models that correspond to the ids given.
##   All models goes to the same file.
sub extract_models {
  my $pdb_file  = shift;
  my $model_ids = shift;
  my $out_file  = shift;
  my $stream    = Bio::Structure::IO->new( -file => $pdb_file,
    -format => 'PDB' );
  while ( my $structure = $stream->next_structure ) {
    my @models    = $structure->get_models();
    my %model_map = map { $ARG->id() => $ARG } @models;
    my $entry     = Bio::Structure::Entry->new( -id => 'default' );
    for my $model ( @{$model_ids} ) {
      my $current_model = $model_map{$model};
      my $out_model = Bio::Structure::Model->new( -id => $current_model->id() );
      $entry->add_model($out_model);
      my @chains = $structure->get_chains($current_model);
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
  return;
}
################################################################################

################################################################################
## extract_chains
##   Retrieve from a pdb the chains that correspond to the ids given.
##   All chain goes to the same file.
sub extract_chains {
  my $pdb_file   = shift;
  my $chains_ids = shift;
  my $out_file   = shift;
  my $stream     = Bio::Structure::IO->new( -file => $pdb_file,
    -format => 'PDB' );
  while ( my $structure = $stream->next_structure ) {
    my @models = $structure->get_models();
    my $entry = Bio::Structure::Entry->new( -id => 'default' );
    for my $model (@models) {
      my $out_model = Bio::Structure::Model->new( -id => $model->id() );
      $entry->add_model($out_model);
      my @chains                 = $structure->get_chains($model);
      my %chain_map              = map { $ARG->id() => $ARG } @chains;
      my @chains_to_be_extracted = ( defined $chains_ids ) ? @{$chains_ids}
        : map { $ARG->id() } @chains;
      for my $chain (@chains_to_be_extracted) {
        my $current_chain = $chain_map{$chain};
        my $out_chain = Bio::Structure::Chain->new( -id => $current_chain->id() );
        $entry->add_chain( $out_model, $out_chain );
        my @residues = $structure->get_residues($current_chain);
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
  return;
}
################################################################################
