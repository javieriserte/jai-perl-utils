#!perl
################################################################################
## This script reads a pdb file and retrieves de amino acido sequence of one of
## the chains.
## Author: Javier Iserte
## Usage:
##   use Jai::PDB::Seq;
##   use Jai::PDB::Seq qw(get_chain_res gen_aa_map);
################################################################################

################################################################################
## Module name
package JAI::PDB::Seq;
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
use English qw( -no_match_vars );
use Data::Dumper::Simple;
use List::Util qw(reduce);
use Exporter 'import';
use warnings::register;
################################################################################

################################################################################
## Exports
our @EXPORT_OK = qw(get_chain_res gen_aa_map);
our %EXPORT_TAGS = ( DEFAULT => [qw(&get_chain_res)],
  Both => [qw(&get_chain_res &gen_aa_map)] );
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
## export_fasta
##   Writes a fasta sequente into a file
sub export_fasta {
  my $description = shift;
  my $sequence    = shift;
  my $outfile     = shift;
  open my $file_handler, '>', $outfile;
  exit if not print {$file_handler} '>' . $description . "\n" . $sequence . "\n";
  return close $file_handler;
}
################################################################################

################################################################################
## get_chain_res
##   Reads a pdb file, extract atom from a single chain and returns the
##   corresponding sequence. Missing residues are returned as 'X' or other
##   char (given by the user).
sub get_chain_res {
  my $pdb_file = shift;
  my $chain    = shift;
  my $missing  = shift;

  ##############################################################################
  ## Read pdb file
  open my $file_handler, '<', $pdb_file;
  my @model_data = get_chain_res_aux($file_handler);
  close $file_handler;
  ##############################################################################

  ##############################################################################
  ## Filter pdb text file lines
  ##   1. Lines, that contains ATOM data
  @model_data = grep { m/^ATOM/msx } @model_data;
  ##   2. Lines, contains from a given Chain
  @model_data = grep { m/^.{21}$chain/msx } @model_data;
  ##############################################################################
  ##############################################################################
  ## Get residue number, and number. And create a map with this data
  @model_data = map { get_chain_res_aux_2($ARG) } @model_data;
  my %model_res = map { $ARG->[1] => $ARG->[0] } @model_data;
  ##############################################################################

  ##############################################################################
  ## The the maximum residue number
  my $max_residue_number = reduce { $a > $b ? $a : $b } keys %model_res;
  ##############################################################################

  ##############################################################################
  ## Get three letter-code to one-letter code amino acid map.
  my $aa_map = gen_aa_map();
  ##############################################################################

  ##############################################################################
  ## Construct the aminoacid sequence. Missing residue are filled with a special
  ## char.
  my @seq = map {
    get_one_letter_code( $ARG, \%model_res, $aa_map, $missing )
  } 1 .. $max_residue_number;
  my $sequence = reduce { $a . $b } @seq;
  ##exists $model_res{$ARG} ? $aa_map->{$model_res{$ARG}} : $missing
  ##############################################################################

  return $sequence;
}

sub get_one_letter_code {
  my $index    = shift;
  my $residues = shift;
  my $map      = shift;
  my $missing  = shift;
  if ( exists $residues->{$index} ) {
    if ( exists $map->{ $residues->{$index} } ) {
      return $map->{ $residues->{$index} }
    } else {
      if ( warnings::enabled() ) {
        warnings::warn("Warning: A nonstandard three letter code for aminoacids was found: $residues->{$index}. Mapped to Missing: $missing.\n");
      }
      return $missing;
    }
  } else {
    return $missing;
  }
}

sub get_chain_res_aux {
  my $file_handler = shift;
  my @model_data   = ();
  ##############################################################################
  ## Get data from the first model in the PDB file
  my $model_found = 0;
  while ( my $line = <$file_handler> ) {
    if ( $line =~ m/^MODEL/msx ) {
      $model_found++;
    }
    if ( $model_found <= 1 ) {
      chomp $line;
      push @model_data, $line;
    } else {
      last;
    }
  }
  ##############################################################################
  return @model_data;
}

sub get_chain_res_aux_2 {
  my $line = shift;
  ##############################################################################
  ## get residue number and amino acid fields from PDB ATOM record,
  $line =~ s/^.{17}(...)..(....).+$/$1_$2/mgsx;
  $line =~ s/\s//msxg;
  my @res = split /_/msx, $line;
  ##############################################################################
  return \@res;
}
################################################################################

################################################################################
## gen_aa_map
##   Generates a map from amino acid written in three letter code to amino acids
##   written in one letter code. All chars are assumed to be in upper case.
sub gen_aa_map {
  my %result = ();
  $result{'GLN'} = 'Q';
  $result{'TRP'} = 'W';
  $result{'GLU'} = 'E';
  $result{'ARG'} = 'R';
  $result{'THR'} = 'T';
  $result{'TYR'} = 'Y';
  $result{'ILE'} = 'I';
  $result{'PRO'} = 'P';
  $result{'ALA'} = 'A';
  $result{'SER'} = 'S';
  $result{'ASP'} = 'D';
  $result{'PHE'} = 'F';
  $result{'GLY'} = 'G';
  $result{'HIS'} = 'H';
  $result{'LYS'} = 'K';
  $result{'LEU'} = 'L';
  $result{'CYS'} = 'C';
  $result{'VAL'} = 'V';
  $result{'ASN'} = 'N';
  $result{'MET'} = 'M';
  return \%result;
}
################################################################################

1;
