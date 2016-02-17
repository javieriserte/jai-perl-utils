#!perl
################################################################################
## This module contains many functiosn to read and write fasta files
## Author: Javier Iserte
## Usage:
##   use JAI::Fasta::Io
################################################################################

################################################################################
## Module name
package JAI::Fasta::Io;
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
use List::MoreUtils;
use Exporter 'import';
use English qw(-no_match_vars);
################################################################################

################################################################################
## Exports
our @EXPORT_OK = qw(export_fasta_single export_fasta_multiple read_first
  id_map_from_fasta_text load_many_files_as_map export_fasta_map);
our %EXPORT_TAGS = (
  DEFAULT => [ qw(&export_fasta_single
      &load_many_files_as_map
      &export_fasta_map) ],
  Both => [ qw(&export_fasta_single
      &export_fasta_multiple
      &read_first
      &id_map_from_fasta_text
      &load_many_files_as_map
      &export_fasta_map) ] );
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
## read_first
sub read_first {
  my $file = shift;
  my $seqio_obj = Bio::SeqIO->new( -file => $file,
    -format => 'fasta' );
  if ( my $first = $seqio_obj->next_seq() ) {
    return ( $first->id, $first->seq() );
  }
  return;
}
################################################################################

################################################################################
## load_many_files_as_map : Gets a list of fasta files, reads them all, and
## returns a hash with all of them
sub load_many_files_as_map {
  my $seqs = {};
  my $desc;
  for my $file (@ARG) {
    open my $fh, '<', $file;
    while ( my $line = <$fh> ) {
      chomp $line;
      load_many_files_as_map_aux( $line, \$desc, $seqs );
    }
    close $fh;
  }
  return $seqs;
}

sub load_many_files_as_map_aux {
  my $line = shift;
  my $desc = shift;
  my $seqs = shift;
  if ( ( substr $line, 0, 1 ) eq q(>) ) {
    ${$desc} = substr $line, 1;
  } else {
    if ( defined $desc ) {
      if ( !exists $seqs->{ ${$desc} } ) {
        $seqs->{ ${$desc} } = q();
      }
      $seqs->{ ${$desc} } .= $line;
    }
  }
  return;
}
################################################################################

################################################################################
## export_fasta_single
sub export_fasta_single {
  my $description = shift;
  my $sequence    = shift;
  my $outfile     = shift;
  open my $file_handler, '>', $outfile;
  exit if not print {$file_handler} '>' . $description . "\n" . $sequence . "\n";
  return close $file_handler;
}
################################################################################

################################################################################
## export_fasta_multiple
sub export_fasta_multiple {
  my $description = shift;
  my $sequence    = shift;
  my $outfile     = shift;
  my $n_desc      = scalar @{$description};
  my $n_seq       = scalar @{$sequence};
  if ( $n_seq != $n_desc ) {
    return;
  } else {
    open my $file_handler, '>', $outfile;
    for my $i ( 0 .. $n_seq - 1 ) {
      exit if not print {$file_handler} '>' . $description->[$i] . "\n" .
        $sequence->[$i] . "\n";
    }
    return close $file_handler;
  }
}
################################################################################

################################################################################
## export_fasta_map : Gets a hash representing a collection of sequences and
## writes them to a file in fasta format.
sub export_fasta_map {
  my $sequences = shift;
  my $outfile   = shift;
  open my $fh, '>', $outfile;
  for my $desc ( sort keys %{$sequences} ) {
    exit if !print {$fh} ">$desc\n" . $sequences->{$desc} . "\n";
  }
  close $fh;
  return;
}
################################################################################

################################################################################
## id_map_from_fasta_text
##   Reads a string containg fasta formatted text. Return a map from sequence id
##   to sequence.
sub id_map_from_fasta_text {
  my $fasta_string = shift;
  my %result       = ();
  my $seqio_obj    = Bio::SeqIO->new(
    -string => $fasta_string,
    -format => 'fasta' );
  while ( my $seq_io = $seqio_obj->next_seq() ) {
    $result{ $seq_io->id } = $seq_io->seq();
  }
  return %result;
}
################################################################################
1;
