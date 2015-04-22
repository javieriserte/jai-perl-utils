#!perl
################################################################################
## This script ....
## Author: ....
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
################################################################################

################################################################################
## Exports
our @EXPORT_OK = qw(export_fasta_single export_fasta_multiple read_first
  id_map_from_fasta_text);
our %EXPORT_TAGS = ( DEFAULT => [qw(&export_fasta_single)],
  Both => [ qw(&export_fasta_single
      &export_fasta_multiple
      &read_first &id_map_from_fasta_text) ] );
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
