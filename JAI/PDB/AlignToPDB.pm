#!perl
################################################################################
## This script allows to align the sequence of a pdb file with a given fasta
##   sequence. The aligning is made with mafft.
## Author: Javier Iserte
## Usage:
##   align_seq_to_pdb.pl <options>
##   use JAI::PDB::AlignToPDB
## Options:
##   -pdb=s           Path to the input pdb file.
##   -chain=s         The chain in the pdb whose sequence will be aligned.
##   -seq=s           Path to the input fasta file. If file contains multiple
##                     sequences, the first one is used. Sequences are assumed
##                     to be non-gapped.
##   -out=s           Path to store the generated MSA.
##   -map=s           Optional. Path to store a map the refers positions of the
##                     input sequence that corresponds to positions of the pdb
##                     file using pdb numeration.
##   -help            Shows this help.
################################################################################

################################################################################
## Module name
package JAI::PDB::AlignToPDB;
################################################################################

################################################################################
## Define Modules Needed
use strict;
use warnings;
use version;
use IO::File;
use Bio::SeqIO;
use Readonly;
use autodie;
use IPC::Open3;
use Symbol;
use List::Util;
use JAI::PDB::Seq 'get_chain_res';
use JAI::Fasta::Io qw(read_first id_map_from_fasta_text);
use JAI::Lists::Reduced 'concat';
use IO::Select;
use Data::Dumper::Simple;
use Exporter 'import';
use Symbol;
################################################################################

################################################################################
## Exports
our @EXPORT_OK   = qw(align_to_pdb);
our %EXPORT_TAGS = ( DEFAULT => [qw(&align_to_pdb)],
                     Both    => [qw(&align_to_pdb)]);
################################################################################

################################################################################
## Constants
Readonly my $BUFFER_SIZE => 10;
################################################################################

################################################################################
## Version
our $VERSION = 1.0000;
################################################################################

################################################################################
###                                Functions                                ####
################################################################################
sub align_to_pdb {
  ##############################################################################
  ## Read function parameters
  my $pdb_file = shift;
  my $seq_file = shift;
  my $out_file = shift;
  my $chain = shift;
  my $map_file = shift;
  ##############################################################################

  ##############################################################################
  ## Read input fasta sequence
  my @fasta_seq = read_first($seq_file);
  ##############################################################################

  ##############################################################################
  ## Get sequence from pdb file
  my $pdb_seq = get_chain_res( $pdb_file, $chain, 'X' );
  ##############################################################################

  ##############################################################################
  ## align sequences
  my @aln_sequences = align( $fasta_seq[1], $pdb_seq );
  my $seq_string = concat(@aln_sequences);
  ##############################################################################

  ##############################################################################
  ## export msa
  send_to_msa( $out_file, $seq_string );
  ##############################################################################

  ##############################################################################
  ## export seq to pdb position map
  if ( defined $map_file ) {
    my %seqs = id_map_from_fasta_text($seq_string);
    my ( $seq_i, $seq_j ) = ( 0, 0 );
    my %map = ();
    for my $i ( 0 .. ( length $seqs{'pdb'} ) - 1 ) {
      my $pdb_is_not_gap   = ( substr $seqs{'pdb'},   $i, 1 ) ne q(-);
      my $pdb_is_not_x     = ( substr $seqs{'pdb'},   $i, 1 ) ne q(X);
      my $fasta_is_not_gap = ( substr $seqs{'fasta'}, $i, 1 ) ne q(-);
      $seq_i = $seq_i + ( $pdb_is_not_gap           ? 1 : 0 );
      $seq_j = $seq_j + ( $fasta_is_not_gap ne q(-) ? 1 : 0 );
      if ( $pdb_is_not_gap && $fasta_is_not_gap && $pdb_is_not_x ) {
        $map{$seq_j} = $seq_i;
      }
    }
    open my $file_hander, '>', $map_file;
    foreach my $fasta_pos ( sort { $a <=> $b } keys %map ) {
      exit if not print {$file_hander} $fasta_pos . "\t" .
        $map{$fasta_pos} . "\n";
    }
    close $file_hander;
  }
  ##############################################################################
  return;
}

################################################################################
## align
##   Align pdb and fasta sequences using mafft.
sub align {
  my $fasta_seq = shift;
  my $pdb_seq   = shift;
  my $outfile   = shift;

  my @results = ();

  ##############################################################################
  ## Create command line to align and sequence data to be aligned
  my $cmdline  = 'mafft -';
  my $seq_data = q(>) . 'fasta' . "\n" . $fasta_seq . "\n" .
    q(>) . 'pdb' . "\n" . $pdb_seq;
  ##############################################################################

  ##############################################################################
  ## Launch aligning process
  my ( $in, $out, $err );
  $err = gensym();
  my $pid = open3( $in, $out, $err, $cmdline );
  ##############################################################################

  ##############################################################################
  ## write sequence data to aligning process via stdin
  exit if not print {$in} $seq_data;
  close $in;
  sleep 1;
  ##############################################################################

  ##############################################################################
  ## Read output from the aligning process.
  my $sel = IO::Select->new();
  $sel->add( $out, $err );
  ## Get those buffers that can be read
  while ( my @ready = $sel->can_read ) {
    ## Iterate over each readable buffer
    foreach my $fh (@ready) {
      my $text;
      ## Read from current buffer
      my $len = sysread $fh, $text, $BUFFER_SIZE;
      ## Check if the read was succesful
      if ( !defined $len ) {
        die "Error from child: $pid\n";
        ## Check if the buffer was entirely read, if so close it.
      } elsif ( $len == 0 ) {
        $sel->remove($fh);
        next;
        ## Otherwise, the buffer must contain data
      } else {
        my %handlers = ( 'current' => $fh, 'out' => $out, 'err' => $err, );
        ## Parse data read from buffer
        process_output_text( $text, \@results, \%handlers );
      }
    }
  }
  ## wait the process to end
  waitpid $pid, 0;
  ##############################################################################
  return @results;
}
################################################################################

################################################################################
## process_output_text
##   process data from a buffer and process it according to the buffer it came
##   from.
sub process_output_text {
  my $text     = shift;
  my $results  = shift;
  my $handlers = shift;
  ## get current buffer and both 'parent buffers'.
  my $current = $handlers->{'current'};
  my $out     = $handlers->{'out'};
  my $err     = $handlers->{'err'};
  ## Check if the current buffer is the stdout buffer, if so add data to results
  if ( $current == $out ) {
    push @{$results}, $text;
    ## Check if the current buffer is the stderr buffer, if so, discard data
  } elsif ( $current == $err ) {
    return;
    ## otherwise is an error!
  } else {
    die "Shouldn't be here\n";
  }
  return;
}
################################################################################

################################################################################
## send_to_msa
##   Send aligned sequences to a file
sub send_to_msa {
  my $out_file      = shift;
  my $aln_sequences = shift;
  open my $file_hander, '>', $out_file;
  exit if not print {$file_hander} $aln_sequences;
  close $file_hander;
  return ();
}
################################################################################
