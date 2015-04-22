#!perl
################################################################################
## This script ....
## Author: ....
## Usage:
##   use JAI::Cov::AlignToPDB
################################################################################

################################################################################
## Module name
package JAI::Cov::AlignToPDB;
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
use JAI::Cov::XmlRepr;
use JAI::Fasta::Io qw(read_first id_map_from_fasta_text);
use JAI::Lists::Reduced qw(concat);
use JAI::PDB::Seq qw(get_chain_res);
use IO::Select;
use Symbol;
################################################################################

################################################################################
## Exports
our @EXPORT_OK = qw(align_to_pdb);
our %EXPORT_TAGS = ( DEFAULT => [qw(&align_to_pdb)],
  Both => [qw(&align_to_pdb)] );
################################################################################

################################################################################
## Constants
Readonly my $FALSE       => 0;
Readonly my $TRUE        => 1;
Readonly my $BUFFER_SIZE => 4096;
################################################################################

################################################################################
## Version
our $VERSION = 1.0000;
################################################################################

sub align_to_pdb {
  my $repr = shift;
  ##############################################################################
  ## extract msa and pdb files from repr
  my $tmp_msa_file = 'tmp_msa.fas';
  my $tmp_pdb_file = 'tmp_pdb.pdb';
  JAI::Cov::XmlRepr::tag_data_to_file( $repr, 'msa', $TRUE, $tmp_msa_file );
  JAI::Cov::XmlRepr::tag_data_to_file( $repr, 'pdb', $TRUE, $tmp_pdb_file );
  ##############################################################################

  ##############################################################################
  ## Read input fasta sequence,
  my $len_first  = JAI::Cov::XmlRepr::get_tag_data( $repr, 'first_length' );
  my $len_second = JAI::Cov::XmlRepr::get_tag_data( $repr, 'second_length' );
  ## An (id, sequence) tuple is returned from function
  my @fasta_seq  = read_first($tmp_msa_file);
  my $first_seq  = substr $fasta_seq[1], 0, $len_first;
  my $second_seq = substr $fasta_seq[1], $len_first, $len_second;
  ##############################################################################

  ##############################################################################
  ## Get sequence from pdb file
  my $first_chain  = JAI::Cov::XmlRepr::get_tag_data( $repr, 'first_chain' );
  my $second_chain = JAI::Cov::XmlRepr::get_tag_data( $repr, 'second_chain' );
  my $pdb_seq_first  = get_chain_res( $tmp_pdb_file, $first_chain,  'X' );
  my $pdb_seq_second = get_chain_res( $tmp_pdb_file, $second_chain, 'X' );
  ##############################################################################

  ##############################################################################
  ## align sequences
  my @aln_sequences_first  = align( $first_seq, $pdb_seq_first );
  my $seq_string_first     = concat(@aln_sequences_first);
  my @aln_sequences_second = align( $second_seq, $pdb_seq_second );
  my $seq_string_second    = concat(@aln_sequences_second);
  ##############################################################################

  ##############################################################################
  ## export msa
  my $tmp_pdb_align_first  = 'pdb_align_first.fas';
  my $tmp_pdb_align_second = 'pdb_align_second.fas';
  send_to_msa( $tmp_pdb_align_first,  $seq_string_first );
  send_to_msa( $tmp_pdb_align_second, $seq_string_second );
  JAI::Cov::XmlRepr::tag_data_from_file(
    $repr, 'pdb_align_first', $TRUE, $tmp_pdb_align_first );
  JAI::Cov::XmlRepr::tag_data_from_file(
    $repr, 'pdb_align_second', $TRUE, $tmp_pdb_align_second );
  ##############################################################################

  ##############################################################################
  ## export seq to pdb position map
  my $tmp_map_file = 'tmp_map.map';
  my $current_map  = create_map($seq_string_first);
  $current_map = create_map( $seq_string_second, $len_first, $current_map );
  open my $file_hander, '>', $tmp_map_file;
  foreach my $fasta_pos ( sort { $a <=> $b } keys %{$current_map} ) {
    exit if not print {$file_hander} $fasta_pos . "\t" .
      $current_map->{$fasta_pos} . "\n";
  }
  close $file_hander;
  JAI::Cov::XmlRepr::tag_data_from_file(
    $repr, 'msa_to_pdb_map', $TRUE, $tmp_map_file );
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
## export_map
sub create_map {
  my $seq_string  = shift;
  my $num_offset  = shift;
  my $current_map = shift;
  my %seqs        = id_map_from_fasta_text($seq_string);
  my ( $seq_i, $seq_j ) = ( 0, 0 );
  $num_offset = defined $num_offset ? $num_offset : 0;
  if ( not defined $current_map ) {
    my %map = ();
    $current_map = \%map;
  }
  for my $i ( 0 .. ( length $seqs{'pdb'} ) - 1 ) {
    my $pdb_is_not_gap   = ( substr $seqs{'pdb'},   $i, 1 ) ne q(-);
    my $pdb_is_not_x     = ( substr $seqs{'pdb'},   $i, 1 ) ne q(X);
    my $fasta_is_not_gap = ( substr $seqs{'fasta'}, $i, 1 ) ne q(-);
    $seq_i = $seq_i + ( $pdb_is_not_gap   ? 1 : 0 );
    $seq_j = $seq_j + ( $fasta_is_not_gap ? 1 : 0 );
    if ( $pdb_is_not_gap && $fasta_is_not_gap && $pdb_is_not_x ) {
      $current_map->{ $seq_j + $num_offset } = $seq_i;
    }
  }
  return $current_map;
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
