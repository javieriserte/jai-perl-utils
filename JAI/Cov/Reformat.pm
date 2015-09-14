#!perl
################################################################################
## This script changes the format of a covariation data file into a common
## standard format
## Author: Javier Iserte
## Usage:
##   use JAI::Cov::Reformat
################################################################################

################################################################################
## Module name
package JAI::Cov::Reformat;
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
use Set::Scalar;
use English '-no_match_vars';
################################################################################

################################################################################
## Exports
our @EXPORT_OK = qw(reformat);
our %EXPORT_TAGS = ( DEFAULT => [qw(&reformat)],
  Both => [qw(&reformat)] );
################################################################################

################################################################################
## Constants
Readonly my $EMPTY              => q();
Readonly my $MAX_SCORE          => 1_000_000;
Readonly my $MIN_SCORE          => -900;
Readonly my $UNDEFINED          => -999;
Readonly my $MI_P1_INDEX        => 1;
Readonly my $MI_AA1_INDEX       => 2;
Readonly my $MI_P2_INDEX        => 4;
Readonly my $MI_AA2_INDEX       => 5;
Readonly my $MI_SCORE_INDEX     => 11;
Readonly my $DCA_P1_INDEX       => 0;
Readonly my $DCA_AA1_INDEX      => 1;
Readonly my $DCA_P2_INDEX       => 2;
Readonly my $DCA_AA2_INDEX      => 3;
Readonly my $DCA_SCORE_INDEX    => 5;
Readonly my $P1_COLUMN_INDEX            => 0;
Readonly my $P2_COLUMN_INDEX            => 1;
Readonly my $SCORE_COLUMN_INDEX         => 2;
Readonly my $PSICOV_P1_INDEX    => 0;
Readonly my $PSICOV_P2_INDEX    => 1;
Readonly my $PSICOV_SCORE_INDEX => 4;
Readonly my $MIN_FIELDS_SIZE_FOR_MATRIX => 15;
Readonly my $FLOAT_REGEX        => '[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?';
Readonly my $THREE_COLUMNS_REGEX => '^[0-9]+\s[0-9]+\s' . $FLOAT_REGEX . q($);
Readonly my $PSICOV_REGEX       => '^[0-9]+\s[0-9]+\s[0-9]+\s[0-9]+\s' .
  $FLOAT_REGEX . q($);
Readonly my $DCA_REGEX => '^[0-9]+\s[A-Za-z]\s[0-9]+\s[A-Za-z]\s' .
  $FLOAT_REGEX . '\s' . $FLOAT_REGEX . q($);
Readonly my $MI_REGEX => '^MI\[\s[0-9]+\s[A-Za-z]\s\]\[\s[0-9]+\s' .
  '[A-Za-z]\s\]\s=\s' . $FLOAT_REGEX . '\s' .
  $FLOAT_REGEX . '\s' . $FLOAT_REGEX . '\s' .
  $FLOAT_REGEX . q($);
################################################################################

################################################################################
## Version
our $VERSION = 1.0000;
################################################################################

################################################################################
## reformat
sub reformat {
  my $infile_path  = shift;
  my $out_path     = shift;
  my $first_length = shift;
  my $in_msa       = shift;

  ##############################################################################
  ## Check if the input data is matrix formatted
  if ( check_if_matrix($infile_path) ) {
    $infile_path = matrix_to_paired($infile_path);
  }
  ##############################################################################


  ##############################################################################
  ## Tries to get a function that reads text lines in the format of the input
  ## get_format_reader, returns a function reference.
  my $reader = get_format_reader($infile_path);
  ##############################################################################

  if ($reader) {
    ############################################################################
    ## Get label map and the total length from input msa.
    my ( $input_length, $ref_seq, $label_map ) =
      get_descriptions_from_msa( $in_msa, $first_length );
    ############################################################################

    ############################################################################
    ## Get the descriptions of the covariation score file that will be in the
    ## header of the output file
    my ( $min, $max, $npair ) = get_descriptions( $infile_path, $reader );
    ############################################################################

    ############################################################################
    ## Convert format
    my $header = get_header( $npair, $input_length, $first_length, $min, $max,
      $ref_seq, $label_map );
    export_data( $header, $infile_path, $out_path, $reader, $input_length );
    ############################################################################
  } else {
    exit 1;
  }
  exit 0;
}

################################################################################
###                                Functions                                ####
################################################################################

################################################################################
##  get_descriptions_from_msa
sub get_descriptions_from_msa {
  my $in_msa_file  = shift;
  my $first_length = shift;
  my $seq_in       = Bio::SeqIO->new(
    -file   => '<' . $in_msa_file,
    -format => 'fasta', );
  if ( my $first_seq = $seq_in->next_seq() ) {
    my $ref_seq        = $first_seq->seq();
    my %labels         = ();
    my $ref_seq_length = length $ref_seq;
    for my $i ( 1 .. length $ref_seq ) {
      $labels{$i} = ( $i <= $first_length ? $i : ( $i - $first_length ) );
    }
    return ( $ref_seq_length, $ref_seq, join q(,),
      map { $labels{$ARG} } sort { $a <=> $b } keys %labels );
  }
  return ( 0, q(), q() );
}
################################################################################

################################################################################
##  get_header
##     Returns a the header for the output file.
sub export_data {
  my $header       = shift;
  my $infile       = shift;
  my $outfile      = shift;
  my $reader       = shift;
  my $input_length = shift;
  open my $file_handler_in,  '<', $infile;
  open my $file_handler_out, '>', $outfile;
  exit if !print {$file_handler_out} $header;
  export_data_aux( $file_handler_in, $file_handler_out, $reader, $input_length );
  close $file_handler_in;
  close $file_handler_out;
  return;
}

sub export_data_aux {
  my $file_handler_in  = shift;
  my $file_handler_out = shift;
  my $reader           = shift;
  my $input_length     = shift;
  my %read_pairs       = ();
  while ( my $line = <$file_handler_in> ) {
    chomp $line;
    if ( $line =~ m/^[^#]/msx ) {
      my ( $pos1, $aa1, $pos2, $aa2, $score ) = $reader->($line);
      my $key = $pos1 . q(_) . $pos2;
      $read_pairs{$key} = $score;
    }
  }
  for my $i ( 1 .. $input_length - 1 ) {
    for my $j ( $i + 1 .. $input_length ) {
      my $key = $i . q(_) . $j;
      my $score = defined $read_pairs{$key} ? $read_pairs{$key} : $UNDEFINED;
      exit if !print {$file_handler_out}
        $i . q( ) . $j . q( ) . $score . "\n";
    }
  }
  return;
}
################################################################################

################################################################################
##  get_header
##     Returns a the header for the output file.
sub get_header {
  my @ids = qw(number_of_pairs total_length first_protein_length
    min_value max_value reference_sequence label_map);
  my $header = $EMPTY;
  my $lower_index = scalar @ARG > scalar @ids ? scalar @ids : scalar @ARG;
  for my $i ( 0 .. $lower_index - 1 ) {
    $header .= q(#) . q( ) . $ids[$i] . q( ) . $ARG[$i] . "\n";
  }
  return $header;
}
################################################################################

################################################################################
##  get_descriptions
##     Reads the covariation score file and gets description data
sub get_descriptions {
  my $infile = shift;
  my $reader = shift;
  my $min    = $MAX_SCORE;
  my $max    = $MIN_SCORE;
  my $npair  = 0;
  open my $file_handler, '<', $infile;

  while ( my $line = <$file_handler> ) {
    chomp $line;
    if ( $line =~ m/^[^#]/msx ) {
      get_descriptions_aux( $line, $reader, \$max, \$min, \$npair );
    }
  }
  close $file_handler;
  return ( $min, $max, $npair );
}

sub get_descriptions_aux {
  my ( $line, $reader, $max_ref, $min_ref, $npair_ref ) = @ARG;
  my ( $pos1, $aa1,    $pos2,    $aa2,     $score )     = $reader->($line);
  if ( $score > ${$max_ref} ) { ${$max_ref} = $score; }
  if ( $score < ${$min_ref} && $score > $MIN_SCORE ) { ${$min_ref} = $score; }
  ${$npair_ref} = ${$npair_ref} + 1;
  return;
}
################################################################################

################################################################################
## Function to read covariation data from morten's mi soft
sub mi_parser {
  my $current_line = shift;
  my @data = split /\s/msx, $current_line;
  return ( $data[$MI_P1_INDEX], $data[$MI_AA1_INDEX], $data[$MI_P2_INDEX],
    $data[$MI_AA2_INDEX], $data[$MI_SCORE_INDEX] );
}
################################################################################

################################################################################
## Function to read covariation data from evfold format files
sub dca_parser {
  my $current_line = shift;
  my @data = split /\s/msx, $current_line;
  return ( $data[$DCA_P1_INDEX], $data[$DCA_AA1_INDEX], $data[$DCA_P2_INDEX],
    $data[$DCA_AA2_INDEX], $data[$DCA_SCORE_INDEX] );
}
################################################################################

################################################################################
## Function to read covariation data from psicov format files
sub psicov_parser {
  my $current_line = shift;
  my @data = split /\s/msx, $current_line;
  return ( $data[$PSICOV_P1_INDEX], q(-), $data[$PSICOV_P2_INDEX],
    q(-), $data[$PSICOV_SCORE_INDEX] );
}
################################################################################

################################################################################
## Function to read covariation data from data in three columns format
sub three_columns_parser {
  my $current_line = shift;
  my @data = split /\s/msx, $current_line;
  return ( $data[$P1_COLUMN_INDEX], q(-), $data[$P2_COLUMN_INDEX],
    q(-), $data[$SCORE_COLUMN_INDEX] );
}
################################################################################


################################################################################
##  get_format_reader
##     Reads the covariation score file and gets a function able to read the
##     data correctly.
sub get_format_reader {
  my $infile = shift;
  open my $file_handler, '<', $infile;
  my $line = <$file_handler>;
  close $file_handler;
  if ($line) {
    chomp $line;
    if ( $line =~ m/$DCA_REGEX/msx ) {
      return \&dca_parser;
    }
    if ( $line =~ m/$MI_REGEX/msx ) {
      return \&mi_parser;
    }
    if ( $line =~ m/$PSICOV_REGEX/msx ) {
      return \&psicov_parser;
    }
    if ( $line =~ m/$THREE_COLUMNS_REGEX/msx ) {
      return \&three_columns_parser;
    }
  }
  return 0;
}
################################################################################

################################################################################
## check_if_matrix
##   reads a covariation file and checks if the content is a matrix table
sub check_if_matrix {
  my $file = shift;
  open my $file_handler, '<', $file;
  my $line = <$file_handler>;
  close $file_handler;
  chomp $line;
  if ( scalar( split /\s/msx, $line ) >= $MIN_FIELDS_SIZE_FOR_MATRIX ) {
    return 1;
  } else {
    return 0;
  }
}
################################################################################

################################################################################
## Convert a matrix-like format of ccmpred to a three column format like any
## other method
sub matrix_to_paired {
  my $file        = shift;
  my $row_counter = 0;
  my %data        = ();
  open my $file_handler, '<', $file;
  while ( my $line = <$file_handler> ) {
    $row_counter++;
    chomp $line;
    my @columns = split /\s/msx, $line;
    $data{$row_counter} = \@columns;
  }
  close $file_handler;
  my $new_file = $file . '.tmp';
  open my $file_handler_2, '>', $new_file;
  matrix_to_paired_aux( $file_handler_2, $row_counter, \%data );
  close $file_handler_2;
  return $new_file;
}

sub matrix_to_paired_aux {
  my $file_handler = shift;
  my $row_counter  = shift;
  my $data         = shift;
  for my $i ( 1 .. $row_counter ) {
    for my $j ( $i + 1 .. $row_counter ) {
      exit if not print {$file_handler}
        $i . "\t" . $j . "\t" . $data->{$i}->[ $j - 1 ] . "\n";
    }
  }
  return;
}
################################################################################
