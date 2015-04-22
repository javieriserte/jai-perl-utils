#!perl
################################################################################
## This script is used to create simple random covarying msas.
## Author: Javier Iserte
## Usage:
##   use JAI::Cov::GenRndMsa
################################################################################

################################################################################
## Module name
package JAI::Cov::GenRndMsa;
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
################################################################################

################################################################################
## Exports
our @EXPORT_OK = qw(new);
our %EXPORT_TAGS = ( DEFAULT => [qw(&new)],
  Both => [qw(&new)] );
################################################################################

################################################################################
## Constants
Readonly my $AMINO_ACIDS                 => 'QWERTYIPASDFGHKLCVNM';
Readonly my $NUMBER_OF_AMINO_ACIDS       => 20;
Readonly my $DEF_NUMBER_OF_SEQUENCES     => 1000;
Readonly my $DEF_MSA_WIDDTH              => 100;
Readonly my $DEF_COVARIATION_PROBABILITY => 0.3;
Readonly my $DEF_MUTATION_PROBABILITY    => 0.2;
Readonly my $DEF_NUMBER_OF_CONTACTS      => 100;
################################################################################

################################################################################
## Version
our $VERSION = 1.0000;
################################################################################

################################################################################
## Constructor
sub new {
  my $class = shift;
  my $self  = {
    _cov_prob    => $DEF_COVARIATION_PROBABILITY,
    _mut_prob    => $DEF_MUTATION_PROBABILITY,
    _nseq        => $DEF_NUMBER_OF_SEQUENCES,
    _width       => $DEF_MSA_WIDDTH,
    _initial_seq => q(),
    _custom_map  => q(),
    _n_contacts  => $DEF_NUMBER_OF_CONTACTS,
    _msa_out     => q(),
    _map_out     => q(),
  };
  bless $self, $class;
  return $self;
}
################################################################################

################################################################################
###                                Methods                                  ####
################################################################################

################################################################################
## gen_random_msa
sub gen_random_msa {
  my $self = shift;
  my @map;
  if ( $self->{_custom_map} eq q() ) {
    @map = $self->create_concat_map();
  } else {
    @map = $self->load_contact_map();
  }
  if ( $self->{_initial_seq} eq q() ) {
    $self->{_initial_seq} = $self->gen_random_seq( $self->{_width} );
  }
  my $generated_seq = $self->new_mutated_seqs( \@map );
  $self->export_msa( $generated_seq, $self->{_msa_out} );
  $self->export_contacts( \@map, $self->{_map_out} );
  return;
}
################################################################################

################################################################################
## new_mutated_seqs
##   Generate new mutated sequences
sub new_mutated_seqs {
  my $self          = shift;
  my $contact_map   = shift;
  my @generated_seq = ();
  for my $index ( 1 .. $self->{_nseq} ) {
    my $new_seq = $self->{_initial_seq};
    $new_seq = $self->mutate($new_seq);
    $new_seq = $self->cov( $new_seq, $contact_map );
    push @generated_seq, $new_seq;
  }
  return \@generated_seq;
}
################################################################################

################################################################################
## load_contact_map
sub load_contact_map {
  my $self      = shift;
  my @all_pairs = ();
  open my $file_handler, '<', $self->{_custom_map};
  while ( my $line = <$file_handler> ) {
    chomp $line;
    my @data = split /\s/msx, $line;
    push @all_pairs, [ $data[0], $data[1] ];
  }
  close $file_handler;
  return @all_pairs;
}
################################################################################

################################################################################
## create_concat_map
sub create_concat_map {
  my $self      = shift;
  my @map       = ();
  my @all_pairs = ();
  ##############################################################################
  ## Initialize contact map
  for my $i ( 0 .. $self->{_width} - 1 ) {
    my @row = ();
    for my $j ( 0 .. $self->{_width} - 1 ) {
      $row[$j] = 0;
    }
    $map[$i] = [@row];
  }
  ##############################################################################

  ##############################################################################
  ## Get all pairs
  for my $i ( 0 .. $self->{_width} - 1 ) {
    for my $j ( $i + 1 .. $self->{_width} - 1 ) {
      push @all_pairs, [ $i, $j ];
    }
  }
  ##############################################################################

  ##############################################################################
  ## Shuffle all pairs and pick some contacts
  for my $i ( 0 .. scalar @all_pairs - 1 ) {
    my $j   = int( rand( scalar @all_pairs - $i ) ) + $i;
    my $tmp = $all_pairs[$j];
    $all_pairs[$j] = $all_pairs[$i];
    $all_pairs[$i] = $tmp;
  }
  splice @all_pairs, $self->{_n_contacts};
  ##############################################################################
  return @all_pairs;
}
################################################################################

################################################################################
## cov
sub cov {
  my $self  = shift;
  my $seq   = shift;
  my $map_r = shift;
  my @map   = @{$map_r};
  for my $m (@map) {
    my $r = rand;
    if ( $r < $self->{_cov_prob} ) {
      substr $seq, $m->[0], 1, $self->get_random_aa();
      substr $seq, $m->[1], 1, $self->get_random_aa();
    }
  }
  return $seq;
}
################################################################################

################################################################################
## mutate
sub mutate {
  my $self = shift;
  my $seq  = shift;
  for my $i ( 0 .. ( length $seq ) - 1 ) {
    my $r = rand;
    if ( $r < $self->{_mut_prob} ) {
      substr $seq, $i, 1, $self->get_random_aa();
    }
  }
  return $seq;
}
################################################################################

################################################################################
## gen_random_seq
sub gen_random_seq {
  my $self = shift;
  my $seq = pack 'A' . $self->{_width}, q();
  for my $index ( 0 .. ( length $seq ) - 1 ) {
    substr $seq, $index, 1, $self->get_random_aa();
  }
  return $seq;
}
################################################################################

################################################################################
## get_random_aa
sub get_random_aa {
  my $pos = int rand $NUMBER_OF_AMINO_ACIDS;
  return substr $AMINO_ACIDS, $pos, 1;
}
################################################################################

################################################################################
## Export sequences
sub export_msa {
  my $self          = shift;
  my $generated_seq = shift;
  open my $file_handler, '>', $self->{_msa_out};
  my $counter = 1;
  for my $seq ( @{$generated_seq} ) {
    exit if !print {$file_handler} '>seq_' . $counter . "\n" . $seq . "\n";
    $counter++;
  }
  return close $file_handler;
}
################################################################################

################################################################################
## export_contacts
sub export_contacts {
  my $self = shift;
  my $map  = shift;
  if ( $self->{_map_out} ne q() ) {
    open my $file_handler, '>', $self->{_map_out};
    for my $pair ( @{$map} ) {
      exit if !print {$file_handler} $pair->[0] . q( ) . $pair->[1] . "\n";
    }
    return close $file_handler;
  }
}
################################################################################

################################################################################
## set_cov_prob
sub set_cov_prob {
  my $self = shift;
  my $new_value= shift;
  $self->{_cov_prob} = $new_value;
}

sub set_mut_prob {
  my $self = shift;
  my $new_value= shift;
  $self->{_mut_prob} = $new_value;
}

sub set_nseq {
  my $self = shift;
  my $new_value= shift;
  $self->{_nseq} = $new_value;
}

sub set_width {
  my $self = shift;
  my $new_value= shift;
  $self->{_width} = $new_value;
}

sub set_initial_seq {
  my $self = shift;
  my $new_value= shift;
  $self->{_initial_seq} = $new_value;
}

sub use_random_initial_seq {
  my $self = shift;
  my $new_value = shift;
  $self->{_initial_seq} = $new_value;
}

sub set_custom_map {
  my $self = shift;
  my $new_value= shift;
  $self->{_custom_map} = $new_value;
}

sub set_n_contacts {
  my $self = shift;
  my $new_value= shift;
  $self->{_n_contacts} = $new_value;
}

sub set_msa_out {
  my $self = shift;
  my $new_value= shift;
  $self->{_msa_out} = $new_value;
}

sub set_map_out {
  my $self = shift;
  my $new_value= shift;
  $self->{_map_out} = $new_value;
}

sub use_random_map {
  my $self = shift;
  $self->{_map_out} = q();
}
