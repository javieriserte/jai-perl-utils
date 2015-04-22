#!perl
################################################################################
## This script creates a contact map from a pdb file.
## Uses contactmap.jar program.
## Author: Javier Iserte
## Usage:
##   use JAI::PDB::Contacts
################################################################################

################################################################################
## Module name
package JAI::PDB::Contacts;
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
use IO::Select;
use Symbol;
################################################################################

################################################################################
## Exports
our @EXPORT_OK = qw(new);
our %EXPORT_TAGS = ( DEFAULT => [qw(&new)],
  Both => [qw(&new)] );
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

################################################################################
## Constructor
sub new {
  my $class = shift;
  my $self  = {
    _executable_jar => shift,
    _criteria       => q(),
    _chain_opt      => q(),
  };
  bless $self, $class;
  return $self;
}
################################################################################

################################################################################
## set criteria
sub set_criteria {
  my $self     = shift;
  my $criteria = shift;
  $self->{_criteria} = $criteria;
  return;
}
################################################################################

################################################################################
## show_chains
sub show_chains {
  my $self = shift;
  $self->{_chain_opt} = '--withChains';
  return;
}
################################################################################

################################################################################
## Create contact map
sub calculate_contact_map {
  my $self     = shift;
  my $pdb_file = shift;

  my @results = ();

  ##############################################################################
  ## Create command line to align and sequence data to be aligned
  my $cmdline = $self->{_executable_jar} . q( ) .
    '--infile=' . $pdb_file . q( ) .
    '--criteria=' . $self->{_criteria} . q( ) .
    $self->{_chain_opt};
  ##############################################################################

  ##############################################################################
  ## Launch contact_map process
  my ( $in, $out, $err );
  $err = gensym();
  my $pid = open3( $in, $out, $err, $cmdline );
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
## process_output_text
sub process_output_text {
  my $text_line = shift;
  my $results   = shift;
  my $handlers  = shift;
  ## get current buffer and both 'parent buffers'.
  my $current = $handlers->{'current'};
  my $out     = $handlers->{'out'};
  my $err     = $handlers->{'err'};
  ## Check if the current buffer is the stdout buffer, if so add data to results
  if ( $current == $out ) {
    push @{$results}, $text_line;
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
