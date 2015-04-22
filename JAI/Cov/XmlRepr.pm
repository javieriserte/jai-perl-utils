#!perl
################################################################################
## This script contains a data structure that represents covariation data
## calculation pipeline and functions to alter it.
## Author: Javier Iserte
## Usage:
##   use JAI::Cov::XmlRepr
################################################################################

################################################################################
## Module name
package JAI::Cov::XmlRepr;
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
use MIME::Base64 qw(encode_base64 decode_base64);
use Compress::Zlib qw(compress uncompress);
use XML::Simple;
use File::Slurp;
use English qw( -no_match_vars );
use List::Util qw(reduce);
################################################################################

################################################################################
## Exports
our @EXPORT_OK = qw(new add_tag repr_to_xml_file xml_file_to_repr get_tag_data);
our %EXPORT_TAGS = ( DEFAULT => [
    qw(&new &add_tag &repr_to_xml_file &xml_file_to_repr &get_tag_data) ],
  Both => [qw(&new &add_tag &repr_to_xml_file &xml_file_to_repr &get_tag_data)] );
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
## new
sub new {
  my %new_hash = ();
  my %repr = ( 'root' => 'cov_data',
    'data' => \%new_hash );
  return \%repr;
}
################################################################################

################################################################################
## contains_tag
sub contains_tag {
  my $repr     = shift;
  my $tag_name = shift;
  return exists $repr->{'data'}{$tag_name};
}
################################################################################

################################################################################
## remove_tag
sub remove_tag {
  my $repr     = shift;
  my $tag_name = shift;
  return delete $repr->{'data'}{$tag_name};
}
################################################################################

################################################################################
## add data tag
sub add_tag {
  my $xml_repr    = shift;
  my $tag_name    = shift;
  my $tag_data    = shift;
  my $encode_flag = shift;
  if ( defined $encode_flag && $encode_flag ) {
    $tag_data = compress_and_encode($tag_data);
  }
  $xml_repr->{'data'}{$tag_name} = $tag_data;
  return;
}
################################################################################

################################################################################
## tag_data_from_file
sub tag_data_from_file {
  my $xml_repr    = shift;
  my $tag_name    = shift;
  my $encode_flag = shift;
  my $infile      = shift;
  my $tag_data    = read_file($infile);
  add_tag( $xml_repr, $tag_name, $tag_data, $encode_flag );
  return;
}
################################################################################

################################################################################
## tag_data_to_file
sub tag_data_to_file {
  my $xml_repr    = shift;
  my $tag_name    = shift;
  my $encode_flag = shift;
  my $outfile     = shift;
  open my $file_handler, '>', $outfile;
  exit if not print {$file_handler}
    get_tag_data( $xml_repr, $tag_name, $encode_flag );
  close $file_handler;
  return;
}
################################################################################

################################################################################
## get data from tag
sub get_tag_data {
  my $xml_repr    = shift;
  my $tag_name    = shift;
  my $encode_flag = shift;
  my $tag_data    = $xml_repr->{'data'}{$tag_name};
  if ( defined $encode_flag && $encode_flag ) {
    $tag_data = decode_and_uncompress($tag_data);
  }
  return $tag_data;
}
################################################################################

################################################################################
## repr_to_xml_file
sub repr_to_xml_file {
  my $xml_repr  = shift;
  my $file_name = shift;
  my $xml = XMLout( $xml_repr->{'data'}, NoAttr => 1, RootName => $xml_repr->{'root'} );
  open my $file_handler, '>', $file_name;
  exit if not print {$file_handler} $xml;
  close $file_handler;
  return;
}
################################################################################

################################################################################
## xml_file_to_repr
sub xml_file_to_repr {
  my $file_name = shift;
  my $xml       = XMLin($file_name);
  my $repr      = new();
  $repr->{'data'} = $xml;
  return $repr;
}
################################################################################

################################################################################
## compress_and_encode
sub compress_and_encode {
  my $data = shift;
  return encode_base64( compress($data) );
}
################################################################################

################################################################################
## decode and uncompress
sub decode_and_uncompress {
  my $data = shift;
  return uncompress( decode_base64($data) );
}
################################################################################

################################################################################
## next_cov_data_index
##  searches a covariation Repr data for all tags that correspond to
##  calculation data. Get the index value of each one and returns the
##  following of the larger index.
sub next_cov_data_index {
  my $repr = shift;
  my @indexes = map { $ARG => 1 }
    map { ( split /_/msx, $ARG )[1] }
    grep { m/^calc/msx } keys $repr->{'data'};
  my $last_index = reduce { $a > $b ? $a : $b } @indexes;
  return defined $last_index ? $last_index + 1 : 1;
}
################################################################################
