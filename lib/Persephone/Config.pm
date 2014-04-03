use strict;

package Persephone::Config;

use base 'Exporter';
our @EXPORT = ('import_experiment_definition_file');


sub import_experiment_definition_file {
#Imports the settings definted in the provided file and returns as hash reference
my ($experiment_definition_file) = @_;
my %setting;

if (defined $experiment_definition_file)
{
  print "Experiment definition file = ", $experiment_definition_file, "\n";
} else {
  die "Error: No experiment definition file provided\n";
}

open RC_FILE, "<", $experiment_definition_file or die "Error: Unable to open experiment definition file";
while (my $line = <RC_FILE>) {
  chomp $line;
  if ($line =~ /=/) {
    my @settings = split '=', $line;
    $setting{$settings[0]} = $settings[1]; 
  }
}
close RC_FILE;

return \%setting;
}