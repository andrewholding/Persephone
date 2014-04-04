use strict;

package Persephone::Config;

use base 'Exporter';
our @EXPORT = ('import_experiment_definition_file', 'load_constants', 'no_of_threads');


sub no_of_threads {
  return 12;
}


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


sub load_constants {

my %constants;

$constants{'mass_of_deuterium'} = 2.01410178;
$constants{'mass_of_hydrogen'}  = 1.00783;
$constants{'mass_of_proton'}    = 1.00728;
$constants{'mass_of_carbon13'}  = 13.00335;
$constants{'mass_of_carbon12'}  = 12;

return \%constants;
}

1;