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


   my %ms2_fragmentation;

  my $scored_ions = '';
  if ( defined $setting{'score_aions'} ) {
    $scored_ions = $scored_ions . 'A-ions, ';
    $ms2_fragmentation{'aions-score'} = '1';
  } else {
        $ms2_fragmentation{'aions-score'} = '0';
  }
  if ( defined $setting{'score_bions'} ) {
    $ms2_fragmentation{'bions-score'} = '1';
    $scored_ions = $scored_ions . 'B-ions, ';
  } else {
        $ms2_fragmentation{'bions-score'} = '0';
  }
  if ( defined $setting{'score_yions'} ) {
    $ms2_fragmentation{'yions-score'} = '1';
    $scored_ions = $scored_ions . 'Y-ions, ';
  } else {
        $ms2_fragmentation{'yions-score'} = '0';
  }
  if ( defined $setting{'score_waterions'} ) {
    $ms2_fragmentation{'waterloss-score'} = '1';
    $scored_ions = $scored_ions . 'Water-loss ions, ';
  } else {
        $ms2_fragmentation{'waterloss-score'} = '0';
  }
  if ( defined $setting{'score_ammoniaions'} ) {
    $ms2_fragmentation{'ammonialoss-score'} = '1';
    $scored_ions = $scored_ions . 'ammonia-loss ions, ';
  }else {
        $ms2_fragmentation{'ammonialoss-score'} = '0';
    }
  
    if   (defined $setting{'calc_aions'}) { $ms2_fragmentation{'aions'} = '1' }
    else                                  { $ms2_fragmentation{'aions'} = '0' }
    if   (defined $setting{'calc_bions'}) { $ms2_fragmentation{'bions'} = '1' }
    else                                  { $ms2_fragmentation{'bions'} = '0' }
    if   (defined $setting{'calc_yions'}) { $ms2_fragmentation{'yions'} = '1' }
    else                                  { $ms2_fragmentation{'yions'} = '0' }
    if   (defined $setting{'calc_waterloss'}) { $ms2_fragmentation{'waterloss'} = '1'; }
    else                                      { $ms2_fragmentation{'waterloss'} = '0'; }
    if   (defined $setting{'calc_ammonialoss'}) { $ms2_fragmentation{'ammonialoss'} = '1'; }
    else                                        { $ms2_fragmentation{'ammonialoss'} = '0'; }

   $setting{'scored_ions'} = $scored_ions;
   $setting{'ms2_fragmentation'} = \%ms2_fragmentation;

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