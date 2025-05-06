#!/usr/global/perl-5.8.8/bin/perl

use Getopt::Std;
use Imager;

# Get the argument for the -w switch. Terminate unless the argument was given.

getopt("wso");
die "Error: You must provide the pixel width in nm.\n" unless defined($opt_w);
die "Error: You must provide the name of an STK file.\n" unless defined($opt_s);
$opt_o = "imageout.txt" unless defined($opt_o);

# Extract the TIFF images from the STK file.

my $result;
my $image = Imager->new;
for (my $j = 0; ; $j++) {
    $result = $image->read(file => "$opt_s", page => $j);
    if ($result) {
	    $image->write(file => "$j.tif");
	}
	else {
		last;
	}
}

# Process the TIFF images.

opendir(DIR, ".") or die $!;            # Open the current directory.
my $file;
my $output;
while ($file = readdir(DIR)) {          # Loop through the files in DIR.
    if ($file =~ /.*\.tif/) {           # If the current file is a TIFF file,
    
        # run the beads program on it and collect the output.
    
        open(README, "./beads $opt_w $file |") or die $!;
        while (<README>) {
            $output .= $_;
        }
        close README;
        $output .= "END\n";
    }
}
open(OF, ">$opt_o") or die $!;          # Open the output file.
print OF $output;                       # Dump the results.
close OF;

unlink glob("*.tif");  # Delete the TIFF files.
