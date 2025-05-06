use Getopt::Std;
use Imager;

# Get the argument for the -s switch. Terminate unless the argument was given.

getopt("s");
die "Error: You must provide the name of an STK file.\n" unless defined($opt_s);

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

