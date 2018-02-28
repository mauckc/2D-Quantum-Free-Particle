DEPENDENCIES:

  Install ffmpeg for converting images to video
______________
Linux (Ubuntu)

sudo apt-get install ffmpeg

_______________________
Mac OS X with homebrew:

  To install FFmpeg with all available options without disabling anything execute:

brew install ffmpeg --with-fdk-aac --with-ffplay --with-freetype --with-frei0r --with-libass --with-libvo-aacenc --with-libvorbis --with-libvpx --with-opencore-amr --with-openjpeg --with-opus --with-rtmpdump --with-schroedinger --with-speex --with-theora --with-tools



CREATE VISUALIZATIONS:
______
STEP 1 PICS: Convert 2Dparticle output (.dat) to images using python2 and matplotlib

python convert-data-to-video.py
______
STEP 2 VIDEO: Convert plot images to video using FFmpeg library

./compilepngs2video.sh

  debug:
    (ffmpeg shell script may need permissions.  If you run into problems when running the command above.
    To give permissions to the shell script for running within quantum-visualize-2D-scalar-field we must give it permissions):

    chmod u+rx compilepngs2video.sh
