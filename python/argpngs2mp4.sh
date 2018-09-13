ffmpeg -r 30 -f image2 -s 640x480 -i outfields/field%04d.png -vcodec libx264 -crf 25 -pix_fmt yuv420p fields_out_$1_$2_$3_$4.mp4
