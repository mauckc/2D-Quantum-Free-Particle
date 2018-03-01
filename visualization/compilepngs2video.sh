ffmpeg -r 30 -f image2 -s 800x800 -i figures/prob2D_%01d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p output-prob2D.mp4
