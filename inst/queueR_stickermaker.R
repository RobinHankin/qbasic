## This file makes the hex sticker in file qbasic.png, which lives in
## the man/figures directory of the package.  The letter 'q' is in the
## Luminari font, with a modified tick that matches the nearest edge
## of the hexagon.

library("hexSticker")

sticker("queueR_modified_tick.png", package="queueR", p_size=30, s_x=1, s_y=0.91,
        s_width=0.6,asp=sqrt(3)/2, white_around_sticker=TRUE, h_fill="#7733FF",
        h_color="#000000", filename="queueR.png")
