## code to prepare hex
library(magick)
library(bunny)
w <- 3500
h <- 3500

pkglogo_msk <- image_read_svg("data-raw/cupid_wings.svg", width = w, height = h) %>%
  image_flop()
pkglogo_msk

bg_color <- "#2B2D42" #darker
bg_color2 <- "#F7934C" # lighter
fg_color <- "#9FA2B2" #border
txt_color <- "#F7B05B"
bg_color3 <- "#F4F4F6" #extra light

pkglogo <- image_blank(w,h, color = bg_color2) %>%
  image_composite(pkglogo_msk, operator = "CopyOpacity")

pkgtxt <- image_blank(1200,650, "transparent") %>%  #insert different background to troubleshoot
  image_annotate("qpd", gravity = "center", location = "+0-100",
                 color=txt_color, size=600, font="Aller", weight = 400)

pkgtxt %>% image_scale("20%")

qpd_hex <-
  image_canvas_hex(fill_color = bg_color, border_color = fg_color) %>%
  image_composite(pkglogo, gravity = "center", offset = "-230+200") %>%
  image_composite(pkgtxt, gravity = "center", offset = "+0+0",
                 operator = "Screen") %>%
  image_composite(image_canvas_hexborder(border_color = fg_color, border_size = 10),
                  gravity = "center", operator = "Atop")

qpd_hex %>% image_scale("20%")


qpd_hex %>%
  image_scale("1200x1200") %>%
  image_write("data-raw/qpd_hex.png", density = 600)

qpd_hex %>%
  image_scale("200x200") %>%
  image_write("man/figures/logo.png", density = 600)

qpd_hex_gh <- qpd_hex %>%
  image_scale("400x400")

gh_logo <- bunny::github %>%
  #image_negate() %>%
  image_scale("40x40")

qpd_ghcard <-
  image_canvas_ghcard(fill_color = bg_color3) %>%
  image_composite(qpd_hex_gh, gravity = "East", offset = "+100+0") %>%
  image_annotate("Flexible. Intuitive.", gravity = "West", location = "+100-60",
                 color=bg_color, size=55, font="Aller", weight = 400) %>%
  image_annotate("Quantile.", gravity = "West", location = "+100+10",
                 color=bg_color, size=55, font="Aller", weight = 400) %>%
  image_compose(gh_logo, gravity="West", offset = "+100+90") %>%
  image_annotate("dmi3kno/qpd", gravity="West", color=fg_color,
                 location="+160+90", size=45, font="Ubuntu Mono") %>%
  image_border_ghcard(bg_color3)

qpd_ghcard %>% image_scale("30%")


qpd_ghcard %>%
  image_write("data-raw/qpd_ghcard.png", density = 600)
