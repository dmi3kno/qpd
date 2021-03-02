## code to prepare `make_logo` dataset goes here
library(magick)
library(bunny)

bg_color <- "white"
fg_color <- "firebrick3"

qpr_hex <-
	image_canvas_hex(fill_color = fg_color, border_color = fg_color) %>% 
	image_composite(image_blank(500, 2000, "white"), operator = "atop", gravity = "north") %>% 
	image_composite(image_blank(50, 2000, "blue"), operator = "atop", gravity = "north", offset = "+200+0") %>% 
	image_composite(image_blank(50, 2000, "green"), operator = "atop", gravity = "north", offset = "-100+0")
qpr_hex %>% image_scale("50%")

clrmd_hex %>%
	image_scale("1200x1200") %>%
	image_write("data-raw/clrmd_hex.png", density = 600)

clrmd_hex %>%
	image_scale("200x200") %>%
	image_write("man/figures/logo.png", density = 600)

clrmd_hex_gh <- clrmd_hex %>%
	image_scale("480x480")

gh_logo <- bunny::github %>%
	image_negate() %>%
	image_scale("50x50")

clrmd_ghcard <-
	image_canvas_ghcard(fill_color = bg_color2) %>%
	image_composite(clrmd_hex_gh, gravity = "East", offset = "+100+0") %>%
	image_annotate("What color are you?", gravity = "West", location = "+60-30",
																color="white", size=55, font="Aller", weight = 400) %>%
	image_compose(gh_logo, gravity="West", offset = "+60+40") %>%
	image_annotate("dmi3kno/colormind", gravity="West", color=fg_color,
																location="+120+45", size=50, font="Ubuntu Mono") %>%
	image_border_ghcard(bg_color2)

clrmd_ghcard %>%
	image_write("data-raw/clrmd_ghcard.png", density = 600)
