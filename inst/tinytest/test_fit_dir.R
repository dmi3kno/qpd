df <- data.frame(id=rep(1:3, each=3),
                spt=rep(0.25*1:3, 3),
                values=c(7, 9, 12, 34, 41, 50, 5, 7.5, 15)/100)
d <- fit_dir(df, "id", "spt", "values")
gd <- fit_gendir(df, "id", "spt", "values")

# checking dirichlet
expect_equal(round(d[["a"]],4), c(3.7679, 12.8603,  2.7001, 10.7198))

#checking generalized dirichlet
expect_equal(round(gd[["a"]],4), c(5.8711, 6.7698, 1.3229))
expect_equal(round(gd[["b"]],4), c(55.2429,  8.0052,  5.2521))
