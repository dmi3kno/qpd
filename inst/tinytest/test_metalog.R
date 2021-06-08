
# Placeholder with simple test
p <- c(0.1, 0.5, 0.85)
q <- c(4,9,17)
a <- fit_metalog(p,q)
al <- fit_metalog(p,q,0)
au <- fit_metalog(p,q, bu=50)
ab <- fit_metalog(p,q, 0, 50)

# test fit_metalog
expect_equal(round(a,4), c(9, 3.5217, 3.1152))
expect_equal(round(al,4), c(2.1972, 0.3678, -0.0032))
expect_equal(round(au,4), c(-3.7136, 0.0912, 0.0970))
expect_equal(round(ab,4), c(-1.5163, 0.4590, 0.0938))

#test qmetalog
expect_equal(round(qmetalog(0.12, a),4), 4.3419)
expect_equal(round(qmetalog(log(0.12), a, log.p=TRUE),4), 4.3419) ## test that log.p works as well
expect_equal(round(qmetalog(0.12, al,0),4), 4.3146)
expect_equal(round(qmetalog(0.12, au, bu=50),4), 4.3144)
expect_equal(round(qmetalog(0.12, ab, 0, 50),4), 4.3146)

#test dqmetalog
expect_equal(round(dqmetalog(0.12, a),4), 0.0628)
expect_equal(round(dqmetalog(log(0.12), a, log.p=TRUE),4), 0.0628)
expect_equal(round(dqmetalog(log(0.12), a, log.p=TRUE, log=TRUE),4), -2.7684)
expect_equal(round(dqmetalog(0.12, a, log=TRUE),4), -2.7684)
expect_equal(round(dqmetalog(0.12, al,0),4), 0.0662)
expect_equal(round(dqmetalog(0.12, au, bu=50),4), 0.0682)
expect_equal(round(dqmetalog(0.12, ab, 0, 50),4), 0.0664)

#test pmetalog
expect_equal(round(pmetalog(7.11, a),4), 0.3501)
expect_equal(round(pmetalog(7.11, a, log.p = TRUE),4), -1.0496)
expect_equal(round(pmetalog(7.11, al,0),4), 0.3452)
expect_equal(round(pmetalog(0.2, al,0, maxiter=1e4, log.p = TRUE),4), -10.3057)
expect_equal(round(pmetalog(7.11, au, bu=50),4), 0.3585)
expect_equal(round(pmetalog(7.11, ab, 0, 50),4), 0.3472)




