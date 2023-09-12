# Placeholder with simple test
expect_equal(qnorm(0.4), sqrt(2)*qpd:::qerf(2*0.4-1, lower.tail = TRUE, log.p = FALSE))
expect_equal(qnorm(0.6), sqrt(2)*qpd:::qerf(2*0.4-1, lower.tail = FALSE, log.p = FALSE))

