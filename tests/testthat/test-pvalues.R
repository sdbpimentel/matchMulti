
# library(matchMulti)
# library( testthat )

context("calculating pvalues")

data(catholic_schools)

# Trim data to speed up example
catholic_schools = dplyr::filter( catholic_schools, 
                                  female_mean > 0.45,
                                  female_mean < 0.60,
                                  acad > quantile( acad, 0.25 ),
                                  acad < quantile( acad, 0.75 ) )

nrow( catholic_schools )
table(  catholic_schools$sector, catholic_schools$school  )

p.tt <- matchMulti:::ttest.balance( "ses", treatment = "sector", orig.data = catholic_schools )

# Handle clustering
p.CRVE <- matchMulti:::CRVE.balance( "ses", treatment = "sector", school.id = "school", data = catholic_schools )

p.agg <- matchMulti:::agg.balance( "ses", treatment = "sector", school.id = "school", data = catholic_schools )

expect_true( p.CRVE > p.tt )
expect_true( p.agg > p.tt )


