
# library(matchMulti)
# library( testthat )

context("balance tables")

data(catholic_schools)

# Trim data to speed up example
catholic_schools <- catholic_schools[catholic_schools$female_mean >.45 &
 catholic_schools$female_mean < .60,]
nrow( catholic_schools )
head( catholic_schools )

student.cov <- c('mathach', 'minority', 'female')

# Does school.fb work?
catholic_schools = dplyr::mutate( catholic_schools,
                                  size_cut = cut( size, 2 ),
                                  discrm_cut = cut( discrm, 2 ) )

catholic_schools = dplyr::filter( catholic_schools, 
                                  acad > quantile( acad, 0.25 ) & acad < quantile( acad, 0.75 ) )

nrow( catholic_schools )
table(  catholic_schools$sector, catholic_schools$school  )

school.fb <- c( "discrm_cut", "size_cut" )
school.cov = c( "discrm", "size" )

match.simpleA <- matchMulti(catholic_schools, treatment = 'sector',
                            school.id = 'school', match.students = FALSE,
                            school.fb = school.fb )
match.simpleA


test_that( "balanceMulti works", {

btab_split = balanceMulti( match.simpleA,
                           school.cov = school.cov, 
                           student.cov = student.cov )
btab_split

expect_true( nrow( btab_split$schools ) == 2 )
expect_true( nrow( btab_split$students ) == 3 )
} )



test_that( "single.table works", {

btab = balanceMulti( match.simpleA, single.table = TRUE)
btab

btab2 = balanceMulti( match.simpleA, single.table = TRUE, include.tests = FALSE )
btab2

expect_equal( nrow(btab), nrow(btab2) )
expect_equal( dim(btab2), c( 14, 6 ) )

} )


