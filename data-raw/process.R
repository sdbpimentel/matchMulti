
# Make small Catholic dataset for faster illustrations

# Trim data to speed up example
data( "catholic_schools" )
minischool = dplyr::filter( catholic_schools, 
                                  female_mean > 0.45,
                                  female_mean < 0.60,
                                  acad > quantile( acad, 0.25 ),
                                  acad < quantile( acad, 0.75 ) )



usethis::use_data(minischool, overwrite = TRUE)


if ( FALSE ) {
  sinew::makeOxygen(minischool, add_fields = "source")
}