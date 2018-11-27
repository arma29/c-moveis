library(geosphere)

locations <- data.frame(
    latitude = c(
        -8.068361111, -8.075916667, -8.076361111, -8.0593334,
        -8.066, -8.06458333333333
    ),
    longitude = c(
        -34.892722222, -34.894611111, -34.908, -34.9073435,
        -34.8894444444444, -34.8945833333333

    ),
    distance = c(124.4, 188.6, 709, 2350, 806.8, 706.2)
)

# Use average as the starting point
fit <- nls(
    distance ~ distm(
        data.frame(longitude, latitude),
        c(fitLongitude, fitLatitude)
    ),
    data = locations,
    start = list(
        fitLongitude=mean(locations$longitude),
        fitLatitude=mean(locations$latitude)
    ),
    control = list(maxiter = 1000, tol = 1e-02)
)

# Result
latitude <- summary(fit)$coefficients[2]
longitude <- summary(fit)$coefficients[1]

paste(latitude, longitude, sep=",")

#
