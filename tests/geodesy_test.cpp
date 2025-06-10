#include <gtest/gtest.h>
#include "geodesy.h"

TEST(GeodesyTest, RoundTrip) {
    double lat = 10.0, lon = -20.0, h = 1000.0;
    double x, y, z;
    geodeticToGeocentric(lat, lon, h, x, y, z);
    double lat2, lon2, h2;
    geocentricToGeodetic(x, y, z, lat2, lon2, h2);
    EXPECT_NEAR(lat, lat2, 1e-6);
    EXPECT_NEAR(lon, lon2, 1e-6);
    EXPECT_NEAR(h, h2, 1e-3);
}