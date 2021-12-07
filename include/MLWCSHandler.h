//
//  MLWCSHandler.h
//  VAMaximumLikelihood
//
//  Created by Josh Cardenzana on 12/4/14.
//  Copyright (c) 2014 Josh Cardenzana. All rights reserved.
//

#ifndef MLWCSHANDLER_H
#define MLWCSHANDLER_H

// C++ HEADERS
#include <cstring>
#include <stdio.h>
#include <string>
#include "wcs.h"

    
    // Method for converting from gnomonic projected coordinates to world coordinates
    void Projected2Coordinate(double xproj, double yproj,
                            double xTangentPoint, double yTangentPoint,
                            double *xcoord, double *ycoord,
                            const std::string& coordSys="J2000") ;
    // Slightly faster way when needing to call this method MANY times
    void Projected2Coordinate_Quick(double xproj, double yproj,
                            WorldCoor* wcs,
                            double *xcoord, double *ycoord) ;

    // Method for converting from world coordinates to gnomonic projected coordinates
    void Coordinate2Projected(double xcoord, double ycoord,
                              double xTangentPoint, double yTangentPoint,
                              double *xproj, double *yproj,
                              const std::string& coordSys="J2000") ;
    // Slightly faster way when needing to call this method MANY times
    void Coordinate2Projected_Quick(double xcoord, double ycoord,
                              WorldCoor* wcs,
                              double *xproj, double *yproj) ;

    // Method for converting from celestial (J2000) to galactic coordinates
    void Celestial2Galactic(double Ra, double Dec,
                            double *GLon, double *GLat) ;
    
    // Method for converting from galactic to celestial (J2000) coordinates
    void Galactic2Celestial(double GLon, double GLat,
                            double *Ra, double *Dec) ;

    // Method for converting from one system to another
    void Coordinate2Coordinate(double xcoord1, double ycoord1,
                               double *xcoord2, double *ycoord2,
                               const std::string& input_system,
                               const std::string& output_system) ;

    // Method for calculating the angular separation between two points on the celestial sphere
    // Specifically assumes coordinates passed and angular distance returned are all in degrees.
    double AngularSeparation_Deg(double x1coord, double y1coord,
                               double x2coord, double y2coord) ;
    
    // Method for calculating the angular separation between two points on the celestial sphere
    // Specifically assumes coordinates passed and angular distance returned are all in radians.
    double AngularSeparation_Rad(double x1coord, double y1coord,
                               double x2coord, double y2coord) ;
    
    // Method for initializing the wcs object which is needed by some of the above coordinate
    // conversion methods
    WorldCoor* initWCS(double xTangentPoint, double yTangentPoint, const std::string& coordSys) ;
    
    // Methods adapted from ROOT's TMath class
    inline double Pi()       { return 3.14159265358979323846; }
    inline double RadToDeg() { return 180.0 / Pi(); }
    inline double DegToRad() { return Pi() / 180.0; }

#endif /* defined(MLWCSHANDLER_H) */
