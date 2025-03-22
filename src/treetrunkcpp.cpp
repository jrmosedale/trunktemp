#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;
// ** Calculates Astronomical Julian day ** //
int juldayCpp(int year, int month, int day)
{
    double dd = day + 0.5;
    int madj = month + (month < 3) * 12;
    int yadj = year + (month < 3) * -1;
    double j = std::trunc(365.25 * (yadj + 4716)) + std::trunc(30.6001 * (madj + 1)) + dd - 1524.5;
    int b = 2 - std::trunc(yadj / 100) + std::trunc(std::trunc(yadj / 100) / 4);
    int jd = static_cast<int>(j + (j > 2299160) * b);
    return jd;
}
// ** Calculates solar time ** //
double soltimeCpp(int jd, double lt, double lon)
{

    double m = 6.24004077 + 0.01720197 * (jd - 2451545);
    double eot = -7.659 * sin(m) + 9.863 * sin(2 * m + 3.5932);
    double st = lt + (4 * lon + eot) / 60;
    return st;
}
// ** Calculates solar position ** //
std::vector<double> solpositionCpp(double lat, double lon, int year, int month, int day, double lt)
{
    int jd = juldayCpp(year, month, day);
    double st = soltimeCpp(jd, lt, lon);
    // Calculate solar zenith (degrees)
    double latr = lat * M_PI / 180;
    double tt = 0.261799 * (st - 12);
    double dec = (M_PI * 23.5 / 180) * cos(2 * M_PI * ((jd - 159.5) / 365.25));
    double coh = sin(dec) * sin(latr) + cos(dec) * cos(latr) * cos(tt);
    double z = acos(coh) * (180 / M_PI);
    // Calculate solar azimuth (degrees)
    double sh = sin(dec) * sin(latr) + cos(dec) * cos(latr) * cos(tt);
    double hh = atan(sh / sqrt(1 - sh * sh));
    double sazi = cos(dec) * sin(tt) / cos(hh);
    double cazi = (sin(latr) * cos(dec) * cos(tt) - cos(latr) * sin(dec)) /
        sqrt(pow(cos(dec) * sin(tt), 2) + pow(sin(latr) * cos(dec) * cos(tt) - cos(latr) * sin(dec), 2));
    double sqt = 1 - sazi * sazi;
    if (sqt < 0) sqt = 0;
    double azi = 180 + (180 * atan(sazi / sqrt(sqt))) / M_PI;
    if (cazi < 0) {
        if (sazi < 0) {
            azi = 180 - azi;
        }
        else {
            azi = 540 - azi;
        }
    }
    // Define and return output variable
    std::vector<double> solpos(2, 0.0);
    solpos[0] = z;
    solpos[1] = azi;
    return solpos;
}
// ** Calculates solar position vector ** //
// [[Rcpp::export]]
List solpositionCppv(double lat, double lon,
    IntegerVector year, IntegerVector month, IntegerVector day,
    NumericVector hour)
{
    int n = year.size();
    NumericVector zen(n);
    NumericVector azi(n);
    for (int i = 0; i < n; ++i) {
        std::vector<double> sp = solpositionCpp(lat, lon, year[i],
            month[i], day[i], hour[i]);
        zen[i] = sp[0];
        azi[i] = sp[1];  
    }
    Rcpp::List out;
    out["zen"] = Rcpp::wrap(zen);
    out["azi"] = Rcpp::wrap(azi);
    return out;
}
// **  Saturated vapour pressure ** //
double satvapCpp(double tc)
{
    double es;
    if (tc > 0) {
        es = 0.61078 * exp(17.27 * tc / (tc + 237.3));
    }
    else {
        es = 0.61078 * exp(21.875 * tc / (tc + 265.5));
    }
    return es;
}
// **  Dewpoint temperature ** // 
double dewpointCpp(double tc, double ea)
{
    double e0;
    double L;
    double it;
    // Dew point
    if (tc >= 0) {
        e0 = 611.2 / 1000;
        L = (2.501 * pow(10, 6)) - (2340 * tc);
        it = 1 / 273.15 - (461.5 / L) * log(ea / e0);
    }
    // Frost point
    else {
        e0 = 610.78 / 1000;
        L = 2.834 * pow(10, 6);
        it = 1 / 273.16 - (461.5 / L) * log(ea / e0);
    }
    double Tdew = 1 / it - 273.15;
    return Tdew;
}
// ** Calculate absorbed radiation
double radabsoptioncpp(double Rdirdown, double Rdifdown, double Rswup,
    double Rlwup, double Rlwdown, double zen, double azi, double aspect,
    double refl, double em)
{
    // Calculate fraction of radiation absorbed by tree trunk surface
    double index = sin(zen * M_PI / 180) * cos((azi - aspect) * M_PI / 180);
    if (index < 0.0) index = 0.0;
    // Calculate radiation absorbed by tree trunk surface
    double Rabs_sw = (1 - refl) * (index * Rdirdown + 0.5 * (Rdifdown + Rswup));
    double Rabs_lw = em * 0.5 * (Rlwup + Rlwdown);
    double Rabs = Rabs_sw + Rabs_lw;
    return Rabs;
}
// ** Calculates steady state temperature of tree trunk surface ** //
double PenmonMonteithcpp(double Rabs, double tair, double windspeed, 
    double relhum, double pk, double em, double treeradius, double dT, double surfwet, bool cap)
{
    // Calculate Pseudo - emmited radiation 
    double sb = 5.67 * pow(10.0, -8.0);
    double Rema = em * sb * pow(tair + 273.15, 4.0);
    // Calculate convective conductance
    double d = treeradius * 2.0; // characteristic dimension
    double gHa = 0.135 * sqrt(windspeed / d); // forced convection
    double gHfr = 0.033 * pow(abs(dT) / d, 0.25); // free convection
    if (gHfr > gHa) gHa = gHfr; // set forced to free if greater
    // Calculate radiative conductance
    double cp = pow(2.0, -5.0) * tair * tair + 0.002 * tair + 29.119;
    double gr = (4.0 / cp) * sb * pow(tair + 273.15, 3.0);
    double gHr = gr + gHa;
    // Calculate vapour pressure deficit
    double la = 45068.7 - 42.8428 * tair;
    if (tair < 0) la = 51078.69 - 4.338 * tair - 0.06367 * tair * tair;
    double es = satvapCpp(tair);
    double Da = es * (1 - 0.01 * relhum);
    // Calculate slope of saturated vapour pressure curve
    double De = satvapCpp(tair + 0.5) - satvapCpp(tair - 0.5);
    // Calculate steady-state temperature
    double Ts = tair + ((Rabs - Rema - la * (gHa / pk) * Da * surfwet)
        / (cp * gHr + la * (gHa / pk) * De * surfwet));
    if (cap) {
        double Tdew = dewpointCpp(tair, es * relhum / 100.0);
        if (Ts < Tdew) Ts = Tdew;
    }
    // return temperature difference
    double dTS = Ts - tair;
    return dTS;
}
// Run steady-state model through time
// [[Rcpp::export]]   
// ** Calculates steady state temperature of tree trunk surface ** //
NumericVector SteadyState(DataFrame microclim, double refl,
    double em, double treeradius, double surfwet, double aspect, 
    bool cap)
{
    // Extract data from data.frame
    IntegerVector year = microclim["year"];
    IntegerVector month = microclim["month"];
    IntegerVector day = microclim["day"];
    NumericVector hour = microclim["hour"];
    NumericVector tair = microclim["tair"];
    NumericVector relhum = microclim["relhum"];
    NumericVector windspeed = microclim["windspeed"];
    NumericVector Rdirdown = microclim["Rdirdown"];
    NumericVector Rdifdown = microclim["Rdifdown"];
    NumericVector Rlwdown = microclim["Rlwdown"];
    NumericVector Rswup = microclim["Rswup"];
    NumericVector Rlwup = microclim["Rlwup"];
    NumericVector pk = microclim["pk"];
    NumericVector zen = microclim["zen"];
    NumericVector azi = microclim["azi"];
    int n = year.size();
    NumericVector ts(n);
    double dT = 0.0;
    for (int i = 0; i < n; ++i) {
        // Calculate absorbed radiation
        double Rabs = radabsoptioncpp(Rdirdown[i], Rdifdown[i], Rswup[i],
            Rlwup[i], Rlwdown[i], zen[i], azi[i], aspect, refl, em);
        dT = PenmonMonteithcpp(Rabs, tair[i], windspeed[i],
            relhum[i], pk[i], em, treeradius, dT, surfwet, cap);
        ts[i] = tair[i] + dT;
    }
    return(ts);
}
// ** Calculates energy balance of tree trunk surface ** //
double energybalancecpp(double Rabs, double tair, double tsurf, double windspeed,
    double relhum, double pk, double em, double treeradius, double dT,
    double surfwet)
{
    // Calculate emitted radiation
    double sb = 5.67 * pow(10.0, -8.0);
    double Rem = em * sb * pow(tsurf + 273.15, 4.0);
    // Calculate sensible heat flux
    double cp = pow(2.0, -5.0) * tair * tair + 0.002 * tair + 29.119;
    double d = treeradius * 2.0; // characteristic dimension
    double gHa = 0.135 * sqrt(windspeed / d); // forced convection
    double gHfr = 0.033 * pow(abs(dT) / d, 0.25); // free convection
    if (gHfr > gHa) gHa = gHfr; // set forced to free if greater
    double H = cp * gHa * (tsurf - tair);
    // Calculate latent heat flux
    double la = 45068.7 - 42.8428 * tair;
    if (tair < 0) la = 51078.69 - 4.338 * tair - 0.06367 * tair * tair;
    double es = satvapCpp(tsurf);
    double ea = satvapCpp(tair) * (relhum / 100.0);
    double L = (la * gHa / pk) * (es - ea) * surfwet;
    // Calculate energy balance
    double Ba = Rabs - Rem - H - L;
    return Ba; // W/m^2
}
// Calculate tree trunk temperatures in one step
// [[Rcpp::export]]   
NumericMatrix onestepcpp(double Rdirdown, double Rdifdown, double Rswup,
    double Rlwup, double Rlwdown, double zen, double azi,
    double tair, double windspeed, double relhum, double pk,
    double refl, double em, double treeradius, double surfwet,
    NumericMatrix ptemps, NumericVector kwood, NumericVector ldist,
    NumericVector sdist, NumericVector heatcap, int timestep)
{
    int nsegs = ptemps.ncol();
    int nlyrs = ptemps.nrow();
    NumericMatrix temps = ptemps;
    // Calculate temperature out outer layer
    for (int i = 0; i < nsegs; ++i) {
        // ** Calculate radiation absorbed by outer layer
        double aspect = (static_cast<double>(i) / static_cast<double>(nsegs)) * 360.0;
        double Rabs = radabsoptioncpp(Rdirdown, Rdifdown, Rswup,
            Rlwup, Rlwdown, zen, azi, aspect, refl, em);
        // ** Energy balance of outer layer in W/m^2
        double dTp = ptemps(i, 0) - tair;
        double eb = energybalancecpp(Rabs, tair, ptemps(i, 0), windspeed,
            relhum, pk, em, treeradius, dTp, surfwet);
        // ** Convert energy balance to Joules
        eb = static_cast<double>(timestep) * eb * 2 *
            M_PI * treeradius / static_cast<double>(nsegs);
        // ** Calculate transient temperature change
        double dTT = eb / heatcap[0];
        temps(0, i) = ptemps(0, i) + dTT;
        // ** Calculate steady-state temperature change
        double dTS = PenmonMonteithcpp(Rabs, tair, windspeed, relhum, pk,
            em, treeradius, dTp, surfwet, true);
        // ** Cap at steady state if necessary
        double tdif = temps(0, i) - tair;
        if (abs(tdif) > abs(dTS)) temps(0, i) = tair + dTS;
    }
    // exchange heat between layers
    for (int i = 0; i < (nlyrs - 1); ++i) {
        for (int j = 0; j < nsegs; ++j) {
            // initial temperature difference between outer layer and layer below
            double dTl = temps(i, j) - temps(i + 1, j);
            // energy passed between outer layer and layer below
            double ep = (kwood[i] / ldist[i]) * timestep * dTl;
            // transient temperature change
            double dT1 = ep / heatcap[i];
            double dT2 = ep / heatcap[i + 1];
            // steady state temperature change
            double dS1 = (dTl * heatcap[i + 1]) / (heatcap[i] + heatcap[i + 1]);
            double dS2 = (dTl * heatcap[i]) / (heatcap[i] + heatcap[i + 1]);
            // set transient to steady state if greater
            if (abs(dT1) > abs(dS1)) dT1 = dS1;
            if (abs(dT2) > abs(dS2)) dT2 = dS2;
            // swap heat
            temps(i, j) = temps(i, j) - dT1;
            temps(i + 1, j) = temps(i + 1, j) + dT2;
        }
    }
    // exchange heat between segments
    for (int i = 0; i < nlyrs; ++i) {
       for (int j = 0; j < nsegs; ++j) {
            double dTs = 0.0;
            if (j < (nsegs - 1)) {
                // initial temperature difference between segment and segment to right
                dTs = temps(i, j) - temps(i, j + 1);
            }
            else {
                // initial temperature difference between segment and segment to right
                dTs = temps(i, j) - temps(i, 0);
            }
            // energy passed between segment and segment to right
            double ep = (kwood[i] / sdist[i]) * timestep * dTs;
            // transient temperature change
            double dTr = ep / heatcap[i];
            // steady-state temperature change
            double dTS = 0.5 * dTs;
            // set transient to steady state if greater
            if (abs(dTr) > abs(dTS)) dTr = dTS;
            // swap heat
            temps(i, j) = temps(i, j) - dTr;
            if (j < (nsegs - 1)) {
                temps(i, j + 1) = temps(i, j + 1) + dTr;
            }
            else {
                temps(i, 0) = temps(i, 0) + dTr;
            }
       }
    }
    return temps;
}
// Run first day repeatedly X times
// [[Rcpp::export]]   
NumericMatrix burnincpp(DataFrame microclim, int reps,
    double refl, double em, double treeradius, double surfwet,
    NumericMatrix ptemps, NumericVector kwood, NumericVector ldist,
    NumericVector sdist, NumericVector heatcap, int timestep)
{
    // Extract Data.Frame variables
    IntegerVector year = microclim["year"];
    IntegerVector month = microclim["month"];
    IntegerVector day = microclim["day"];
    NumericVector hour = microclim["hour"];
    NumericVector tair = microclim["tair"];
    NumericVector relhum = microclim["relhum"];
    NumericVector windspeed = microclim["windspeed"];
    NumericVector Rdirdown = microclim["Rdirdown"];
    NumericVector Rdifdown = microclim["Rdifdown"];
    NumericVector Rlwdown = microclim["Rlwdown"];
    NumericVector Rswup = microclim["Rswup"];
    NumericVector Rlwup = microclim["Rlwup"];
    NumericVector pk = microclim["pk"];
    NumericVector zen = microclim["zen"];
    NumericVector azi = microclim["azi"];
    int tms = 24 * 3600 / timestep;
    NumericMatrix temps = ptemps;
    for (int rep = 0; rep < reps; ++rep) {
        for (int i = 0; i < tms; ++i) {
            temps = onestepcpp(Rdirdown[i], Rdifdown[i], Rswup[i],
                Rlwup[i], Rlwdown[i], zen[i], azi[i],
                tair[i], windspeed[i], relhum[i], pk[i],
                refl, em, treeradius, surfwet, temps, kwood, ldist,
                sdist, heatcap, timestep);

        }
    }
    return (temps);
}
// Run model
// [[Rcpp::export]]   
NumericMatrix runmodelcpp(DataFrame microclim, int n,
    double refl, double em, double treeradius, double surfwet,
    NumericMatrix ptemps, NumericVector kwood, NumericVector ldist,
    NumericVector sdist, NumericVector heatcap, int timestep)
{
    // Extract Data.Frame variables
    IntegerVector year = microclim["year"];
    IntegerVector month = microclim["month"];
    IntegerVector day = microclim["day"];
    NumericVector hour = microclim["hour"];
    NumericVector tair = microclim["tair"];
    NumericVector relhum = microclim["relhum"];
    NumericVector windspeed = microclim["windspeed"];
    NumericVector Rdirdown = microclim["Rdirdown"];
    NumericVector Rdifdown = microclim["Rdifdown"];
    NumericVector Rlwdown = microclim["Rlwdown"];
    NumericVector Rswup = microclim["Rswup"];
    NumericVector Rlwup = microclim["Rlwup"];
    NumericVector pk = microclim["pk"];
    NumericVector zen = microclim["zen"];
    NumericVector azi = microclim["azi"];
    NumericMatrix temps = ptemps;
    for (int i = 0; i < n; ++i) {
        temps = onestepcpp(Rdirdown[i], Rdifdown[i], Rswup[i],
            Rlwup[i], Rlwdown[i], zen[i], azi[i],
            tair[i], windspeed[i], relhum[i], pk[i],
            refl, em, treeradius, surfwet, temps, kwood, ldist,
            sdist, heatcap, timestep);
    }
    return (temps);
}
// Run model and extract data for one segment and layer
// [[Rcpp::export]]   
NumericVector runmodelthroughtime(DataFrame microclim, int seg, int lyr,
    double refl, double em, double treeradius, double surfwet,
    NumericMatrix ptemps, NumericVector kwood, NumericVector ldist,
    NumericVector sdist, NumericVector heatcap, int timestep)
{
    // Extract Data.Frame variables
    IntegerVector year = microclim["year"];
    IntegerVector month = microclim["month"];
    IntegerVector day = microclim["day"];
    NumericVector hour = microclim["hour"];
    NumericVector tair = microclim["tair"];
    NumericVector relhum = microclim["relhum"];
    NumericVector windspeed = microclim["windspeed"];
    NumericVector Rdirdown = microclim["Rdirdown"];
    NumericVector Rdifdown = microclim["Rdifdown"];
    NumericVector Rlwdown = microclim["Rlwdown"];
    NumericVector Rswup = microclim["Rswup"];
    NumericVector Rlwup = microclim["Rlwup"];
    NumericVector pk = microclim["pk"];
    NumericVector zen = microclim["zen"];
    NumericVector azi = microclim["azi"];
    NumericMatrix temps = ptemps;
    // How many time periods
    int n = year.size();
    NumericVector tout(n);
    for (int i = 0; i < n; ++i) {
        temps = onestepcpp(Rdirdown[i], Rdifdown[i], Rswup[i],
            Rlwup[i], Rlwdown[i], zen[i], azi[i],
            tair[i], windspeed[i], relhum[i], pk[i],
            refl, em, treeradius, surfwet, temps, kwood, ldist,
            sdist, heatcap, timestep);
        tout[i] = temps(lyr, seg);
    }
    return tout;
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ******************** Data processing functions ******************* //
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// [[Rcpp::export]]
NumericVector satvapCpp(NumericVector tc) {
    IntegerVector dims = tc.attr("dim");
    // Get dimensions
    int rows = dims[0];
    int cols = dims[1];
    int hrs = dims[2];
    int n = rows * cols * hrs;
    NumericVector es(n);
    for (int i = 0; i < n; ++i) {
        if (!NumericVector::is_na(tc[i])) {
            if (tc[i] > 0) {
                es[i] = 0.61078 * exp(17.27 * tc[i] / (tc[i] + 237.3));
            }
            else {
                es[i] = 0.61078 * exp(21.875 * tc[i] / (tc[i] + 265.5));
            }
        }
        else {
            es[i] = NA_REAL;
        }
    }
    es.attr("dim") = IntegerVector::create(rows, cols, hrs);
    return es;
}
// ** Calculates solar zenith in radians ** //
double solzenCpp(double jd, double st, double lat, double lon)
{
    double latr = lat * M_PI / 180;
    double tt = 0.261799 * (st - 12);
    double dec = (M_PI * 23.5 / 180) * cos(2 * M_PI * ((jd - 159.5) / 365.25));
    double coh = sin(dec) * sin(latr) + cos(dec) * cos(latr) * cos(tt);
    double z = acos(coh);
    return z;
}
// [[Rcpp::export]]
NumericVector solarindexarray(IntegerVector year, IntegerVector month,
    IntegerVector day, NumericVector hour, NumericMatrix lats, NumericMatrix lons)
{
    // Get dims
    int nrow = lats.nrow();
    int ncol = lats.ncol();
    int n = year.size();
    // Calculate julian day vector
    NumericVector jd(n);
    for (int k = 0; k < n; ++k) jd[k] = juldayCpp(year[k], month[k], day[k]);
    // Get solar index
    NumericVector si(nrow * ncol * n);
    for (int i = 0; i < nrow; ++i) {
        for (int j = 0; j < ncol; ++j) {
            for (int k = 0; k < n; ++k) {
                int idx = i + nrow * (j + ncol * k);
                double st = soltimeCpp(jd[k], hour[k], lons(i, j));
                double z = solzenCpp(jd[k], st, lats(i, j), lons(i, j));
                si[idx] = cos(z);
                if (si[idx] < 0.0) si[idx] = 0.0;
            }
        }
    }
    si.attr("dim") = IntegerVector::create(nrow, ncol, n);
    return si;
}
// ** Blend met office and e.g. era5 climate data
// [[Rcpp::export]]
NumericVector blendtempCpp(NumericVector tasmin, NumericVector tasmax,
    NumericVector temp)
{
    // Do dimension and length checks
    IntegerVector dims1 = tasmin.attr("dim");
    IntegerVector dims2 = tasmax.attr("dim");
    IntegerVector dims3 = temp.attr("dim");
    for (int i = 0; i < 3; ++i) {
        if (dims1[i] != dims2[i]) stop("Dimensions of tasmin and tasmax don't match");
    }
    for (int i = 0; i < 2; ++i) {
        if (dims1[i] != dims3[i]) stop("xy dimensions of tasmin and temp don't match");
    }
    if (dims1[2] * 24 != dims3[2]) stop("temp must have 24 times more entries than tasmin");
    // Get dimensions
    int rows = dims1[0];
    int cols = dims1[1];
    int days = dims1[2];
    int hours = days * 24;
    NumericVector tcout(rows * cols * hours);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            if (!NumericVector::is_na(temp[i + rows * (j + cols * 0)])) {
                // Extract vector of hourly temperatures
                NumericVector tempv(hours);
                for (int k = 0; k < hours; ++k) {
                    tempv[k] = temp[i + rows * (j + cols * k)];
                }
                // Extract vector of daily mins and maxes temperatures
                NumericVector tmxd(days);
                NumericVector tmnd(days);
                for (int k = 0; k < days; ++k) {
                    tmnd[k] = tasmin[i + rows * (j + cols * k)];
                    tmxd[k] = tasmax[i + rows * (j + cols * k)];
                }
                int index1 = 0;
                int index2 = 0;
                for (int day = 0; day < days; ++day) {
                    // For each day, calculate fraction of dtr in temp
                    double tmx = -273.15;
                    double tmn = 273.15;
                    for (int hr = 0; hr < 24; ++hr) {
                        if (tempv[index1] < tmn) tmn = tempv[index1];
                        if (tempv[index1] > tmx) tmx = tempv[index1];
                        ++index1;
                    }
                    double dtr1 = tmxd[day] - tmnd[day]; // dtr of met office
                    double dtr2 = tmx - tmn; // dtr of era5
                    for (int hr = 0; hr < 24; ++hr) {
                        double tfrac = (tempv[index2] - tmn) / dtr2;
                        double tc = tfrac * dtr1 + tmnd[day];
                        tcout[i + rows * (j + cols * index2)] = tc;
                        ++index2;
                    } // end hours
                } // end days
            } // NA check
            else {
                for (int k = 0; k < hours; ++k) {
                    tcout[i + rows * (j + cols * k)] = NA_REAL;
                } // end hours
            } // end NA
        } // end col
    } // end row
    tcout.attr("dim") = IntegerVector::create(rows, cols, hours);
    return tcout;
}