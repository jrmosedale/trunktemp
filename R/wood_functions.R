
#' Calculate trunk parameters from species table for use with trunk emulator model
#'
#' @param spp_params = dataframe of parameters by species
#' @param sp = species
#' @param tradius = trunk radius
#' @param outerlayers =  vector of the width of outerlayers
#'
#' @return
#' @export
#'
#' @examples
trunk_parameters<-function(spp_params,sp,tradius, outerlayers=c(0.015,0.02,0.03,0.035)){
  # Calculate thermal properties from moisture content etc
  sw_mc<-spp_params['sw_mc',sp]
  hw_mc<-spp_params['hw_mc',sp]
  bk_mc<-spp_params['bk_mc',sp]
  fsp<-spp_params['fsp',sp]
  bk_fsp<-spp_params['bk_fsp',sp]

  # Calculate thermal conductivity
  swK<-calc_thcond(Gb=spp_params['Gb',sp],x=sw_mc)
  hwK<-calc_thcond(Gb=spp_params['Gb',sp],x=hw_mc)
  bkK<-calc_thcond(Gb=spp_params['bk_Gb',sp],x=bk_mc)

  # Calculate spec heat capacity (sapwood outer and inner other layers) kJ/Kg/K
  swcph<-calc_sph(sw_mc,cp0=NA,T=283,fsat=fsp)
  hwcph<-calc_sph(hw_mc,cp0=NA,T=283,fsat=fsp)
  bkcph<-calc_sph(bk_mc,cp0=NA,T=283,fsat=bk_fsp)

  # Calculate wood density (green)
  swrho<-calc_density(Gb=spp_params['Gb',sp],x=sw_mc)
  hwrho<-calc_density(Gb=spp_params['Gb',sp],x=hw_mc)
  bkrho<-calc_density(Gb=spp_params['bk_Gb',sp],x=bk_mc)

  # Calculate layer parameters etc
  nlyr<-8
  nsegs<-16

  # Width of layers - probes shallow = 2 (2.5 cm) deep = layer 3 (5cm)
  # Fixed outer layers then calculate others
  #outerlayers<-c(0.015,0.02,0.03,0.035)
  #outerlayers<-c(0.005,0.025,0.03,0.035)
  thick<-sqrt(c(5:8))
  thick <- (thick/ sum(thick)) * (tradius-sum(outerlayers))
  layer_widths<-c(outerlayers,thick) # sum to tradius
  if( abs(sum(layer_widths)-tradius) > 1e-5) stop("Layer widths do NOT sum to tree radius!!!")

  if(sp=='nsp'){
    cs <- c(bkcph,swcph,swcph,rep(hwcph,(nlyr-3)))
    rho <- c(bkrho,swrho,swrho,rep(hwrho,(nlyr-3)))
    k <- c(bkK,swK,swK,rep(hwK,(nlyr-3)))
  }
  if(sp %in% c('oak','sycamore')){
    cs <- c(bkcph,swcph,rep(hwcph,(nlyr-2)))
    rho <- c(bkrho,swrho,rep(hwrho,(nlyr-2)))
    k <- c(bkK,swK,rep(hwK,(nlyr-2)))
  }

  treeparams <- list(
    tradius = tradius,
    layer_widths=layer_widths,
    cs = cs,
    rho = rho,
    refl = spp_params[["bk_swref",sp]],
    em = spp_params[["bk_tem",sp]],
    k = k,
    surfwet = 1 )

  return(treeparams)
}


#' Calculate thermal conductivity
#'
#' @param Gb = wood basic specific gravity (oven dry mass, green volume in g/cm3?)
#' @param x = moisture content of wood in %
#' @param MCfs = moisture content at fibre saturation as %
#'
#' @return
#' @export
#'
#' @examples
#' calc_thcond(Gb=0.6,x=25); calc_thcond(Gb=0.37,x=25)
calc_thcond<-function(Gb,x=25,MCfs=30){
  if(x>25) warning("Calculation of thermal conductivity unreliable at moisture contents > 25%")
  A<-0.01864
  B<-0.1941
  C<-0.004064
  # Convert Spec Gravity to moisture content x
  # Assumes total vol shrinkage S0 estimated from Gb
  Gx<-Gb/(1-0.265*Gb*(1-x/MCfs)) #Eq 4.13
  # Calc conductivity from Spec gravity at MCx
  k<-Gx*(B+C*x)+A
  return(k)
}


#' Calculate Specific heat capacity of wood at different moisture contents
#' The heat capacity of wood depends on the temperature and moisture content
#'  of the wood but is practically independent of density or species. (Glass & Zelinka)
#' @param x = Moisture content of wood as %
#' @param cp0 = spec heat cap of dry wood
#' @param T = temperature
#' @param fsat = moisture content at fibre saturation as %
#'
#' @return
#' @export
#'
#' @examples
calc_sph<-function(x,cp0=NA,T=283,fsat=30){
  if (is.na(cp0)) cp0<-0.1031+0.003867*T
  # Calculate correction factor Ac for <fibre sat moisture content
  cpw<-4.18 # Water Sp H Cap
  b1<- -0.06191
  b2<-0.000236
  b3<- -0.000133
  # Moisture content < fibre sat
  if(x<=fsat){
    Ac<-x*(b1+b2*T+b3*x) # Correction factor for <fibre sat moisture content
    cpx<-(cp0+cpw*x/100) / (1+x/100) + Ac
  }
  if(x>fsat){
    Ac<-fsat*(b1+b2*T+b3*fsat)
    cpfs<-(cp0+cpw*fsat/100) / (1+fsat/100) + Ac
    cpx<- cpfs*(fsat/x) + (cpw*(1-fsat/x))
  }
  return(cpx)
}

#' Calculate wood density
#'
#' @param Gb = wood basic specific gravity (oven dry mass, green volume in g/cm3?)
#' @param x = = Moisture content of wood as %
#' @param MCfs = moisture content at fibre sat (%)
#' @param rhow in cm/g3 =1
#'
#' @returns density of wood at given moisture content in kg/m3 units
#' @export
calc_density<-function(Gb,x,MCfs=30, rhow=1){
  Gx<-Gb/(1-0.265*Gb*(1-x/MCfs))
  rho<-rhow*Gx*(1+x/100)
  return(rho*1000)
}

