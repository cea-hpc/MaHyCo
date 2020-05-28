
KOKKOS_INLINE_FUNCTION
double fluxLimiter(int projectionLimiterId, double r) 
{
  if (projectionLimiterId == options->minmod) 
    {
      return MathFunctions::max(0.0, MathFunctions::min(1.0, r));
    }
  else if (projectionLimiterId == options->superBee) 
    {
      //if (r !=0.) std::cout << " r et Limiter  " << r << "  "<< MathFunctions::max(0.0, MathFunctions::max(MathFunctions::min(2.0 * r, 1.0), MathFunctions::min(r, 2.0)))  << std::endl;
      return MathFunctions::max(0.0, MathFunctions::max(MathFunctions::min(2.0 * r, 1.0), MathFunctions::min(r, 2.0)));
    }
  else if (projectionLimiterId == options->vanLeer) 
    {
      if (r <= 0.0) 
	return 0.0;
      else 
	return 2.0 * r / (1.0 + r);
    }
  else return 0.0; // ordre 1
}
KOKKOS_INLINE_FUNCTION
double fluxLimiterPP(int projectionLimiterId, double gradplus, double gradmoins, double y0, double yplus, double ymoins, double h0, double hplus , double hmoins) 
{
  double grady,gradM, gradMplus,gradMmoins;
  // limitation rupture de pente (formule 16 si on utilise pas le plateau pente) 
  if (gradplus * gradmoins < 0.0) return 0.;
	  
  if (projectionLimiterId == options->minmodG) // formule 9c
    {		  
      if ((yplus-ymoins) > 0.) grady = MathFunctions::min(fabs(gradplus), fabs(gradmoins));
      else grady = - MathFunctions::min(fabs(gradplus), fabs(gradmoins));
    }
  else if (projectionLimiterId == options->superBeeG) // formule 9g
    {
      if ((yplus-ymoins) > 0.) grady = MathFunctions::max(fabs(gradplus), fabs(gradmoins));
      else grady = - MathFunctions::max(fabs(gradplus), fabs(gradmoins));
    }
  else if (projectionLimiterId == options->vanLeerG) // formule 9e
    {
      double lambdaplus = (h0/2+hplus) / (h0+hplus+hmoins);
      double lambdamoins = (h0/2+hmoins) / (h0+hplus+hmoins);
      if ((lambdaplus * gradplus + lambdamoins * gradmoins) != 0.) 
	{
	  grady = gradplus * gradmoins / (lambdaplus * gradplus + lambdamoins * gradmoins);
	}		  
      else grady = 0.;
    }
  else if (projectionLimiterId == options->arithmeticG)
    {
      double lambdaplus = (h0/2+hplus) / (h0+hplus+hmoins);
      double lambdamoins = (h0/2+hmoins) / (h0+hplus+hmoins);
      grady = lambdamoins * gradplus + lambdaplus * gradmoins;
    }
	  
  // limitation simple-pente (formule 10)
  gradMplus = gradplus * (h0 + hplus)/h0;
  gradMmoins = gradmoins * (h0 + hmoins)/h0;
  gradM = MathFunctions::min(fabs(gradMplus), fabs(gradMmoins));
  if ((yplus-ymoins) > 0.) grady = MathFunctions::min(fabs(gradM), fabs(grady));
  else grady = - MathFunctions::min(fabs(gradM), fabs(grady));
	  
  return grady;
}
KOKKOS_INLINE_FUNCTION
double computeY0(int projectionLimiterId, double y0, double yplus, double ymoins, double h0, double hplus, double hmoins, int type) 
{
  // retourne {{y0plus, y0moins}}
  double y0plus=0., y0moins=0.;
  if (projectionLimiterId == options->minmodG || projectionLimiterId == options->minmod) // minmod
    {
      y0plus = yplus;
      y0moins = ymoins;  
    }
  else if (projectionLimiterId == options->superBeeG || projectionLimiterId == options->superBee) // superbee
    {
      y0plus =  ((h0+hmoins)*yplus+h0*ymoins)/(2*h0+hmoins);
      y0moins = ((h0+hplus)*ymoins+h0*yplus)/(2*h0+hplus);
    }
  else if (projectionLimiterId == options->vanLeerG || projectionLimiterId == options->vanLeer) // vanleer
    {
      double a = MathFunctions::min(yplus, ymoins);
      double b = MathFunctions::max(yplus, ymoins);
      double xplus = (h0*h0+3*h0*hmoins+2*hmoins*hmoins)*yplus;
      double xmoins = (h0*h0+3*h0*hplus+2*hplus*hplus)*ymoins;
      xplus += (h0*h0-h0*hplus-2*hplus*hplus+2*h0*hmoins)*ymoins;
      xmoins += (h0*h0-h0*hmoins-2*hmoins*hmoins+2*h0*hplus)*yplus;
      xplus /= (2*h0*h0+5*h0*hmoins+2*hmoins*hmoins-h0*hplus-2*hplus*hplus);
      xmoins /= (2*h0*h0+5*h0*hplus+2*hplus*hplus-h0*hmoins-2*hmoins*hmoins);

      y0plus = MathFunctions::min(MathFunctions::max(xplus,a), b);
      y0moins = MathFunctions::min(MathFunctions::max(xmoins,a), b);		  		  
    }
  else if (projectionLimiterId == options->arithmeticG)
    {
      y0plus =  ((h0+hmoins+hplus)*yplus+h0*ymoins)/(2*h0+hmoins+hplus);
      y0moins = ((h0+hmoins+hplus)*ymoins+h0*yplus)/(2*h0+hmoins+hplus);
    }
  else if (projectionLimiterId == 3000) 
    {
      y0plus = yplus;
      y0moins = ymoins;  
    }
  if (type == 0) return y0plus;
  else if (type == 1) return y0moins;
  else return 0.0; // lancer forcement avec type 0 ou 1 mais warning compile
}
KOKKOS_INLINE_FUNCTION
double computexgxd(double y0, double yplus, double ymoins, double h0, double y0plus, double y0moins, int type) 
{
  // retourne {{xg, xd}}
  double xd = 0., xg=0.;
  double xplus = 1.;
  if (abs(y0plus - yplus) > options->threshold) xplus = (y0-yplus)/(y0plus - yplus) -1./2.;
  double xmoins = 1.;
  if (abs(y0moins - ymoins) > options->threshold) xmoins = (y0-ymoins)/(y0moins - ymoins) -1./2.;
  xd = +h0 * MathFunctions::min(MathFunctions::max(xplus,-1./2.), 1./2.);
  xg = -h0 * MathFunctions::min(MathFunctions::max(xmoins,-1./2.), 1./2.);
  if (type == 0)  return xg;
  else if (type == 1)  return xd;
  else return 0.0; // lancer forcement avec type 0 ou 1 mais warning compile
}
// KOKKOS_INLINE_FUNCTION
// double computeygyd(double y0, double yplus, double ymoins, double h0, double grady, int type) 
// {
//   // retourne {{yg, yd}}
//   double yd, yg;
//   double xtd = y0+h0/2*grady;
//   double xtg = y0-h0/2*grady;
//   double a = MathFunctions::min(yplus, ymoins);
//   double b = MathFunctions::max(yplus, ymoins);
//   yd = MathFunctions::min(MathFunctions::max(xtd,a),b);
//   yg = MathFunctions::min(MathFunctions::max(xtg,a),b);
//   if (type == 0)  return yg;
//   else if (type == 1)  return yd;
//   else return 0.0; // lancer forcement avec type 0 ou 1 mais warning compile
// }
KOKKOS_INLINE_FUNCTION
double computeygyd(double y0, double yplus, double ymoins, double h0, double y0plus, double y0moins, double grady, int type) 
{
  // retourne {{yg, yd}}
  double yd, yg;
  double xtd = y0+h0/2*grady;
  double xtg = y0-h0/2*grady;
  double ad = MathFunctions::min(yplus, 2.*y0moins-ymoins);
  double bd = MathFunctions::max(yplus, 2.*y0moins-ymoins); 
  double ag = MathFunctions::min(ymoins, 2.*y0plus-yplus);
  double bg = MathFunctions::max(ymoins, 2.*y0plus-yplus);
  yd = MathFunctions::min(MathFunctions::max(xtd,ad),bd);
  yg = MathFunctions::min(MathFunctions::max(xtg,ag),bg);
  if (type == 0)  return yg;
  else if (type == 1)  return yd;
  else return 0.0; // lancer forcement avec type 0 ou 1 mais warning compile
}  
KOKKOS_INLINE_FUNCTION
double INTY(double X, double x0, double y0, double x1, double y1) 
{
  double flux=0.;
  double Xbar = MathFunctions::min(MathFunctions::max(x0, X),x1);
  //std::cout << " Xbar  " << Xbar << std::endl; 
  if (abs(x1-x0) > 1.e-14) flux = (y0 + 0.5* ((Xbar - x0)/(x1-x0)) * (y1-y0) ) * (Xbar -x0);
  return flux;
}
KOKKOS_INLINE_FUNCTION
double INT2Y(double X, double x0, double y0, double x1, double y1) 
{
  double flux=0.;
  //std::cout << " x0 " << x0 << std::endl;
  //std::cout << " x1 " << x1 << std::endl;
  if (abs(x1-x0) > 1.e-14) {
    double eta = MathFunctions::min(MathFunctions::max(0., (X-x0)/(x1-x0)),1.);
    // std::cout << " eta " << eta << std::endl;
    flux = (y0 + 0.5 * eta * (y1-y0) ) * (x1 -x0) * eta;
  } 
  return flux;
}
template<size_t d>
KOKKOS_INLINE_FUNCTION
RealArray1D<d> computeAndLimitGradPhi(int projectionLimiterId, RealArray1D<d> gradphiplus, RealArray1D<d> gradphimoins,
				      RealArray1D<d> phi, RealArray1D<d> phiplus, RealArray1D<d> phimoins, double h0, double hplus, double hmoins) 
{
  RealArray1D<d> res;
  if (projectionLimiterId < options->minmodG) 
    {
      // std::cout << " Passage gradient limite Classique " << std::endl;
      for (size_t i=0; i<d; i++) {
	res[i] = (fluxLimiter(projectionLimiterId, divideNoExcept(gradphiplus[i], gradphimoins[i])) * gradphimoins[i] + fluxLimiter(projectionLimiterId, divideNoExcept(gradphimoins[i], gradphiplus[i])) * gradphiplus[i]) / 2.0;		    
      }
      return res;
    }
  else
    {
		  
      // std::cout << " Passage gradient limite Genéralisé " << std::endl;
      for (size_t i=0; i<d; i++) {
	res[i] = fluxLimiterPP(projectionLimiterId, gradphiplus[i], gradphimoins[i], phi[i], phiplus[i], phimoins[i], h0, hplus, hmoins);
      }
      return res;
    }
}
template<size_t d>
KOKKOS_INLINE_FUNCTION
RealArray1D<d> computeIntersectionPP(RealArray1D<d> gradphi,
				     RealArray1D<d> phi, RealArray1D<d> phiplus, RealArray1D<d> phimoins, double h0, double hplus, double hmoins,
				     double face_normal_velocity, double deltat_n, int type, int cell, double flux_threhold) 
{
  RealArray1D<d> Flux = options->Uzero;
  double y0plus, y0moins, xd, xg, yd, yg;
  double flux1, flux2, flux3, flux1m, flux2m, flux3m;	  
  double partie_positive_v = 0.5*(face_normal_velocity + abs(face_normal_velocity))*deltat_n;
  if (partie_positive_v == 0.) return Flux;
  int cas_PP =0;
  for (size_t i=0; i<nbmatmax; i++) {
    // calcul des seuils y0plus, y0moins pour cCells
    y0plus  = computeY0(options->projectionLimiterId, phi[i], phiplus[i], phimoins[i], h0, hplus, hmoins, 0);
    y0moins = computeY0(options->projectionLimiterId, phi[i], phiplus[i], phimoins[i], h0, hplus, hmoins, 1);

    // calcul des points d'intersections xd,xg
    xg = computexgxd(phi[i], phiplus[i], phimoins[i], h0, y0plus, y0moins, 0);
    xd = computexgxd(phi[i], phiplus[i], phimoins[i], h0, y0plus, y0moins, 1);

		    
    // calcul des valeurs sur ces points d'intersections
    //yg = computeygyd(phi[i], phiplus[i], phimoins[i], h0, gradphi[i], 0);
    //yd = computeygyd(phi[i], phiplus[i], phimoins[i], h0, gradphi[i], 1);
    
    yg = computeygyd(phi[i], phiplus[i], phimoins[i], h0, y0plus, y0moins, gradphi[i], 0);
    yd = computeygyd(phi[i], phiplus[i], phimoins[i], h0, y0plus, y0moins, gradphi[i], 1);
		    
    // if (cell == -1009) {
    //   std::cout << " -------------" << std::endl;
    //   std::cout << " type = " << type << " I =" << i << " du mat 1"<< std::endl;
    //   std::cout << " y0moins " << y0moins << std::endl;
    //   std::cout << " y0plus " << y0plus << std::endl;
    //   std::cout << " phimoins " << phimoins[i] << std::endl;
    //   std::cout << " phi " << phi[i] << std::endl;
    //   std::cout << " phiplus " << phiplus[i] << std::endl;		      
    //   std::cout << " h0 " << h0 << " donc " << h0/2. << std::endl;		      
    //   std::cout << " hmoins " << hmoins << std::endl;		      
    //   std::cout << " hplus " << hplus << std::endl;		      
    //   std::cout << " grady " << gradphi[i] << std::endl;
    //   std::cout << " xg " << xg << std::endl;
    //   std::cout << " xd " << xd << std::endl;
    //   std::cout << " yg " << yg << std::endl;
    //   std::cout << " yd " << yd << std::endl;
    //   std::cout << " partie_positive_v " << partie_positive_v << std::endl;
    //   cas_PP = 1;
    //  }
		    
		      
    if (type == 0) // flux arriere ou en dessous de cCells, integration entre -h0/2. et -h0/2.+abs(face_normal_velocity)*deltat_n
      {
	// Flux1m : integrale -inf,  -h0/2.+partie_positive_v
	flux1m = INT2Y( -h0/2.+partie_positive_v, -h0/2., phimoins[i], xg, yg);		 
	// Flux1m : integrale -inf,  -h0/2.
	flux1  = INT2Y( -h0/2., -h0/2., phimoins[i], xg, yg);
	//
	// if (abs(flux1-flux1m) > flux_threhold && abs(xg) != abs(xd)) {
	//   std::cout << i << " Activation Plateau Pente  AR Flux -h0/2, xg non nul " << cell << " = " << abs(flux1-flux1m) << std::endl;
		   
	//   std::cout << " xg " << xg << " et h0/2 " << h0/2. << std::endl;
	//   std::cout << " xd " << xd << " et h0/2 " << h0/2. << std::endl;
	//   cas_PP = 1;
	//  }
	//
	// Flux2m : integrale -inf,  -h0/2.+partie_positive_v
	flux2m = INT2Y( -h0/2.+partie_positive_v, xg, yg, xd, yd);
	// Flux2 : integrale -inf,  -h0/2.
	flux2  = INT2Y( -h0/2., xg, yg, xd, yd);		 
	//
	// Flux3m : integrale -inf,  -h0/2.+partie_positive_v
	flux3m = INT2Y( -h0/2.+partie_positive_v, xd, yd, h0/2., phiplus[i]);
	// Flux3 : integrale -inf,  -h0/2.
	flux3  = INT2Y( -h0/2., xd, yd, h0/2., phiplus[i]);
	//
	// if (abs(flux3-flux3m) > flux_threhold && abs(xg) != abs(xd)) {
	//   std::cout << i << " Activation Plateau Pente  AR Flux xd , h0/2, non nul " << cell << " = " << abs(flux3-flux3m) << std::endl;
	//   std::cout << " xg " << xg << " et h0/2 " << h0/2. << std::endl;
	//   std::cout << " xd " << xd << " et h0/2 " << h0/2. << std::endl;
	//   cas_PP = 1;
	// }
	//
	// integrale positive
	Flux[i] = MathFunctions::max(((flux1m-flux1)+(flux2m-flux2)+(flux3m-flux3)), 0.);
	// formule 16
	if (((phiplus[i]-phi[i])*(phimoins[i]-phi[i])) >= 0.) Flux[i] = phi[i]*partie_positive_v;
	//
	//if (cas_PP == 1) std::cout << " type 1 Flux1 " << flux1 << " Flux1m " << flux1m << " Flux2 " << flux2
	//			    << " Flux2m " << flux2m << " Flux3 " << flux3 << " Flux3m " << flux3m << " soit " << Flux[i] << std::endl;
      } else if (type == 1) // flux devant ou au dessus de cCells, integration entre h0/2.-abs(face_normal_velocity)*deltat_n et h0/2.
      {
		 
	// Flux1 : integrale -inf,  h0/2.-partie_positive_v
	flux1  = INT2Y( h0/2.-partie_positive_v, -h0/2., phimoins[i], xg, yg);
	// Flux1m : integrale -inf,  -h0/2.
	flux1m = INT2Y( h0/2., -h0/2., phimoins[i], xg, yg);
	//
	// if (abs(flux1-flux1m) > flux_threhold && abs(xg) != abs(xd)) {
	//   std::cout << i << " Activation Plateau Pente  AV Flux -h0/2, xg non nul " << cell << " = " << abs(flux1-flux1m) << std::endl;
	//   std::cout << " xg " << xg << " et h0/2 " << h0/2. << std::endl;
	//   std::cout << " xd " << xd << " et h0/2 " << h0/2. << std::endl;
	//   cas_PP = 1;
	//  }
	//
	// Flux2 : integrale -inf,  h0/2.-partie_positive_v
	flux2  = INT2Y( h0/2.-partie_positive_v, xg, yg, xd, yd);
	// Flux2m : integrale -inf,  -h0/2.
	flux2m = INT2Y( h0/2., xg, yg, xd, yd);
	//
	// Flux3 : integrale -inf,  h0/2.-partie_positive_v
	flux3  = INT2Y( h0/2.-partie_positive_v, xd, yd, h0/2., phiplus[i]);
	// Flux3m : integrale -inf,  -h0/2.
	flux3m = INT2Y( h0/2., xd, yd, h0/2., phiplus[i]);
	// if (abs(flux3-flux3m) > flux_threhold && abs(xg) != abs(xd)) {
	//    std::cout << i << " Activation Plateau Pente  AV Flux xd , h0/2, non nul " << cell << " = " << abs(flux3-flux3m) << std::endl;
	//  std::cout << " xg " << xg << " et h0/2 " << h0/2. << std::endl;
	//  std::cout << " xd " << xd << " et h0/2 " << h0/2. << std::endl;
	//  cas_PP = 1;
	// }
	//  
	// integrale positive
	Flux[i] = MathFunctions::max(((flux1m-flux1)+(flux2m-flux2)+(flux3m-flux3)), 0.);
	// formule 16
	if (((phiplus[i]-phi[i])*(phimoins[i]-phi[i])) >= 0.) Flux[i] = phi[i]*partie_positive_v;
		 
	//if (cas_PP == 1) std::cout << " type 1 Flux1 " << flux1 << " Flux1m " << flux1m << " Flux2 " << flux2
	//			    << " Flux2m " << flux2m << " Flux3 " << flux3 << " Flux3m " << flux3m << " soit " << Flux[i] << std::endl;
      }
		    
    // Calcul du flux a deplacer peut etre dans computeUremap1 et 2
    // FluxFace[i] = computeflux(faceNormalVelocity, faceNormal, faceLength, deltat_n, h0, phiplus[i], phimoins[i], xd[i], xg[i], yd[i], yg[i], outerFaceNormal, exy);
  }
  //std::cout << " Flux " << Flux << std::endl;
  // les flux de masse, de quantité de mouvement et d'energie massique se deduisent des flux de volumes
  double somme_flux_masse=0.;
  for (imat=0; imat < nbmatmax; imat++) {
    Flux[nbmatmax+imat]     =  phi[nbmatmax+imat]*Flux[imat]; // flux de masse de imat
    Flux[2*nbmatmax+2+imat] =  phi[2*nbmatmax+2+imat]*Flux[nbmatmax+imat]; // flux de masse energy de imat
    somme_flux_masse += Flux[nbmatmax+imat];
  }
  Flux[2*nbmatmax] =  phi[2*nbmatmax]*somme_flux_masse; // flux de quantité de mouvement x
  Flux[2*nbmatmax+1] =  phi[2*nbmatmax+1]*somme_flux_masse; // flux de quantité de mouvement y
  Flux[3*nbmatmax+2] =  phi[3*nbmatmax+2]*somme_flux_masse; // flux d'energie cinetique
   
  return Flux;
}
template<size_t d>
KOKKOS_INLINE_FUNCTION
RealArray1D<d> computeIntersectionPPPure(RealArray1D<d> gradphi,
				     RealArray1D<d> phi, RealArray1D<d> phiplus, RealArray1D<d> phimoins, double h0, double hplus, double hmoins,
				     double face_normal_velocity, double deltat_n, int type, int cell, double flux_threhold) 
{
  RealArray1D<d> Flux = options->Uzero;
  double y0plus, y0moins, xd, xg, yd, yg;
  double flux1, flux2, flux3, flux1m, flux2m, flux3m;	  
  double partie_positive_v = 0.5*(face_normal_velocity + abs(face_normal_velocity))*deltat_n;
  if (partie_positive_v == 0.) return Flux;
  int cas_PP =0;
  // on ne fait que la projection des volumes et masses 
  for (size_t i=0; i<2*nbmatmax; i++) {
    // calcul des seuils y0plus, y0moins pour cCells
    y0plus  = computeY0(options->projectionLimiterIdPure, phi[i], phiplus[i], phimoins[i], h0, hplus, hmoins, 0);
    y0moins = computeY0(options->projectionLimiterIdPure, phi[i], phiplus[i], phimoins[i], h0, hplus, hmoins, 1);

    // calcul des points d'intersections xd,xg
    xg = computexgxd(phi[i], phiplus[i], phimoins[i], h0, y0plus, y0moins, 0);
    xd = computexgxd(phi[i], phiplus[i], phimoins[i], h0, y0plus, y0moins, 1);

		    
    // calcul des valeurs sur ces points d'intersections
    //yg = computeygyd(phi[i], phiplus[i], phimoins[i], h0, gradphi[i], 0);
    //yd = computeygyd(phi[i], phiplus[i], phimoins[i], h0, gradphi[i], 1);
    
    yg = computeygyd(phi[i], phiplus[i], phimoins[i], h0, y0plus, y0moins, gradphi[i], 0);
    yd = computeygyd(phi[i], phiplus[i], phimoins[i], h0, y0plus, y0moins, gradphi[i], 1);
		    
    // if (cell == -1009) {
    //   std::cout << " -------------" << std::endl;
    //   std::cout << " type = " << type << " I =" << i << " du mat 1"<< std::endl;
    //   std::cout << " y0moins " << y0moins << std::endl;
    //   std::cout << " y0plus " << y0plus << std::endl;
    //   std::cout << " phimoins " << phimoins[i] << std::endl;
    //   std::cout << " phi " << phi[i] << std::endl;
    //   std::cout << " phiplus " << phiplus[i] << std::endl;		      
    //   std::cout << " h0 " << h0 << " donc " << h0/2. << std::endl;		      
    //   std::cout << " hmoins " << hmoins << std::endl;		      
    //   std::cout << " hplus " << hplus << std::endl;		      
    //   std::cout << " grady " << gradphi[i] << std::endl;
    //   std::cout << " xg " << xg << std::endl;
    //   std::cout << " xd " << xd << std::endl;
    //   std::cout << " yg " << yg << std::endl;
    //   std::cout << " yd " << yd << std::endl;
    //   std::cout << " partie_positive_v " << partie_positive_v << std::endl;
    //   cas_PP = 1;
    //  }
		    
		      
    if (type == 0) // flux arriere ou en dessous de cCells, integration entre -h0/2. et -h0/2.+abs(face_normal_velocity)*deltat_n
      {
	// Flux1m : integrale -inf,  -h0/2.+partie_positive_v
	flux1m = INT2Y( -h0/2.+partie_positive_v, -h0/2., phimoins[i], xg, yg);		 
	// Flux1m : integrale -inf,  -h0/2.
	flux1  = INT2Y( -h0/2., -h0/2., phimoins[i], xg, yg);
	//
	// if (abs(flux1-flux1m) > flux_threhold && abs(xg) != abs(xd)) {
	//   std::cout << i << " Activation Plateau Pente  AR Flux -h0/2, xg non nul " << cell << " = " << abs(flux1-flux1m) << std::endl;
		   
	//   std::cout << " xg " << xg << " et h0/2 " << h0/2. << std::endl;
	//   std::cout << " xd " << xd << " et h0/2 " << h0/2. << std::endl;
	//   cas_PP = 1;
	//  }
	//
	// Flux2m : integrale -inf,  -h0/2.+partie_positive_v
	flux2m = INT2Y( -h0/2.+partie_positive_v, xg, yg, xd, yd);
	// Flux2 : integrale -inf,  -h0/2.
	flux2  = INT2Y( -h0/2., xg, yg, xd, yd);		 
	//
	// Flux3m : integrale -inf,  -h0/2.+partie_positive_v
	flux3m = INT2Y( -h0/2.+partie_positive_v, xd, yd, h0/2., phiplus[i]);
	// Flux3 : integrale -inf,  -h0/2.
	flux3  = INT2Y( -h0/2., xd, yd, h0/2., phiplus[i]);
	//
	// if (abs(flux3-flux3m) > flux_threhold && abs(xg) != abs(xd)) {
	//   std::cout << i << " Activation Plateau Pente  AR Flux xd , h0/2, non nul " << cell << " = " << abs(flux3-flux3m) << std::endl;
	//   std::cout << " xg " << xg << " et h0/2 " << h0/2. << std::endl;
	//   std::cout << " xd " << xd << " et h0/2 " << h0/2. << std::endl;
	//   cas_PP = 1;
	// }
	//
	// integrale positive
	Flux[i] = MathFunctions::max(((flux1m-flux1)+(flux2m-flux2)+(flux3m-flux3)), 0.);
	// formule 16
	if (((phiplus[i]-phi[i])*(phimoins[i]-phi[i])) >= 0.) Flux[i] = phi[i]*partie_positive_v;
	//
	//if (cas_PP == 1) std::cout << " type 1 Flux1 " << flux1 << " Flux1m " << flux1m << " Flux2 " << flux2
	//			    << " Flux2m " << flux2m << " Flux3 " << flux3 << " Flux3m " << flux3m << " soit " << Flux[i] << std::endl;
      } else if (type == 1) // flux devant ou au dessus de cCells, integration entre h0/2.-abs(face_normal_velocity)*deltat_n et h0/2.
      {
		 
	// Flux1 : integrale -inf,  h0/2.-partie_positive_v
	flux1  = INT2Y( h0/2.-partie_positive_v, -h0/2., phimoins[i], xg, yg);
	// Flux1m : integrale -inf,  -h0/2.
	flux1m = INT2Y( h0/2., -h0/2., phimoins[i], xg, yg);
	//
	// if (abs(flux1-flux1m) > flux_threhold && abs(xg) != abs(xd)) {
	//   std::cout << i << " Activation Plateau Pente  AV Flux -h0/2, xg non nul " << cell << " = " << abs(flux1-flux1m) << std::endl;
	//   std::cout << " xg " << xg << " et h0/2 " << h0/2. << std::endl;
	//   std::cout << " xd " << xd << " et h0/2 " << h0/2. << std::endl;
	//   cas_PP = 1;
	//  }
	//
	// Flux2 : integrale -inf,  h0/2.-partie_positive_v
	flux2  = INT2Y( h0/2.-partie_positive_v, xg, yg, xd, yd);
	// Flux2m : integrale -inf,  -h0/2.
	flux2m = INT2Y( h0/2., xg, yg, xd, yd);
	//
	// Flux3 : integrale -inf,  h0/2.-partie_positive_v
	flux3  = INT2Y( h0/2.-partie_positive_v, xd, yd, h0/2., phiplus[i]);
	// Flux3m : integrale -inf,  -h0/2.
	flux3m = INT2Y( h0/2., xd, yd, h0/2., phiplus[i]);
	// if (abs(flux3-flux3m) > flux_threhold && abs(xg) != abs(xd)) {
	//    std::cout << i << " Activation Plateau Pente  AV Flux xd , h0/2, non nul " << cell << " = " << abs(flux3-flux3m) << std::endl;
	//  std::cout << " xg " << xg << " et h0/2 " << h0/2. << std::endl;
	//  std::cout << " xd " << xd << " et h0/2 " << h0/2. << std::endl;
	//  cas_PP = 1;
	// }
	//  
	// integrale positive
	Flux[i] = MathFunctions::max(((flux1m-flux1)+(flux2m-flux2)+(flux3m-flux3)), 0.);
	// formule 16
	if (((phiplus[i]-phi[i])*(phimoins[i]-phi[i])) >= 0.) Flux[i] = phi[i]*partie_positive_v;
		 
	//if (cas_PP == 1) std::cout << " type 1 Flux1 " << flux1 << " Flux1m " << flux1m << " Flux2 " << flux2
	//			    << " Flux2m " << flux2m << " Flux3 " << flux3 << " Flux3m " << flux3m << " soit " << Flux[i] << std::endl;
      }
		    
    // Calcul du flux a deplacer peut etre dans computeUremap1 et 2
    // FluxFace[i] = computeflux(faceNormalVelocity, faceNormal, faceLength, deltat_n, h0, phiplus[i], phimoins[i], xd[i], xg[i], yd[i], yg[i], outerFaceNormal, exy);
  }
  //std::cout << " Flux " << Flux << std::endl;
  // les flux de masse, de quantité de mouvement et d'energie massique se deduisent des flux de volumes
  double somme_flux_masse=0.;
  for (imat=0; imat < nbmatmax; imat++) {
    Flux[2*nbmatmax+2+imat] =  phi[2*nbmatmax+2+imat]*Flux[nbmatmax+imat]; // flux de masse energy de imat
    somme_flux_masse += Flux[nbmatmax+imat];
  }
  Flux[2*nbmatmax] =  phi[2*nbmatmax]*somme_flux_masse; // flux de quantité de mouvement x
  Flux[2*nbmatmax+1] =  phi[2*nbmatmax+1]*somme_flux_masse; // flux de quantité de mouvement y
  Flux[3*nbmatmax+2] =  phi[3*nbmatmax+2]*somme_flux_masse; // flux d'energie cinetique
   
  return Flux;
}
template<size_t d>
KOKKOS_INLINE_FUNCTION
RealArray1D<d> computeUpwindFaceQuantities(RealArray1D<dim> face_normal, double face_normal_velocity, double delta_x, RealArray1D<dim> x_f, RealArray1D<d> phi_cb, RealArray1D<d> grad_phi_cb, RealArray1D<dim> x_cb, RealArray1D<d> phi_cf, RealArray1D<d> grad_phi_cf, RealArray1D<dim> x_cf) 
{
  if (face_normal_velocity * delta_x > 0.0) 
    return ArrayOperations::plus(phi_cb, ArrayOperations::multiply(MathFunctions::dot(ArrayOperations::minus(x_f, x_cb), face_normal), grad_phi_cb));
  else 
    return ArrayOperations::plus(phi_cf, ArrayOperations::multiply(MathFunctions::dot(ArrayOperations::minus(x_f, x_cf), face_normal), grad_phi_cf));
}

KOKKOS_INLINE_FUNCTION
RealArray1D<dim> xThenYToDirection(bool x_then_y_) 
{
  if (x_then_y_) 
    return {{1.0, 0.0}};



  else 
    return {{0.0, 1.0}};
}
template<size_t d>
KOKKOS_INLINE_FUNCTION
RealArray1D<d> computeRemapFlux(int projectionAvecPlateauPente, double face_normal_velocity, RealArray1D<dim> face_normal, double face_length, RealArray1D<d> phi_face, RealArray1D<dim> outer_face_normal, RealArray1D<dim> exy, double deltat_n) 
{
  RealArray1D<d> Flux;
  if (projectionAvecPlateauPente == 0) 
    {
      if (MathFunctions::fabs(MathFunctions::dot(face_normal, exy)) < 1.0E-10) 
	return ArrayOperations::multiply(0.0, phi_face);
      return ArrayOperations::multiply(MathFunctions::dot(outer_face_normal, exy) * face_normal_velocity * face_length * deltat_n, phi_face);
    } else 	{
    if (MathFunctions::fabs(MathFunctions::dot(face_normal, exy)) < 1.0E-10) 
      return ArrayOperations::multiply(0.0, phi_face);
    return ArrayOperations::multiply(MathFunctions::dot(outer_face_normal, exy) * face_length, phi_face);
  }
		

}
