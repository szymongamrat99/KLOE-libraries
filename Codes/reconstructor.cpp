#include "reconstructor.h"
#include "../const.h"
#include <cmath>

const double c = c_vel;

void Reconstructor::SetClu(int i, double x, double y, double z, double t, double E){
  _clu[i][0] = x;
  _clu[i][1] = y;
  _clu[i][2] = z;
  _clu[i][3] = t;
  _ene[i] = E;
}

/*
  wartosci w selected liczone od 1
*/
Solution Reconstructor::MySolve(int * selected){
  
  double t[4];
  double x[4];
  double y[4];
  double z[4];
  
  for(int i=0;i<4;i++){
    x[i] = _clu[ selected[i] - 1 ][0];
    y[i] = _clu[ selected[i] - 1 ][1];
    z[i] = _clu[ selected[i] - 1 ][2];
    t[i] = _clu[ selected[i] - 1 ][3];
  }
  
  // tablice robocze
  double T[4];
  double B[4];
  double C[4];
  double D[4];
  double X[4];
  double Y[4];
  double Z[4];
  
  // form linear system
  for(int k=1;k<4;k++){
    T[k] = 2.*c*c*(t[k]-t[0]);
    X[k] = 2.*(x[0]-x[k]);
    Y[k] = 2.*(y[0]-y[k]);
    Z[k] = 2.*(z[0]-z[k]);
    
    B[k] = c*c*(t[k]*t[k]-t[0]*t[0]) + x[0]*x[0]-x[k]*x[k] + y[0]*y[0]-y[k]*y[k] + z[0]*z[0]-z[k]*z[k];
  }
  
  // solution of linear system
  double denom = -1.*X[3]*Y[2]*Z[1]+X[2]*Y[3]*Z[1]+X[3]*Y[1]*Z[2]-X[1]*Y[3]*Z[2]-X[2]*Y[1]*Z[3]+X[1]*Y[2]*Z[3];

  // poszukwanie zrodla NaN
  //assert( denom != 0 );
  
  C[1] = ( B[2]*Y[3]*Z[1]-B[1]*Y[3]*Z[2]+B[3]*(-Y[2]*Z[1]+Y[1]*Z[2])-B[2]*Y[1]*Z[3]+B[1]*Y[2]*Z[3]) / denom; 
  
  C[2] = (-1.*B[2]*X[3]*Z[1]+B[1]*X[3]*Z[2]+B[3]*(X[2]*Z[1]-X[1]*Z[2])+B[2]*X[1]*Z[3]-B[1]*X[2]*Z[3]) / denom; 
  
  C[3] = ( B[2]*X[3]*Y[1]-B[1]*X[3]*Y[2]+B[3]*(-X[2]*Y[1]+X[1]*Y[2])-B[2]*X[1]*Y[3]+B[1]*X[2]*Y[3]) / denom; 
  
  D[1] = ( T[3]*Y[2]*Z[1]-T[2]*Y[3]*Z[1]-T[3]*Y[1]*Z[2]+T[1]*Y[3]*Z[2]+T[2]*Y[1]*Z[3]-T[1]*Y[2]*Z[3] ) / denom;
  
  D[2] =  ( -1.*T[3]*X[2]*Z[1]+T[2]*X[3]*Z[1]+T[3]*X[1]*Z[2]-T[1]*X[3]*Z[2]-T[2]*X[1]*Z[3]+T[1]*X[2]*Z[3] ) / denom;
  
  D[3] = ( T[3]*X[2]*Y[1]-T[2]*X[3]*Y[1]-T[3]*X[1]*Y[2]+T[1]*X[3]*Y[2]+T[2]*X[1]*Y[3]-T[1]*X[2]*Y[3] ) / denom;
  
  // solutions of quadratic equation
  Solution S;
  
  denom = 2.*( c*c - D[1]*D[1] - D[2]*D[2] - D[3]*D[3] );
  
  double p = 2.*( C[1]*D[1] + C[2]*D[2] +C[3]*D[3] + c*c*t[0] - D[1]*x[0] - D[2]*y[0] - D[3]*z[0] ) ;
  double q = -C[1]*C[1]-C[2]*C[2]-C[3]*C[3]+c*c*t[0]*t[0]+2.*C[1]*x[0]-x[0]*x[0]+2.*C[2]*y[0]-y[0]*y[0]+2.*C[3]*z[0]-z[0]*z[0];
  
  //double delta;
  //delta = p*p - 2.*denom*q;
  
  S.sol[0][3] = ( p - sqrt( p*p - 2.*denom*q ) ) / denom;
  S.sol[1][3] = ( p + sqrt( p*p - 2.*denom*q ) ) / denom;

  for(int i=0;i<2;i++){
    S.error[i] = false;
    for(int j=0;j<3;j++){
      S.sol[i][j] = C[j+1]+D[j+1]*S.sol[i][3];
      if( std::isnan(  S.sol[i][j] ) ||
	  std::isfinite( S.sol[i][j] ) == 0 ){
	S.error[i] = true;
      }
    }
  }

  return S;

}

/**************************** Kleusberg *****************************/
Solution Reconstructor::KleusbergSolve(int * selected){

  Solution S;
  return S;
}

/*
  x0 - initial guess (x,y,z,t)
 */
Solution Reconstructor::LeastSquaresSolve(double * x0){

  Solution S;
  return S;
}

/************************ Solve with Minui5Bt *************************/
Reconstructor * r;

double rsq(const double * x){

  double sum = 0;

  for(int i=0;i<6;i++){
    sum += pow( r->ResidualErr(i,x) , 2.);
  }
  return sum;
}

Solution Reconstructor::MinuitSolve(double * x0){
  
  Solution S;
  return S;
}

/******************** Residual Error computation ********************/
// i - counts from 0
double Reconstructor::ResidualErr(int i, const double * x){
  double err = 0;
  err = sqrt(pow(_clu[i][0] - x[0], 2.) + 
	     pow(_clu[i][1] - x[1], 2.) + 
	     pow(_clu[i][2] - x[2], 2.) ) -
    c*( _clu[i][3] - x[3] );
  
  return fabs(err);
}

double Reconstructor::ResidualErrTot(const double * x){
  double R =  0;
  for(int i=0; i<6;i++){
    R += ResidualErr( i, x);
  }
  return R;
}

/************** Total energy of 4-cluster combination ***************/
double Reconstructor::CombEnergy(int * selected){
  
  double E = 0;
  for(int i=0;i<4;i++){
    E += _ene[ selected[i]-1 ];
  }
  return E;
}

/****************** Total energy of all 6 clusters ******************/
double Reconstructor::TotalEnergy()const{

  double E = 0;
  for(int i=0;i<6;i++){
    E += _ene[ i ];
  }
  return E;

}

/****************** Invariant masses of three gg pairs **************/
/*
  Fills Minvgg (double[3]) with gg invariant masses
  Returns pi0pi0pi0 (or 6gamma) invariant mass
 */
double Reconstructor::GetInvMasses(double * sol, int * comb, 
				   double * Minvgg)const{
  
  // calculate momenta of gammas
  double p[6][3];
  double tot;
  for(int i=0;i<6;i++){ // over gammas
    tot = 0;
    for(int j=0;j<3;j++){ // over coordinates
      p[i][j] = _clu[i][j] - sol[j];
      tot += pow( p[i][j] , 2.);
    }
    tot = sqrt( tot );
    for(int j=0;j<3;j++){ // over coordinates
      p[i][j] = p[i][j] * _ene[i] / tot;
    }
  }
  
  // calculate gg invariant masses
  for(int i=0; i<3;i++){
    Minvgg[i] = sqrt( pow( _ene[comb[2*i]-1] + _ene[comb[2*i+1]-1] ,2.) -
		      pow( p[comb[2*i]-1][0] + p[comb[2*i+1]-1][0] ,2.) -
		      pow( p[comb[2*i]-1][1] + p[comb[2*i+1]-1][1] ,2.) -
		      pow( p[comb[2*i]-1][2] + p[comb[2*i+1]-1][2] ,2.) );
  }
  
  // calculate 6g invariant mass
  double M6g = sqrt( pow( _ene[0]+_ene[1]+_ene[2]+_ene[3]+_ene[4]+_ene[5] ,2.) -
		     pow( p[0][0]+p[1][0]+p[2][0]+p[3][0]+p[4][0]+p[5][0] ,2.) -
		     pow( p[0][1]+p[1][1]+p[2][1]+p[3][1]+p[4][1]+p[5][1] ,2.) -
		     pow( p[0][2]+p[1][2]+p[2][2]+p[3][2]+p[4][2]+p[5][2] ,2.) );

  return M6g;
}

/****************** Difference between Mpi and mgg x 3 **************/
double Reconstructor::GetInvMassDiscrepancy(double * sol, int * comb)const{

  const double mPi0 = m_pi0; // Mev/c2  
  double * mgg = new double[3];
  GetInvMasses(sol, comb, mgg);
  
  // calculate sigmas of energy
  double dE[6];
  for(int i=0;i<6;i++){
    dE[i] = 0.057 * _ene[i] / sqrt( _ene[i]/1000. );
  }
  // calculate sigmas of invariant masses
  double dM[3];
  for(int i=0;i<3;i++){
    dM[i] = 0.5*sqrt( (_ene[comb[2*i+1]-1] / _ene[comb[2*i]-1]) * 
		      pow( dE[2*i] , 2.) +
		      (_ene[comb[2*i]-1] / _ene[comb[2*i+1]-1]) * 
		      pow( dE[2*i+1] , 2.) );
  }
  
  // finally calculate discrepancy
  double d = 0;
  for(int i=0;i<3;i++){
    d += abs( mgg[i] - mPi0 ) / dM[i];
    //    d += abs( mgg[i] - mPi0 ) * pow( _ene[comb[2*i+1]-1] + _ene[comb[2*i]-1], -2. );
  }

  return d;
}

/****************** Kaon momentum from gamma momenta   **************/
/*
  fills p[] with K momentum coordinates
  p[0-2] - x,y,z
  p[3] - K energy
  p[4] - total K momentum
  uses solution sol for vertex
 */
void Reconstructor::GetKmomentum(const double * sol, float * p)const{

  const double mK0  = m_k0;  // Mev/c2 

  // calculate momenta of gammas
  double pgam[6][3];
  double tot;
  for(int i=0;i<6;i++){ // over gammas
    tot = 0;
    for(int j=0;j<3;j++){ // over coordinates
      pgam[i][j] = _clu[i][j] - sol[j];
      tot += pow( pgam[i][j] , 2.);
    }
    tot = sqrt( tot );
    for(int j=0;j<3;j++){ // over coordinates
      pgam[i][j] = pgam[i][j] * _ene[i] / tot;
    }
  }
  
  // total momentum of kaon
  p[4] = 0;
  for(int j=0;j<3;j++){
    p[j] = 0;
    for(int i=0;i<6;i++){
      p[j] += pgam[i][j];
    }
    p[4] += p[j]*p[j];
  }
  p[3] = sqrt( p[4] + mK0*mK0 );
  p[4] = sqrt( p[4] );
  
}











































































