//saralamba@gmail.com

#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

  
SEXP getVar(std::string varname){
  Environment glob = Environment::global_env();
  SEXP var = glob[varname];
  return var;
}


Function getFunc(std::string funcname){
  Environment glob = Environment::global_env();
  Function fnc = glob[funcname];
  return fnc;
}


// [[Rcpp::export]]
List measlesmod(double t, arma::vec y, List parms){

  Function tfun = getFunc("tfun");
  Function tfun2 = getFunc("tfun2");
  Function tfun3 = getFunc("tfun3");
  Function birthfunc2 = getFunc("birth.func2");
  
  
  arma::vec S = y.subvec(0,100);
  arma::vec I = y.subvec(101,201);
  arma::vec R = y.subvec(202,302);
  arma::vec V1S = y.subvec(303,403);
  arma::vec V1I = y.subvec(404,504);
  arma::vec V1R = y.subvec(505,605);
  arma::vec V2S = y.subvec(606,706);
  arma::vec V2I = y.subvec(707,807);
  arma::vec V2R = y.subvec(808,908);
  arma::vec V3S = y.subvec(909,1009);
  arma::vec V3I = y.subvec(1010,1110);
  arma::vec V3R = y.subvec(1111,1211);
  arma::vec newcaseindex = y.subvec(1212,1312);
  
    
  //I[I<0] = 0;
  I.elem(find(I<0)).fill(0);
  //R[R<0] = 0;
  R.elem(find(R<0)).fill(0);
  //V1I[V1I<0] = 0;
  V1I.elem(find(V1I<0)).fill(0);
  //V1R[V1R<0] = 0;
  V1R.elem(find(V1R<0)).fill(0);
  //V2I[V2I<0] = 0;
  V2I.elem(find(V2I<0)).fill(0);
  //V2R[V2R<0] = 0;
  V2R.elem(find(V2R<0)).fill(0);
  //V3I[V3I<0] = 0;
  V3I.elem(find(V3I<0)).fill(0);
  //V3R[V3R<0] = 0;
  V3R.elem(find(V3R<0)).fill(0);
  
  List ff = getVar("ff");
  arma::vec mort = as<arma::vec>(ff["m1970"]);
  
  double N = arma::sum(S+I+R+V1S+V1I+V1R+V2S+V2I+V2R+V3S+V3I+V3R);
  
  double inf = as<double>(getVar("inf"));
  //#infectious persons
  arma::vec II = inf*(I+V1I+V2I+V3I);
    
  double newborn = as<double>(birthfunc2(t));
  newborn = (newborn/365)*N;
  
  //NumericMatrix contact = getVar("contact");
  arma::mat contact = as<arma::mat>(getVar("contact"));
  
  arma::mat lambda = contact * II/N; 
  
  //#coverage of measles immunization
    
  double v1cov = as<double>(tfun(1955+t/365));
  v1cov = v1cov/100;
  
  double v2cov = as<double>(tfun2(1955+t/365));
  v2cov = v2cov/100;
  
  //v1cov[Rcpp::is_na(v1cov)] = 0;
  //v1cov.replace(arma::datum::nan,0);
  //v2cov[Rcpp::is_na(v2cov)] = 0;
  //v2cov.replace(arma::datum::nan,0);
  if(R_IsNA(v1cov)) v1cov=0;
  if(R_IsNA(v2cov)) v2cov=0;
  
  double v3cov = as<double>(tfun3(1955+t/365));
  v3cov = v3cov/100;
  
  //v3cov[Rcpp::is_na(v3cov)] = 0;
  //v3cov.replace(arma::datum::nan,0);
  if(R_IsNA(v3cov)) v3cov=0;
  
  
  // NumericVector v1doset = NumericVector::create(0.,-log(1-v1cov)/30,rep(0.0,99));
  // NumericVector v2doset = NumericVector::create(0.,0.,0.,0.,-log(1.-v2cov)/30,rep(0.0,96));
  // NumericVector v3doset = NumericVector::create(rep(-log(1-v2.cov)/30,40),rep(0.0,61));

  arma::vec v1dose = {0,-log(1-v1cov)/30,0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 
                      0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 
                      0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 
                      0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 
                      0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 
                      0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };
  arma::vec v2dose = {0,0,0,0,-log(1.-v2cov)/30,0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 
                      0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 
                      0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 
                      0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 
                      0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 
                      0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
  arma::vec v3dose = {-log(1-v2cov)/30,-log(1-v2cov)/30,-log(1-v2cov)/30,-log(1-v2cov)/30,-log(1-v2cov)/30,-log(1-v2cov)/30,-log(1-v2cov)/30,
                      -log(1-v2cov)/30,-log(1-v2cov)/30,-log(1-v2cov)/30,-log(1-v2cov)/30,-log(1-v2cov)/30,-log(1-v2cov)/30,-log(1-v2cov)/30,
                      -log(1-v2cov)/30,-log(1-v2cov)/30,-log(1-v2cov)/30,-log(1-v2cov)/30,-log(1-v2cov)/30,-log(1-v2cov)/30,-log(1-v2cov)/30,
                      -log(1-v2cov)/30,-log(1-v2cov)/30,-log(1-v2cov)/30,-log(1-v2cov)/30,-log(1-v2cov)/30,-log(1-v2cov)/30,-log(1-v2cov)/30,
                      -log(1-v2cov)/30,-log(1-v2cov)/30,-log(1-v2cov)/30,-log(1-v2cov)/30,-log(1-v2cov)/30,-log(1-v2cov)/30,-log(1-v2cov)/30,
                      -log(1-v2cov)/30,-log(1-v2cov)/30,-log(1-v2cov)/30,-log(1-v2cov)/30,-log(1-v2cov)/30, 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                      0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                      0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
  
  
  arma::mat aging =  as<arma::mat>(getVar("aging"));
  
  double w1 = as<double>(getVar("w1"));
  double w2 = as<double>(getVar("w2"));
  double w3 = as<double>(getVar("w3"));
  
  double eff1 = as<double>(getVar("eff1"));
  double eff2 = as<double>(getVar("eff2"));
  double eff3 = as<double>(getVar("eff3"));
  
  double gamma = as<double>(getVar("gamma"));
  
  arma::vec dS = -lambda%S-mort%S+aging*S-v1dose%S+w1*(V1S+V1I+V1R)+w2*(V2S+V2I+V2R)-v3dose%S+w3*(V3S+V3I+V3R);
  arma::vec dI = lambda%S-mort%I+aging*I-gamma*I-v1dose%I-v3dose%I;
  arma::vec dR = gamma*I-mort%R+aging*R-v1dose%R-v3dose%R;
  arma::vec dV1S = -lambda%V1S-mort%V1S+aging*V1S+(1-eff1)*v1dose%S-v2dose%V1S-w1*V1S-v3dose%V1S;
  arma::vec dV1I = lambda%V1S-mort%V1I+aging*V1I-gamma*V1I+(1-eff1)*v1dose%I-v2dose%V1I-w1*V1I-v3dose%V1I;
  arma::vec dV1R = gamma*V1I-mort%V1R+aging*V1R+eff1*v1dose%(S+I)+v1dose%R-v2dose%V1R-w1*V1R-v3dose%V1R;
  arma::vec dV2S = -lambda%V2S-mort%V2S+aging*V2S+(1-eff2)*v2dose%V1S-w2*V2S-v3dose%V2S;
  arma::vec dV2I = lambda%V2S-mort%V2I+aging*V2I-gamma*V2I+(1-eff2)*v2dose%V1I-v2dose%V1I-w2*V2I-v3dose%V2I;
  arma::vec dV2R = gamma*V2I-mort%V2R+aging*V2R+eff2*v2dose%(V1S+V1I)+v2dose%V1R-w2*V2R-v3dose%V2R;
  arma::vec dV3S = -lambda%V3S-mort%V3S+aging*V3S+(1-eff3)*v3dose%(V1S+V2S+S)-w3*V3S;
  arma::vec dV3I = lambda%V3S-mort%V3I+aging*V3I-gamma*V3I+(1-eff3)*v3dose%(V1I+V2I+I)-w3*V3I;
  arma::vec dV3R = gamma*V3I-mort%V3R+aging*V3R+eff3*v3dose%(V1S+V1I+V2S+V2I+S+I)+v3dose%(V1R+V2R+R)-w3*V3R;
  
    
  arma::vec dnewcase = lambda%(S+V1S+V2S);
  dS[0] = dS[0]+newborn;
  
  List output(5);
  arma::vec outvec = dS;
  outvec.insert_rows(dI.size(),dI);
  outvec.insert_rows(dR.size(),dR);
  outvec.insert_rows(dV1S.size(),dV1S);
  outvec.insert_rows(dV1I.size(),dV1I);
  outvec.insert_rows(dV1R.size(),dV1R);
  
  outvec.insert_rows(dV2S.size(),dV2S);
  outvec.insert_rows(dV2I.size(),dV2I);
  outvec.insert_rows(dV2R.size(),dV2R);
  
  outvec.insert_rows(dV3S.size(),dV3S);
  outvec.insert_rows(dV3I.size(),dV3I);
  outvec.insert_rows(dV3R.size(),dV3R);
  
  outvec.insert_rows(dnewcase.size(),dnewcase);
  
  //output[""] = NumericVector::create(dS,dI,dR,dV1S,dV1I,dV1R,dV2S,dV2I,dV2R,dV3S,dV3I,dV3R,dnewcase);
  output[0] = Rcpp::NumericVector(outvec.begin(),outvec.end());
  output[1]=v1cov;
  output[2]=v2cov;
  output[3]=v3cov;
  output[4]=N;

    //list(c(dS,dI,dR,dV1S,dV1I,dV1R,dV2S,dV2I,dV2R,dV3S,dV3I,dV3R,dnewcase),v1.cov,v2.cov,v3.cov,N)
  return output;

}
  
  

/*** R


*/
