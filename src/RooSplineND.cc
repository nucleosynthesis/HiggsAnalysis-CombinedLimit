#include "../interface/RooSplineND.h"
/*
Use of radial basis functions for interpolation 
between points in multidimensional space.

Produces N ->1 function from TTree

Branch to be considered as F(x) should be passed in fName, 
otherwise it is assumed to be called, "f"
*/


RooSplineND::RooSplineND(const char *name, const char *title, RooArgList &vars, TTree *tree, const char *fName, double eps) :
  RooAbsReal(name,title),
  vars_("vars","Variables", this)
{
  ndim_ = vars.getSize();
  M_    = tree->GetEntries();

  // Interface is from ROOT TTree. Assume that each variable in the vars 
  // corresponds to a (double) branch in the tree
  std::map<int,double> b_map;

  RooAbsReal *rIt;	
  TIterator *iter = vars.createIterator(); int it_c=0;
  while( (rIt = (RooAbsReal*) iter->Next()) ){ 
    vars_.add(*rIt);
    std::vector<double >tmpv(M_,0); 
    v_map.insert(std::pair<int, std::vector<double> >(it_c,tmpv));
    b_map.insert(std::pair<int, double >(it_c,0));
    r_map.insert(std::pair<int, std::pair<double,double> >(it_c,std::pair<double,double>(1.e6,-1e6))); 
    tree->SetBranchAddress(rIt->GetName(),&b_map[it_c]);
    it_c++;
  }
  // Assume the function val (yi) branch is f
  double F;
  tree->SetBranchAddress(fName,&F);
  std::vector<double> F_vec;

  // Run through tree and store points 
  for (int i=0;i<M_;i++){
    tree->GetEntry(i);
    //std::cout << "Adding point, " << i << " ";
    for (int k=0;k<ndim_;k++){
      double cval = b_map[k];
      if (cval < r_map[k].first) r_map[k].first=cval;
      if (cval > r_map[k].second) r_map[k].second=cval;
      std::cout << cval << " ";
      v_map[k][i] = cval;
    }
    F_vec.push_back(F);
    std::cout << "  func = " << F  <<std::endl;
  }

  std::cout << "RooSplineND -- Num Dimensions == " << ndim_ <<std::endl;
  std::cout << "RooSplineND -- Num Samples    == " << M_ << std::endl;
  //eps_  = 2*1./M_; // This is a parameter which should be configurable !
  //eps_=0.5;
  
  axis_pts_ = TMath::Power(M_,1./ndim_);
  //eps_ = 0.1*axis_pts_;
  // distances between neighboruing points are ~ 1
  // hence should define an eps based on that. 
  // Best guess for width of Gaussian is to use 1.5 neighbours on each side, ie eps = 3
  // recommendation is 3
  eps_= eps;
  std::cout << "NPTS_AXIS == " << axis_pts_<<std::endl;
  std::cout << "EPS == " << eps_<<std::endl;
  // Init and solve for weights
  calculateWeights(F_vec); 	
}

//_____________________________________________________________________________
// Copy Constructor
RooSplineND::RooSplineND(const RooSplineND& other, const char *name) :
 RooAbsReal(other, name),vars_("vars",this,RooListProxy())
{
  ndim_ = other.ndim_;
  M_    = other.M_;
  eps_  = other.eps_;
  axis_pts_ = other.axis_pts_;

  RooAbsReal *rIt;	
  TIterator *iter = other.vars_.createIterator();
  while( (rIt = (RooAbsReal*) iter->Next()) ){ 
    vars_.add(*rIt);
  }

  // STL copy constructors
  w_    = other.w_;
  v_map = other.v_map; 
  r_map = other.r_map;

}
//_____________________________________________________________________________
// Clone Constructor
RooSplineND::RooSplineND(const char *name, const char *title, const RooListProxy &vars, 
 int ndim, int M, double eps, std::vector<double> &w, std::map<int,std::vector<double> > &map, std::map<int,std::pair<double,double> > &rmap) :
 RooAbsReal(name, title),vars_("vars",this,RooListProxy()) 
{

  RooAbsReal *rIt;	
  TIterator *iter = vars.createIterator();
  while( (rIt = (RooAbsReal*) iter->Next()) ){ 
    vars_.add(*rIt);
  }

  ndim_ = ndim;
  M_    = M;
  eps_  = eps;
  axis_pts_ = TMath::Power(M_,1./ndim_);

  w_    = w;
  v_map = map;
  r_map = rmap;
  /*
  std::map<int,std::vector<double> >::const_iterator vit = map.begin();
  for (;vit!=map.end();vit++){
   v_map.insert(*vit);
  }
  */
  
}

//_____________________________________________________________________________
TObject *RooSplineND::clone(const char *newname) const 
{
    return new RooSplineND(newname, this->GetTitle(), 
	vars_,ndim_,M_,eps_,w_,v_map,r_map);
}
//_____________________________________________________________________________
RooSplineND::~RooSplineND() 
{
}
//_____________________________________________________________________________
TGraph * RooSplineND::getGraph(const char *xvar, double step){

  TGraph *gr = new TGraph();
  gr->SetLineWidth(2);
  RooRealVar* v = (RooRealVar*) vars_.find(xvar);
  gr->SetTitle(v->GetTitle());
  double vorig = v->getVal();
  int cp=0;

  for (double xv=v->getMin();xv<=v->getMax();xv+=step){
    v->setVal(xv);
    gr->SetPoint(cp,xv,evaluate());
    cp++;
  }

  v->setVal(vorig);
  return gr;
}
//_____________________________________________________________________________
void RooSplineND::calculateWeights(std::vector<double> &f){

  // Solve system of Linear equations for weights vector 
  TMatrixTSym<double> fMatrix(M_);
  //TMatrixF fMatrix(M_,M_);
 
  // Fill the Matrix
  for (int i=0;i<M_;i++){
    fMatrix(i,i)=1.;
    for (int j=i+1;j<M_;j++){
        double d2  = getDistSquare(i,j);
	double rad = radialFunc(d2,eps_);
	//std::cout << "Distance between " << i << " and " << j << " = " << TMath::Sqrt(d2) <<std::endl;
        fMatrix(i,j) = (double) rad;
	fMatrix(j,i) = (double) rad; // it is symmetric	
    }
  }

  std::cout << "RooSplineND -- Solving for Weights" << std::endl;
  TVectorD weights(M_);
  for (int i=0;i<M_;i++){
    weights[i]=(double)f[i];
  }

  TDecompBK decomp(fMatrix);
  //decomp.Print();
  // Invert (not sure how stable this is though!)
  decomp.Solve(weights); // Solution now in weights
  //fMatrix.Invert();
  std::cout << "RooSplineND -- ........ Done" << std::endl;
  for (int i=0;i<M_;i++){
    w_.push_back((double)weights[i]);
  }
  // Store weights
  /*
  for (int i=0;i<M_;i++){
    double wi = 0.;
    for (int j=0;j<M_;j++){
  //    std::cout << "i,j " << i << ", " << j << ", " << fMatrix(i,j) <<std::endl;
      wi+=fMatrix(i,j)*f[j];
    }
    //std::cout << "Weights (calc) = " << i << " " << wi <<std::endl;
    w_.push_back(wi);
  }
  */
}
//_____________________________________________________________________________
double RooSplineND::getDistSquare(int i, int j){
  double D = 0.; 
  for (int k=0;k<ndim_;k++){
    double v_i = v_map[k][i];
    double v_j = v_map[k][j];
    double dk = axis_pts_*(v_i-v_j)/(r_map[k].second-r_map[k].first);
    //std::cout << "dimension - " << k << ", Values at pts " << i <<"," << j << " are "<<v_i<< ","<<v_j<< " Distance " <<dk<<std::endl;
    D += dk*dk;
  }
  return D; // only ever use square of distance!
}
//_____________________________________________________________________________
double RooSplineND::getDistFromSquare(int i) const{
  // Read parameters distance from point i in the sample
  double D = 0.; 
  for (int k=0;k<ndim_;k++){
    double v_i = v_map[k][i];
    RooAbsReal *v = (RooAbsReal*)vars_.at(k);
    double v_j = v->getVal();
    double dk = axis_pts_*(v_i-v_j)/(r_map[k].second-r_map[k].first);
    D += dk*dk;
  }
  return D; // only ever use square of distance!
  
}
//_____________________________________________________________________________
double RooSplineND::radialFunc(double d2, double eps) const{
  double expo = (d2/(eps*eps));
  double retval = TMath::Exp(-1*expo);
  //double retval = (eps*eps*eps*eps)/(d2*d2);
//  if (retval < 1e-3) retval=0.;
  return retval;
}
//_____________________________________________________________________________
Double_t RooSplineND::evaluate() const {
 double ret = 0;
 for (int i=0;i<M_;i++){
 //  std::cout << "Weights (eval) = " << i << " " << w_[i] <<std::endl;
 //  std::cout << "EVAL == "<< i << " " << w_[i] << " " << getDistFromSquare(i) << std::endl;
   if (w_[i]==0) continue;
   ret+=(w_[i]*radialFunc(getDistFromSquare(i),eps_));
 }
 return ret;
}
//_____________________________________________________________________________

ClassImp(RooSplineND)
