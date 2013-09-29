#ifndef HiggsAnalysis_CombinedLimit_RooSplineND_h
#define HiggsAnalysis_CombinedLimit_RooSplineND_h

#include <RooAbsReal.h>
#include <RooRealProxy.h>
#include <RooListProxy.h>
#include "TMath.h"
#include "TMatrixTSym.h"
#include "TMatrix.h"
#include "TMatrixF.h"
#include "TMatrixD.h"
#include "TDecompSVD.h"
#include "TVectorD.h"
#include "TDecompLU.h"
#include "TTree.h"
#include "TGraph.h"
#include "RooRealVar.h"

#include <map>
#include <vector>
 
//_________________________________________________
/*
BEGIN_HTML
<p>
RooSplineND is helper class for producing a smooth function 
(F:N->1) given (potentially sparse) samplings in the form of a TTree
</p>
END_HTML
*/
//

class RooSplineND : public RooAbsReal {

   public:
      RooSplineND() : ndim_(0),M_(0),eps_(1.) {}
      RooSplineND(const char *name, const char *title, RooArgList &vars, TTree *tree ) ;
      RooSplineND(const RooSplineND& other, const char *name) ; 
      RooSplineND(const char *name, const char *title, const RooListProxy &vars, int ndim, int M, double eps, std::vector<double> &w, std::map<int,std::vector<double> > &map) ;
      ~RooSplineND() ;

      TObject * clone(const char *newname) const ;

      TGraph* getGraph(const char *xvar, double step) ;

    protected:
        Double_t evaluate() const;

    private:
        RooListProxy vars_;
 
	mutable std::vector<double> w_;
	mutable std::map<int,std::vector<double> > v_map;
	mutable std::map<int,std::pair<double,double> > r_map;
	
	int ndim_;
	int M_;
	double eps_;
  	double axis_pts_;

	void calculateWeights(std::vector<double> &);
	double getDistSquare(int i, int j);
	double getDistFromSquare(int i) const;
	double radialFunc(double d2, double eps) const;
	

  ClassDef(RooSplineND,1) // MultiDim interpolations
};

#endif
