//

#include "blmr.h"




void Cblmr::set_x(const double *Xdata)
// precalculate working variables that depend only on x values
// and Sigma matrix
{
	int i,j;
	const double zero_eq = ldexp(1.,-40);

	int nstest=0;
	bool Xbad = false;
	for (i=0;i<n;i++) {

		if ( isinf( *(Xdata+i) ) || isnan( *(Xdata+i) ) ) { 
			Rcout << "Invalid x[i] value " << *(Xdata+i) << endl;
			Rcout << "x values not assigned";
			Rcout << endl;
			Xbad = true;
			break;
		}

		if ( i>0 && ( *(Xdata+i-1) > *(Xdata+i) ) ) {
			Rcout << "x[" << i-1 << "] =" << *(Xdata+i-1) 
				<< " and x[" << i << "] =" << *(Xdata+i) << endl;
			Rcout << "x values must be non-decreasing." << endl;
			Rcout << "x values not assigned";
			Rcout << endl;
			Xbad = true;
			break;
		}

		if ( *(Xdata+i) != *(Xdata+i-1) ) nstest++;
	}
	if (Xbad) return;
	nstest++;

	if ( (nstest<4 && model==M1) || (nstest<3 && model==M2) ) {
		Rcout << "Invalid number of seperate x- points: " << nstest << endl;
		Rcout << "Number must be at least 4 in M1 or 3 in M2." << endl;
		Rcout << "x values not assigned";
		Rcout << endl;
		return;
	}



	Vector<double> x(n);
	for (i=0;i<n;i++) {
		Xd[i] = Xdata[i];
		x[i] = Xdata[i];
	}
	*px = x;


	const double minXdiff = (Xd[n-1]-Xd[0])*0.001;

	ns= 0;
	for (i=1;i<n;i++) {
		if ( Xd[i] != Xd[i-1] ) {
			is[ns] = i-1;
			ns++;
			if ( Xd[i] - Xd[i-1] < minXdiff  ) {
				Rcout << "x[" << i-1 << "] =" << Xd[i-1] << " and x[" << i << "] =" << Xd[i] << endl;
				Rcout << "Warning: x values too close for reliable computations." << endl;
				Rcout << "Suggest combining these data into a single x value." << endl;
			}
		}
	}

	is[ns] = n-1;
	ns++;



	Vector<double> v1(n,1.), s1(n), sx(n);
	s1 = *pirS*v1;
	sx = *pirS*x;
	*psig1 = s1;
	*psigx = sx;

	s11 = s1*s1;
	sx1 = sx*s1;
	sxx = sx*sx;
	n1 = sqrt(s11);
	*pv1h= 1./n1 * s1;

	double x1 = sx*(*pv1h);
	*pxh = 1./sqrt( sx*sx - x1*x1) * ( sx - x1*(*pv1h) );



	get_Q();
	Matrix<double> QS(m,n);
	QS = (*pQ)*(*pirS);
	for (i=0;i<m;i++) for (j=0;j<n;j++) if (fabs(QS[i][j])<zero_eq) QS[i][j]=0.;

	Vector<double> e1(n,0.), en(n,0.),sen(n),dummy(m,0.),dummy2(n,0.);
	e1[0] = 1.;
	en[n-1] = 1.;

	*vund = dummy;
	(*vund)[0] = und;
	*vund2 = dummy2;
	(*vund2)[0] = und;
	*pnse1 = -1.*(*pirS*e1);
	se1sq = *pnse1*(*pnse1);
	*pnuse1 = 1./sqrt(se1sq) * (*pnse1);
	*puqe1 = 1./sqrt( (QS*e1)*(QS*e1) ) * (QS*e1);
	sen = *pirS*en;
	*pusen = 1./sqrt(sen*sen) * sen;
	*puqen = 1./sqrt( (QS*en)*(QS*en) ) * (QS*en);
	*puqx = 1./sqrt( (QS*x)*(QS*x) ) * (QS*x);

	
	if (Xs != NULL) set_xs();


	return;
}


