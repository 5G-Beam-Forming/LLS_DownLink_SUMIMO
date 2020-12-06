/*--------------------------------------------------------------------------------------*/
// misc.cpp
// - Miscellaneous
//                                                                           by LIT KAIST
/*--------------------------------------------------------------------------------------*/

#include "stdafx.h"
#include "misc.h"



double seed = 4.0 ;

double rnd()	/* generate a uniform r.v */
{
    seed = fmod(RR * seed, MM);
    return(seed * 4.656612875e-10);
}



double gasdis()	/* generate a normal r.v */
{
    static double t = 0.0;
    double x, v1, v2, r;
    
    if (t == 0) {
	do {
	    v1 = 2.0 * rnd() - 1.0;
	    v2 = 2.0 * rnd() - 1.0;
	    r = v1 * v1 + v2 * v2;
	} while (r >= 1.0);

	r = sqrt(-2.0 * log(r) / r);
        t = v2 * r;
		return(v1 * r);
    }
    else {
        x = t;
        t = 0.0;
		return(x);
    }
}

/*--------- 전달 시 삭제 예정 --------*/