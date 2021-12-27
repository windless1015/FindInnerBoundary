#include "MeshCube.h"

#pragma   warning(push)
#pragma   warning(disable: 4018)
#pragma   warning(disable: 4244)
using namespace std;
//! C1连续插值
/*!
  \return J方向（U方向）插值
*/
void C1_Terp1(unsigned n, double* blen, double * rb, unsigned nrs, unsigned* ndxrb, int* start, int* iend, unsigned* rs,
	double* ptsXX, double* ptsYY, double* ptsZZ, int* itype, unsigned dim10, unsigned dim20, unsigned dim30,
	unsigned protyp, unsigned* proj, unsigned* prof, unsigned* blend);

//! C2 超限插值法
/*!
  \return K方向(V方向)插值
*/
void C2_Terp2(unsigned n1, unsigned n2, double* blen, double * rb, unsigned nrs, unsigned* ndxrb, int* start, int* iend, unsigned* rs,
	double* ptsXX, double* ptsYY, double* ptsZZ, int* itype, unsigned dim10, unsigned dim20, unsigned dim30,
	unsigned protyp, unsigned* proj, unsigned* prof, unsigned* blend);

//! C3超限插值法
/*!
  \return L方向插值
*/
void C3_Terp3(double* blen, double * rb, unsigned nrs, unsigned* ndxrb, int* start, int* iend, unsigned* rs,
	double* ptsXX, double* ptsYY, double* ptsZZ, int* itype, unsigned dim10, unsigned dim20, unsigned dim30,
	unsigned protyp, unsigned* proj, unsigned* prof, unsigned* blend);


//! 超限插值生成网格(基于Coons面)
/*!
  \return
*/
void MeshCube::operator()(double* ptsXX, double* ptsYY, double* ptsZZ, int* start, int* iend, unsigned dim1, unsigned dim2, unsigned dim3, unsigned dimax)
{
	long long dmrb = 0;
	dmrb += (iend[1] - start[1] + 1)*(iend[0] - start[0] + 1);
	dmrb += 2 * (iend[1] - start[1] + iend[0] - start[0])*(iend[2] - start[2] + 1);
	dmrb += (iend[1] - start[1] + 1)*(iend[0] - start[0] + 1);
	unsigned JDim = dim1;
	unsigned KDim = dim2;
	unsigned LDim = dim3;
	unsigned JKMulDim = JDim * KDim;
	unsigned j = 0, k = 0, l = 0;
	unsigned dim123 = JKMulDim * LDim;
	double *facin = (double*)malloc(sizeof(double)*dim123 * 3);
	double *blen = (double*)malloc(sizeof(double)*dim123 * 12);
	unsigned *ndxrb = (unsigned*)malloc(sizeof(unsigned)*dim123);
	int *itype = (int*)malloc(sizeof(int)*dim123);
	double *arc = (double*)malloc(sizeof(double) * 12 * dimax);
	double *ars = (double*)malloc(sizeof(double) * 12 * dimax*dimax);
	double *dar = (double *)malloc(sizeof(double)*dimax * 4);
	double *tars = (double *)malloc(sizeof(double)*dimax * 12);
	double *arsi = (double *)malloc(sizeof(double)*dimax * 2);
	double *rb = (double *)malloc(sizeof(double) * 24 * dmrb);
	double *fi = (double *)malloc(sizeof(double)*dimax * 3);
	memset(facin, 0.0, sizeof(double)*dim123 * 3);
	memset(blen, 0.0, sizeof(double)*dim123 * 12);
	memset(ndxrb, 0, sizeof(unsigned)*dim123);
	memset(itype, 0, sizeof(int)*dim123);
	memset(arc, 0.0, sizeof(double) * 12 * dimax);
	memset(ars, 0.0, sizeof(double) * 12 * dimax*dimax);
	memset(dar, 0.0, sizeof(double)*dimax * 4);
	memset(tars, 0.0, sizeof(double)*dimax * 12);
	memset(arsi, 0.0, sizeof(double)*dimax * 2);
	memset(rb, 0.0, sizeof(double) * 24 * dmrb);
	memset(fi, 0.0, sizeof(double)*dimax * 3);
	unsigned ptsD3Sz = JKMulDim;
	unsigned ptsD2Sz = JDim;
	double zo = 1.0e-15;
	unsigned lfix = 555;
	unsigned cy[3][2] = { 1, 2, 2, 0, 0, 1 }; 
	unsigned d[3][3] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
	double tarc[2][2][3]; 
	unsigned nrs = 0; 
	start[0] = 1;  start[1] = 1;	start[2] = 1; 
	for (l = start[2] - 1; l < iend[2]; l++)
		for (k = start[1] - 1; k < iend[1]; k++)
			for (j = start[0] - 1; j < iend[0]; j++)
				itype[l*JKMulDim + k * JDim + j] = 0;
	if (start[0] != iend[0])
		for (l = start[2] - 1; l < iend[2]; l++)
			for (k = start[1] - 1; k < iend[1]; k++)
				for (j = start[0] - 1; j < iend[0]; j += iend[0] - start[0])
					itype[l*JKMulDim + k * JDim + j] = 555;
	if (start[1] != iend[1])
		for (l = start[2] - 1; l < iend[2]; l++)
			for (k = start[1] - 1; k < iend[1]; k += iend[1] - start[1])
				for (j = start[0] - 1; j < iend[0]; j++)
					itype[l*JKMulDim + k * JDim + j] = 555;
	if (start[2] != iend[2])
		for (l = start[2] - 1; l < LDim; l += iend[2] - start[2])
			for (k = start[1] - 1; k < KDim; k++)
				for (j = start[0] - 1; j < JDim; j++)
					itype[l*JKMulDim + k * JDim + j] = 555;
	unsigned prof[3] = { 0, 0, 0 };
	unsigned proj[3] = { 99, 1, 1 };
	unsigned blend[3] = { 222, 222, 222 };
	unsigned protyp = 912312;
	unsigned rs[3] = { 0, 1, 2 };
	for (unsigned ci = 0; ci < 3; ci++) 
	{
		int i1, i2;
		if (rs[ci] != 9999)	nrs = nrs + 1;
		if (iend[ci] >= start[ci]) {
			i1 = start[ci]; i2 = iend[ci];
		}
		else {
			i2 = start[ci]; i1 = iend[ci];
		}
		start[ci] = i1;
		iend[ci] = i2;
	}
	unsigned dimold = 0; 
	unsigned sam = 4;
	unsigned dir[3];
	if (proj[0] == 99)
	{
		for (unsigned ci = 0; ci < 3; ci++)
		{
			if (start[ci] != iend[ci])
			{
				++dimold;
				dir[dimold - 1] = ci;
			}
			else
			{
				--sam;
				dir[sam - 1] = ci;
			}
		}
		if (dimold == 1)
		{
			proj[0] = dir[0];
			proj[1] = 99;
			proj[2] = 99;
			blend[0] = 143190;
			blend[1] = 143190;
			blend[2] = 143190;
			protyp = 912312;
		}
		else if (dimold == 2)
		{
			if (protyp == 912312) {
				proj[0] = dir[0];
				proj[1] = dir[1];
				proj[2] = 99; 
			}
			else {
				proj[0] = dir[2];
				proj[1] = 99;
				proj[2] = 99;
				protyp = 79102;
			}
		}
		else
		{
			proj[0] = dir[0];
			proj[1] = dir[1];
			proj[2] = dir[2];
		}
	}
	else {
		for (unsigned ci = 0; ci < 3; ci++) {
			if (proj[ci] != 99)  ++dimold;
			else    --sam;
		}
	}
	if (proj[0] == 99)
	{
		return;
	}
	unsigned is1, is2, is3, ie1, ie2, ie3;
	if (iend[0] >= start[0]) { is1 = start[0]; ie1 = iend[0]; }
	else { ie1 = start[0]; is1 = iend[0]; }
	if (iend[1] >= start[1]) { is2 = start[1]; ie2 = iend[1]; }
	else { ie2 = start[1]; is2 = iend[1]; }
	if (iend[2] >= start[2]) { is3 = start[2]; ie3 = iend[2]; }
	else { ie3 = start[2]; is3 = iend[2]; }
	unsigned it1 = ie1 - is1;
	unsigned it2 = ie2 - is2;
	unsigned it3 = ie3 - is3;
	unsigned c[3];
	for (unsigned n = 0; n < 3; n++)
	{
		unsigned n1 = cy[n][0];
		unsigned n2 = cy[n][1];
		unsigned ndir = 0;
		if (start[n1] != iend[n1])  ++ndir;
		if (start[n2] != iend[n2])  ++ndir;
		unsigned isn = start[n] + 1;
		if ((blend[n] == 222) && (start[n] != iend[n]))
		{
			unsigned curBase = (start[n] - 1) * 12 + n;
			arc[curBase] = 0.0;
			arc[curBase + 6] = 0.0;
			arc[curBase + 3] = 0.0;
			arc[curBase + 9] = 0.0;
			for (unsigned ci = isn; ci <= iend[n]; ci++)
			{
				double darc = 0.0;
				c[n] = ci;
				c[n1] = start[n1];
				c[n2] = start[n2];
				unsigned base1 = (c[2] - 1)*ptsD3Sz + (c[1] - 1)*ptsD2Sz + c[0] - 1;
				unsigned base2 = (c[2] - d[n][2] - 1)*ptsD3Sz + (c[1] - d[n][1] - 1)*ptsD2Sz + (c[0] - d[n][0] - 1);
				double tmpD = ptsXX[base1] - ptsXX[base2];				darc += tmpD * tmpD;
				tmpD = ptsYY[base1] - ptsYY[base2];				darc += tmpD * tmpD;
				tmpD = ptsZZ[base1] - ptsZZ[base2];				darc += tmpD * tmpD;
				dar[(ci - 1) * 4] = darc;
				c[n1] = start[n1];
				c[n2] = iend[n2];
				darc = 0.0;
				base1 = (c[2] - 1)*ptsD3Sz + (c[1] - 1)*ptsD2Sz + c[0] - 1;
				base2 = (c[2] - d[n][2] - 1)*ptsD3Sz + (c[1] - d[n][1] - 1)*ptsD2Sz + (c[0] - d[n][0] - 1);
				tmpD = ptsXX[base1] - ptsXX[base2];				darc += tmpD * tmpD;
				tmpD = ptsYY[base1] - ptsYY[base2];				darc += tmpD * tmpD;
				tmpD = ptsZZ[base1] - ptsZZ[base2];				darc += tmpD * tmpD;
				dar[(ci - 1) * 4 + 2] = darc;
				c[n1] = iend[n1];
				c[n2] = start[n2];
				darc = 0.0;
				base1 = (c[2] - 1)*ptsD3Sz + (c[1] - 1)*ptsD2Sz + c[0] - 1;
				base2 = (c[2] - d[n][2] - 1)*ptsD3Sz + (c[1] - d[n][1] - 1)*ptsD2Sz + (c[0] - d[n][0] - 1);
				tmpD = ptsXX[base1] - ptsXX[base2];				darc += tmpD * tmpD;
				tmpD = ptsYY[base1] - ptsYY[base2];				darc += tmpD * tmpD;
				tmpD = ptsZZ[base1] - ptsZZ[base2];				darc += tmpD * tmpD;
				dar[(ci - 1) * 4 + 1] = darc;
				c[n1] = iend[n1];
				c[n2] = iend[n2];
				darc = 0.0;
				base1 = (c[2] - 1)*ptsD3Sz + (c[1] - 1)*ptsD2Sz + c[0] - 1;
				base2 = (c[2] - d[n][2] - 1)*ptsD3Sz + (c[1] - d[n][1] - 1)*ptsD2Sz + (c[0] - d[n][0] - 1);
				tmpD = ptsXX[base1] - ptsXX[base2];				darc += tmpD * tmpD;
				tmpD = ptsYY[base1] - ptsYY[base2];				darc += tmpD * tmpD;
				tmpD = ptsZZ[base1] - ptsZZ[base2];				darc += tmpD * tmpD;
				dar[(ci - 1) * 4 + 3] = darc;
			}
			for (unsigned ci = isn - 1; ci < iend[n]; ci++)
			{
				unsigned curBase = ci * 4;
				dar[curBase] = sqrt(dar[curBase]);
				dar[curBase + 2] = sqrt(dar[curBase + 2]);
				dar[curBase + 1] = sqrt(dar[curBase + 1]);
				dar[curBase + 3] = sqrt(dar[curBase + 3]);
			}
			for (unsigned ci = isn - 1; ci < iend[n]; ci++)
			{
				unsigned curBase = ci * 12 + n;
				unsigned curBaseDr = ci * 4;
				arc[curBase] = arc[curBase - 12] + dar[curBaseDr];
				arc[curBase + 6] = arc[curBase - 6] + dar[curBaseDr + 2];
				arc[curBase + 3] = arc[curBase - 9] + dar[curBaseDr + 1];
				arc[curBase + 9] = arc[curBase - 3] + dar[curBaseDr + 3];
			}

			if (ndir == 2)
			{
				is1 = start[n1] + 1;
				ie1 = iend[n1] - 1;
				double darc = 0.0;
				for (unsigned i1 = is1; i1 <= ie1; i1++)
				{
					c[n1] = i1;
					unsigned curBase = (start[n] - 1)*dimax * 12 + (i1 - 1) * 12 + 3 + n;
					ars[curBase] = 0.0;
					ars[curBase + 6] = 0.0;
					for (unsigned ci = isn; ci <= iend[n]; ci++)
					{
						c[n] = ci; 
						c[n2] = start[n2];
						darc = 0.0;
						unsigned base1 = (c[2] - 1)*ptsD3Sz + (c[1] - 1)*ptsD2Sz + c[0] - 1;
						unsigned base2 = (c[2] - d[n][2] - 1)*ptsD3Sz + (c[1] - d[n][1] - 1)*ptsD2Sz + (c[0] - d[n][0] - 1);
						double tmpD = ptsXX[base1] - ptsXX[base2];				darc += tmpD * tmpD;
						tmpD = ptsYY[base1] - ptsYY[base2];				darc += tmpD * tmpD;
						tmpD = ptsZZ[base1] - ptsZZ[base2];				darc += tmpD * tmpD;
						dar[(ci - 1) * 4] = darc; 
						c[n2] = iend[n2]; 
						darc = 0.0;
						base1 = (c[2] - 1)*ptsD3Sz + (c[1] - 1)*ptsD2Sz + c[0] - 1;
						base2 = (c[2] - d[n][2] - 1)*ptsD3Sz + (c[1] - d[n][1] - 1)*ptsD2Sz + (c[0] - d[n][0] - 1);
						tmpD = ptsXX[base1] - ptsXX[base2];				darc += tmpD * tmpD;
						tmpD = ptsYY[base1] - ptsYY[base2];				darc += tmpD * tmpD;
						tmpD = ptsZZ[base1] - ptsZZ[base2];				darc += tmpD * tmpD;
						dar[(ci - 1) * 4 + 1] = darc;
					}
					for (unsigned ci = isn - 1; ci < iend[n]; ci++) //do 153
					{
						unsigned curBase = ci * 4;
						dar[curBase] = sqrt(dar[curBase]);
						dar[curBase + 1] = sqrt(dar[curBase + 1]);
					}
					for (unsigned ci = isn - 1; ci < iend[n]; ci++) //do 154
					{
						unsigned curBase = ci * dimax * 12 + (i1 - 1) * 12 + 3 + n;   //ars( 3 ,   2, 0:1 , dimax , dimax )
						ars[curBase] = ars[curBase - 12 * dimax] + dar[ci * 4];
						ars[curBase + 6] = ars[curBase - 12 * dimax + 6] + dar[ci * 4 + 1];
					}
				}
				is2 = start[n2] + 1;
				ie2 = iend[n2] - 1;
				for (unsigned i2 = is2; i2 <= ie2; i2++) //do 158
				{
					c[n2] = i2;
					unsigned tBase = (start[n] - 1)*dimax * 12 + (i2 - 1) * 12 + n;
					ars[tBase] = 0.0;
					ars[tBase + 6] = 0.0;
					for (unsigned ci = isn; ci <= iend[n]; ci++) //do 136
					{
						c[n] = ci;
						c[n1] = start[n1];
						darc = 0.0;
						unsigned base1 = (c[2] - 1)*ptsD3Sz + (c[1] - 1)*ptsD2Sz + c[0] - 1;
						unsigned base2 = (c[2] - d[n][2] - 1)*ptsD3Sz + (c[1] - d[n][1] - 1)*ptsD2Sz + (c[0] - d[n][0] - 1);
						double tmpD = ptsXX[base1] - ptsXX[base2];				darc += tmpD * tmpD;
						tmpD = ptsYY[base1] - ptsYY[base2];				darc += tmpD * tmpD;
						tmpD = ptsZZ[base1] - ptsZZ[base2];				darc += tmpD * tmpD;
						dar[(ci - 1) * 4] = darc;   //dar( 0:1 , 0:1 , dimax )
						c[n1] = iend[n1];
						darc = 0.0;
						base1 = (c[2] - 1)*ptsD3Sz + (c[1] - 1)*ptsD2Sz + c[0] - 1;
						base2 = (c[2] - d[n][2] - 1)*ptsD3Sz + (c[1] - d[n][1] - 1)*ptsD2Sz + (c[0] - d[n][0] - 1);
						tmpD = ptsXX[base1] - ptsXX[base2];				darc += tmpD * tmpD;
						tmpD = ptsYY[base1] - ptsYY[base2];				darc += tmpD * tmpD;
						tmpD = ptsZZ[base1] - ptsZZ[base2];				darc += tmpD * tmpD;
						dar[(ci - 1) * 4 + 1] = darc;   //dar( 0:1 , 0:1 , dimax )
					}
					for (unsigned ci = isn - 1; ci < iend[n]; ci++) //do 156
					{
						dar[ci * 4] = sqrt(dar[ci * 4]);
						dar[ci * 4 + 1] = sqrt(dar[ci * 4 + 1]);
					}
					for (unsigned ci = isn - 1; ci < iend[n]; ci++) //do 157
					{
						unsigned curBase = ci * dimax * 12 + (i2 - 1) * 12 + n;   //ars( 3 ,   2, 0:1 , dimax , dimax )
						ars[curBase] = ars[curBase - 12 * dimax] + dar[ci * 4];
						ars[curBase + 6] = ars[curBase - 12 * dimax + 6] + dar[ci * 4 + 1];
					}
				}
			}
			curBase = (isn - 1) * 12 + n;			 
			if ((arc[curBase] + arc[curBase + 6] + arc[curBase + 3] + arc[curBase + 9]) < zo)
			{
				blend[n] = 143190;
				continue;
			}
			if (arc[(iend[n] - 1) * 12 + n] < zo)
			{
				double fac = 0.0;
				if (ndir == 2)
				{
					fac = 1.0;
					unsigned curBase = (iend[n] - 1)*dimax * 12 + n;   //ars( 3 ,   2, 0:1 , dimax , dimax )
					if ((ars[curBase + (is1 - 1) * 12 + 3] > zo) && (ars[curBase + (is2 - 1) * 12] > zo))
						fac = 0.5;
					for (unsigned ci = isn - 1; ci < iend[n]; ci++) //do 41
					{
						unsigned curBase = ci * 12 * dimax + n;  ////ars( 3 ,   2, 0:1 , dimax , dimax )
						arc[ci * 12 + n] = (ars[curBase + (is1 - 1) * 12 + 3] + ars[curBase + (is2 - 1) * 12])*fac;
					}
				}
				else
				{
					for (unsigned ci = isn - 1; ci < iend[n]; ci++) //do 141
					{
						unsigned curBase = ci * 12 + n;
						arc[curBase] = max(arc[curBase + 6], arc[curBase + 3]);
					}
				}
			}
			if (arc[(iend[n] - 1) * 12 + n + 3] < zo)
			{
				double fac = 0.0;
				if (ndir == 2)
				{
					fac = 1.0;
					unsigned curBase = (iend[n] - 1)*dimax * 12 + n; 
					if ((ars[curBase + (ie1 - 1) * 12 + 3] > zo) && (ars[curBase + (is2 - 1) * 12 + 6] > zo))
						fac = 0.5;
					for (unsigned ci = isn - 1; ci < iend[n]; ci++)
					{
						unsigned curBase = ci * 12 * dimax + n;
						arc[ci * 12 + n + 3] = (ars[curBase + (ie1 - 1) * 12 + 3] + ars[curBase + (is2 - 1) * 12 + 6])*fac;
					}
				}
				else
				{
					for (unsigned ci = isn - 1; ci < iend[n]; ci++) 
					{
						unsigned curBase = ci * 12 + n;
						arc[curBase + 3] = max(arc[curBase], arc[curBase + 9]);
					}
				}
			}
			if (arc[(iend[n] - 1) * 12 + n + 6] < zo)
			{
				double fac = 0.0;
				if (ndir == 2)
				{
					fac = 1.0;
					unsigned curBase = (iend[n] - 1)*dimax * 12 + n;   //ars( 3 ,   2, 0:1 , dimax , dimax )
					if ((ars[curBase + (is1 - 1) * 12 + 9] > zo) && (ars[curBase + (ie2 - 1) * 12] > zo))
						fac = 0.5;
					for (unsigned ci = isn - 1; ci < iend[n]; ci++) //do 43   //arc( 3 , 0:1, 0:1 , dimax )	
					{
						unsigned curBase = ci * 12 * dimax + n;  ////ars( 3 ,   2, 0:1 , dimax , dimax )
						arc[ci * 12 + n + 6] = (ars[curBase + (is1 - 1) * 12 + 9] + ars[curBase + (ie2 - 1) * 12])*fac;
					}
				}
				else
				{
					for (unsigned ci = isn - 1; ci < iend[n]; ci++) //do 143
					{
						unsigned curBase = ci * 12 + n;
						arc[curBase + 6] = max(arc[curBase], arc[curBase + 9]);
					}
				}
			}
			if (arc[(iend[n] - 1) * 12 + n + 9] < zo)
			{
				double fac = 0.0;
				if (ndir == 2)
				{
					fac = 1.0;
					unsigned curBase = (iend[n] - 1)*dimax * 12 + n;   //ars( 3 ,   2, 0:1 , dimax , dimax )
					if ((ars[curBase + (ie1 - 1) * 12 + 9] > zo) && (ars[curBase + (ie2 - 1) * 12 + 6] > zo))
						fac = 0.5;
					for (unsigned ci = isn - 1; ci < iend[n]; ci++) //do 44   //arc( 3 , 0:1, 0:1 , dimax )	
					{
						unsigned curBase = ci * 12 * dimax + n + 6;  ////ars( 3 ,   2, 0:1 , dimax , dimax )
						arc[ci * 12 + n + 9] = (ars[curBase + (ie1 - 1) * 12 + 9] + ars[curBase + (ie2 - 1) * 12 + 6])*fac;
					}
				}
				else
				{
					for (unsigned ci = isn - 1; ci < iend[n]; ci++) //do 144
					{
						unsigned curBase = ci * 12 + n;
						arc[curBase + 9] = max(arc[curBase + 3], arc[curBase + 6]);
					}
				}
			}

			if (ndir == 2)
			{
				if (start[n1] == iend[n1]) cout << "ft1--n1: " << n1 << " " << start[n1] << std::endl;
				if (start[n2] == iend[n2]) cout << "ft2--n2: " << n2 << " " << start[n2] << std::endl;
				double ft1 = 1.0 / float(start[n1] - iend[n1]);
				double ft2 = 1.0 / float(start[n2] - iend[n2]);
				double fac = 0.0, fac1 = 0.0;
				for (unsigned m = 0; m <= 1; m++)
				{
					unsigned curBase1 = (iend[n] - 1)*dimax * 12 + 6 * m + 3 + n;
					for (unsigned i = is1; i <= ie1; i++)
					{
						if (ars[curBase1 + (i - 1) * 12] < zo)
						{
							fac = (i - start[n])*ft1;
							fac1 = 1.0 - fac;
							unsigned curBase2 = (i - 1) * 12 + 6 * m + 3 + n;
							for (unsigned ci = isn - 1; ci < iend[n]; ci++) //do 45
								ars[ci*dimax * 12 + curBase2] = fac1 * arc[12 * ci + m * 6] + fac * arc[12 * ci + m * 6 + 3];
						}
					}
				}
				for (unsigned m = 0; m <= 1; m++) {
					unsigned curBase1 = (iend[n] - 1)*dimax * 12 + 6 * m + n;
					for (unsigned i = is2; i <= ie2; i++)
					{
						if (ars[curBase1 + (i - 1) * 12] < zo)
						{
							fac = (i - start[n])*ft2;
							fac1 = 1.0 - fac;
							unsigned curBase2 = (i - 1) * 12 + 6 * m + n;
							for (unsigned ci = isn - 1; ci < iend[n]; ci++)
								ars[ci*dimax * 12 + curBase2] = fac1 * arc[12 * ci + m * 3 + n] + fac * arc[12 * ci + m * 3 + 6 + n];
						}
					} 
				} 
			} 

			for (unsigned m1 = 0; m1 <= 1; m1++) {
				for (unsigned m2 = 0; m2 <= 1; m2++) 
				{
					tarc[m2][m1][n] = arc[(iend[n] - 1) * 12 + m2 * 6 + m1 * 3 + n];
					double arci = 0.0;
					if (tarc[m2][m1][n] == 0.0)
						cout << "zero=arci:Error: n,m1,m2 " << n << " " << m1 << " " << m2 << std::endl;
					else
						arci = 1.0 / tarc[m2][m1][n]; //LFS 2016.10.17
					for (unsigned ci = isn - 1; ci < iend[n]; ci++) //do 40
						arc[12 * ci + m2 * 6 + m1 * 3 + n] *= arci;
				}
			}


			if (ndir == 2)
			{
				for (unsigned i1 = is1; i1 <= ie1; i1++) // do 145 
				{
					unsigned curBase = (iend[n] - 1)*dimax * 12 + (i1 - 1) * 12 + 3 + n;
					tars[(i1 - 1) * 12 + 3 + n] = ars[curBase];           //tars( 3 , 2 , 0:1 , dimax )
					tars[(i1 - 1) * 12 + 9 + n] = ars[curBase + 6];
					if (ars[curBase] == 0.0)
						cout << "By C version---arsi:1i1,n,end " << i1 << " " << n << " " << iend[n] << std::endl;
					if (ars[curBase + 6] == 0.0)
						cout << "By C version---arsi:2i1,n,end " << i1 << " " << n << " " << iend[n] << std::endl;
					arsi[(i1 - 1) * 2] = 1.0 / ars[curBase];          //arsi( 0:1 , dimax )
					arsi[(i1 - 1) * 2 + 1] = 1.0 / ars[curBase + 6];
				}// ends do 145
				for (unsigned i1 = is1 - 1; i1 < ie1; i1++) {             //do 245
					for (unsigned ci = isn - 1; ci < iend[n]; ci++) { //do 245
						unsigned curBase = ci * dimax * 12 + i1 * 12 + 3 + n;
						ars[curBase] *= arsi[i1 * 2];
						ars[curBase + 6] *= arsi[i1 * 2 + 1];
					}
				}// ends do 245

				for (unsigned i2 = is2 - 1; i2 < ie2; i2++) // do 146 
				{
					unsigned curBase = (iend[n] - 1)*dimax * 12 + i2 * 12 + n;
					tars[i2 * 12 + n] = ars[curBase];           //tars( 3 , 2 , 0:1 , dimax )
					tars[i2 * 12 + 6 + n] = ars[curBase + 6];          //ars( 3 ,   2, 0:1 , dimax , dimax )
					if (ars[curBase] == 0.0)
						cout << "By C version---arsi:1i1,n,end " << i2 << " " << n << " " << iend[n] << std::endl;
					if (ars[curBase + 6] == 0.0)
						cout << "By C version---arsi:2i1,n,end " << i2 << " " << n << " " << iend[n] << std::endl;
					arsi[i2 * 2] = 1.0 / ars[curBase];          //arsi( 0:1 , dimax )
					arsi[i2 * 2 + 1] = 1.0 / ars[curBase + 6];
				}// ends do 146
				for (unsigned i2 = is2 - 1; i2 < ie2; i2++) {             //do 246
					for (unsigned ci = isn - 1; ci < iend[n]; ci++) { //do 246
						unsigned curBase = ci * dimax * 12 + i2 * 12 + n;
						ars[curBase] *= arsi[i2 * 2];
						ars[curBase + 6] *= arsi[i2 * 2 + 1];
					}
				}// ends do 246

				for (unsigned ci = start[n] - 1; ci < iend[n]; ci++) //do 147
				{
					unsigned curBase = ci * dimax * 12 + n;
					unsigned curBase2 = ci * 12 + n;
					ars[curBase + (start[n2] - 1) * 12] = arc[curBase2];
					ars[curBase + (iend[n2] - 1) * 12] = arc[curBase2 + 6];
					ars[curBase + (start[n2] - 1) * 12 + 6] = arc[curBase2 + 3];
					ars[curBase + (iend[n2] - 1) * 12 + 6] = arc[curBase2 + 9];
					ars[curBase + (start[n1] - 1) * 12 + 3] = arc[curBase2];
					ars[curBase + (iend[n1] - 1) * 12 + 3] = arc[curBase2 + 3];
					ars[curBase + (start[n1] - 1) * 12 + 9] = arc[curBase2 + 6];
					ars[curBase + (iend[n1] - 1) * 12 + 9] = arc[curBase2 + 9];
				} //ends do 147
				tars[(start[n2] - 1) * 12 + n] = tarc[0][0][n];
				tars[(iend[n2] - 1) * 12 + n] = tarc[1][0][n];
				tars[(start[n2] - 1) * 12 + n + 6] = tarc[0][1][n];
				tars[(iend[n2] - 1) * 12 + n + 6] = tarc[1][1][n];
				tars[(start[n1] - 1) * 12 + n + 3] = tarc[0][0][n];
				tars[(iend[n1] - 1) * 12 + n + 3] = tarc[0][1][n];
				tars[(start[n1] - 1) * 12 + n + 9] = tarc[1][0][n];
				tars[(iend[n1] - 1) * 12 + n + 9] = tarc[1][1][n];
			} 
		}	
	}
	unsigned ndx = 0;
	for (unsigned i2 = start[1] - 1; i2 < iend[1]; i2++) {   //do 70
		unsigned curBase = (start[2] - 1)*JKMulDim + i2 * dim1;
		for (unsigned i1 = start[0] - 1; i1 < iend[0]; i1++)
		{
			ndxrb[curBase + i1] = ++ndx;
		}
	}// ends do 70
	if (start[2] != iend[2])
	{
		for (unsigned i3 = start[2] + 1; i3 <= iend[2] - 1; i3++) //do 75
		{
			unsigned curBase = (i3 - 1)*JKMulDim;
			unsigned curBase1 = curBase + (start[1] - 1)*dim1;
			for (unsigned i1 = start[0] - 1; i1 < iend[0]; i1++)
			{
				ndxrb[curBase1 + i1] = ++ndx;
			}//end do 71
			if (start[1] != iend[1])
			{
				for (unsigned i2 = start[1] + 1; i2 <= iend[1] - 1; i2++)
				{
					ndxrb[curBase + (i2 - 1)*dim1 + start[0] - 1] = ++ndx;
					ndxrb[curBase + (i2 - 1)*dim1 + iend[0] - 1] = ++ndx;
				}//end do 72
				for (unsigned i1 = start[0] - 1; i1 < iend[0]; i1++)
				{
					ndxrb[curBase + (iend[1] - 1)*dim1 + i1] = ++ndx;
				}//end do 73
			}//ends if
		}//ends do 75
		for (unsigned i2 = start[1] - 1; i2 < iend[1]; i2++) //do 76
		{
			unsigned curBase = (iend[2] - 1)*JKMulDim + i2 * dim1;
			for (unsigned i1 = start[0] - 1; i1 < iend[0]; i1++)
			{
				ndxrb[curBase + i1] = ++ndx;
			}//end do 73
		}//ends do 76
	}// ends if start[2] != iend[2]

	for (unsigned nn = 0; nn < dimold; nn++) //do 110
	{
		unsigned n = proj[nn];
		unsigned n1 = cy[n][0];
		unsigned n2 = cy[n][1];
		for (unsigned i = start[n]; i <= iend[n]; i += iend[n] - start[n])//do 110
		{
			c[n] = i;
			//for(unsigned ri=0;ri<nrs;ri++)//do 110
			{
				if (iend[n1] - start[n1] > iend[n2] - start[n2])
				{
					for (unsigned i2 = start[n2]; i2 <= iend[n2]; i2++) //do 1104
					{
						c[n2] = i2;
						if (n1 == 0)
						{
							unsigned curBase1 = (c[2] - 1)*JKMulDim + (c[1] - 1)*dim1;
							unsigned curBase2 = (c[2] - 1)*ptsD3Sz + (c[1] - 1)*ptsD2Sz;
							for (unsigned i1 = start[n1]; i1 <= iend[n1]; i1++) //do 1101
							{
								ndx = ndxrb[curBase1 + i1 - 1];  //points( 3 , 0:dim1 , 0:dim2 , 0:dim3 )
								rb[(ndx - 1) * 24] = ptsXX[curBase2 + i1 - 1];
								rb[(ndx - 1) * 24 + 1] = ptsYY[curBase2 + i1 - 1];
								rb[(ndx - 1) * 24 + 2] = ptsZZ[curBase2 + i1 - 1];

							}//do 1101
						}
						else if (n1 == 1)
						{
							unsigned curBase1 = (c[2] - 1)*JKMulDim + c[0] - 1;
							unsigned curBase2 = (c[2] - 1)*ptsD3Sz + c[0] - 1;
							for (unsigned i1 = start[n1]; i1 <= iend[n1]; i1++)//do 1102
							{
								ndx = ndxrb[curBase1 + (i1 - 1)*dim1];  //points( 3 , 0:dim1 , 0:dim2 , 0:dim3 )
								rb[(ndx - 1) * 24] = ptsXX[curBase2 + (i1 - 1)*ptsD2Sz];
								rb[(ndx - 1) * 24 + 1] = ptsYY[curBase2 + (i1 - 1)*ptsD2Sz];
								rb[(ndx - 1) * 24 + 2] = ptsZZ[curBase2 + (i1 - 1)*ptsD2Sz];
							}
						}
						else if (n1 == 2)
						{
							unsigned curBase1 = (c[1] - 1)*dim1 + c[0] - 1;
							unsigned curBase2 = (c[1] - 1)*ptsD2Sz + c[0] - 1;
							for (unsigned i1 = start[n1]; i1 <= iend[n1]; i1++)//do 1103
							{
								ndx = ndxrb[(i1 - 1)*JKMulDim + curBase1]; 
								rb[(ndx - 1) * 24] = ptsXX[curBase2 + (i1 - 1)*ptsD3Sz];
								rb[(ndx - 1) * 24 + 1] = ptsYY[curBase2 + (i1 - 1)*ptsD3Sz];
								rb[(ndx - 1) * 24 + 2] = ptsZZ[curBase2 + (i1 - 1)*ptsD3Sz];
							}
						}
					}
				}
				else
				{
					for (unsigned i1 = start[n1]; i1 <= iend[n1]; i1++)//do 1114
					{
						c[n1] = i1;
						if (n2 == 0)
						{
							unsigned curBase1 = (c[2] - 1)*JKMulDim + (c[1] - 1)*dim1;
							unsigned curBase2 = (c[2] - 1)*ptsD3Sz + (c[1] - 1)*ptsD2Sz;
							for (unsigned i2 = start[n2]; i2 <= iend[n2]; i2++)//do 1111
							{
								ndx = ndxrb[curBase1 + i2 - 1];
								rb[(ndx - 1) * 24] = ptsXX[curBase2 + (i2 - 1)];
								rb[(ndx - 1) * 24 + 1] = ptsYY[curBase2 + (i2 - 1)];
								rb[(ndx - 1) * 24 + 2] = ptsZZ[curBase2 + (i2 - 1)];
							}//do 1111
						}
						else if (n2 == 1)
						{
							unsigned curBase1 = (c[2] - 1)*JKMulDim + c[0] - 1;
							unsigned curBase2 = (c[2] - 1)*ptsD3Sz + c[0] - 1;
							for (unsigned i2 = start[n2]; i2 <= iend[n2]; i2++)//do 1112
							{
								ndx = ndxrb[curBase1 + (i2 - 1)*dim1];  //points( 3 , 0:dim1 , 0:dim2 , 0:dim3 )
								rb[(ndx - 1) * 24] = ptsXX[curBase2 + (i2 - 1)*ptsD2Sz];
								rb[(ndx - 1) * 24 + 1] = ptsYY[curBase2 + (i2 - 1)*ptsD2Sz];
								rb[(ndx - 1) * 24 + 2] = ptsZZ[curBase2 + (i2 - 1)*ptsD2Sz];
							}//do 1112
						}
						else if (n2 == 2)
						{
							unsigned curBase1 = (c[1] - 1)*dim1 + c[0] - 1;
							unsigned curBase2 = (c[1] - 1)*ptsD2Sz + c[0] - 1;
							for (unsigned i2 = start[n2]; i2 <= iend[n2]; i2++)//do 1113
							{
								ndx = ndxrb[(i2 - 1)*JKMulDim + curBase1];  //points( 3 , 0:dim1 , 0:dim2 , 0:dim3 )
								rb[(ndx - 1) * 24] = ptsXX[curBase2 + (i2 - 1)*ptsD3Sz];
								rb[(ndx - 1) * 24 + 1] = ptsYY[curBase2 + (i2 - 1)*ptsD3Sz];
								rb[(ndx - 1) * 24 + 2] = ptsZZ[curBase2 + (i2 - 1)*ptsD3Sz];
							}//do 1113
						}// ends if n2==
					}//do 1114
				}//ends if iend---start
			}//ends do 110 ri
		}//ends do 110 i
	}//ends do 110 nn
	unsigned n, n1, n2, ndir;
	unsigned itn;
	double ft1, ft2, ftn;
	unsigned id[3];
	int m, l1, m1, l2, m2;
	double dr[4][4][3]; //dr( 3 , 0:3 , 0:3 )
	for (unsigned nn = 0; nn < dimold; nn++) //do 240
	{
		n = proj[nn];
		if (prof[n] != 0)
		{
			n1 = cy[n][0];	n2 = cy[n][1];	ndir = 0;
			if (start[n1] != iend[n1])	++ndir;
			if (start[n2] != iend[n2])	++ndir;
			it1 = iend[n1] - start[n1];
			it2 = iend[n2] - start[n2];
			itn = iend[n] - start[n];
			if (start[n1] != iend[n1]) { is1 = start[n1] + 1;   ie1 = iend[n1] - 1; }
			else { is1 = start[n1];    ie1 = iend[n1]; }
			if (start[n2] != iend[n2]) { is2 = start[n2] + 1;   ie2 = iend[n2] - 1; }
			else { is2 = start[n2];    ie2 = iend[n2]; }
			if (itn == 0) cout << "By C version---ftn:itn " << itn << std::endl;
			if (it1 != 0)  ft1 = 1.0 / float(it1);
			if (it2 != 0)  ft2 = 1.0 / float(it2);
			ftn = 1.0 / float(itn);
			id[n] = 1;	id[n1] = 0;		id[n2] = 0;
			l = 3;		m = -1;
			double dn = 0.0;
			double fal1 = 0.0, fal2 = 0.0;
			double fac1 = 0.0, fac2 = 0.0;
			double ofac1 = 0.0, ofac2 = 0.0;
			for (unsigned i = start[n]; i <= iend[n]; i += iend[n] - start[n]) //do 230
			{
				l = l - 2;	m = m + 1;	c[n] = i;
				for (unsigned i1 = is1; i1 <= ie1; i1++)   //do 230
				{
					c[n1] = i1;
					for (unsigned i2 = is2; i2 <= ie2; i2++)	//do 230
					{
						c[n2] = i2;
						ndx = ndxrb[(c[2] - 1)*JKMulDim + (c[1] - 1)*dim1 + c[0] - 1];
						unsigned curBase1 = (c[2] - 1)*ptsD3Sz + (c[1] - 1)*ptsD2Sz + c[0] - 1;
						unsigned curBase2 = d[n1][2] * ptsD3Sz + d[n1][1] * ptsD2Sz + d[n1][0];
						unsigned curBase3 = d[n2][2] * ptsD3Sz + d[n2][1] * ptsD2Sz + d[n2][0];
						dr[0][n1][0] = ptsXX[curBase1 + curBase2] - ptsXX[curBase1 - curBase2];   // do 210
						dr[0][n2][0] = ptsXX[curBase1 + curBase3] - ptsXX[curBase1 - curBase3];
						dr[0][n1][1] = ptsYY[curBase1 + curBase2] - ptsYY[curBase1 - curBase2];
						dr[0][n2][1] = ptsYY[curBase1 + curBase3] - ptsYY[curBase1 - curBase3];
						dr[0][n1][2] = ptsZZ[curBase1 + curBase2] - ptsZZ[curBase1 - curBase2];
						dr[0][n2][2] = ptsZZ[curBase1 + curBase3] - ptsZZ[curBase1 - curBase3];

						dn = 0.0;
						for (unsigned ri = 0; ri < 3; ri++) // do 220
						{
							dr[0][n][ri] = dr[0][n1][cy[ri][0]] * dr[0][n2][cy[ri][1]]
								- dr[0][n1][cy[ri][1]] * dr[0][n2][cy[ri][0]];
							dn += dr[0][n][ri] * dr[0][n][ri];
						}//ends do 220
						if (blend[n] == 222)  // blend=222 means arc blending interpaltion,else line interpaltion
						{
							if (start[n1] != iend[n1])		fal1 = (c[n1] - start[n1])*ft1;
							else	fal1 = 0.0;
							if (start[n2] != iend[n2])		fal2 = (c[n2] - start[n2])*ft2;
							else	fal2 = 0.0;
							if (blend[n1] == 222) 	//arc( 3 , 0:1, 0:1 , dimax )
								fac1 = (1.0 - fal2)*arc[(c[n1] - 1) * 12 + m * 6 + n1] + fal2 * arc[(c[n1] - 1) * 12 + m * 6 + 3 + n1];
							else fac1 = fal1;
							if (blend[n2] == 222) 	//arc( 3 , 0:1, 0:1 , dimax )
								fac2 = (1.0 - fal1)*arc[(c[n2] - 1) * 12 + m * 3 + n2] + fal1 * arc[(c[n2] - 1) * 12 + 6 + 3 * m + n2];
							else fac2 = fal2;
							ofac1 = 1.0 - fac1;		ofac2 = 1.0 - fac2;
							//此处省略span( c(1),c(2),c(3) ) .ne. 0.0条件分支
							double tal = 0.0;
							if (ndir == 1)
								tal = ofac1 * tars[(i2 - 1) * 12 + n] + fac1 * tars[(i2 - 1) * 12 + n + 6]
								+ ofac2 * tars[(i1 - 1) * 12 + n + 3] + fac2 * tars[(i1 - 1) * 12 + n + 9]
								- ofac1 * ofac2*tarc[0][0][n] - ofac1 * fac2*tarc[1][0][n]
								- fac1 * ofac2*tarc[0][1][n] - fac1 * fac2*tarc[1][1][n];
							else
								tal = ofac1 * ofac2*tarc[0][0][n] + ofac1 * fac2*tarc[1][0][n]
								+ fac1 * ofac2*tarc[0][1][n] + fac1 * fac2*tarc[1][1][n];
							if (dn == 0.0) cout << "dn:dn " << dn << std::endl;
							dn = tal / sqrt(dn);
						} //blend[n] == 222
						else
						{
							dn = 0.0;
							cout << "Maybe error! DN is: " << dn << std::endl;
							; //Since span is not initialized, not implemented. Please see to JavaScipt code
						} //blend[n] != 222
						for (unsigned ri = 0; ri < nrs; ri++) // do 230
							rb[(ndx - 1) * 24 + id[2] * 12 + id[1] * 6 + id[0] * 3 + rs[ri]] = dr[0][n][rs[ri]] * dn;  //rb( 3 , 0:1 , 0:1 , 0:1 , dmrb )
					}
				}
			}// ends do 230
		}// ends if if(prof[n] !=0)
	}// ends do 240
	double dn1, dn2, fac;
	if (dimold > 1)
	{
		for (unsigned n = sam - 1; n < 3; n++) //do 350
		{
			n1 = proj[cy[n][0]];		n2 = proj[cy[n][1]];  //proj数组的值不用减1；赋值时就已经做过处理
			if (prof[n1] * prof[n2] != 0) //JavaScipt的0
			{
				l1 = 3;		m1 = -1;
				for (unsigned i1 = start[n1]; i1 <= iend[n1]; i1 += iend[n1] - start[n1]) // do 340
				{
					c[n1] = i1;		l1 -= 2;	m1 += 1;
					l2 = 3;		m2 = -1;
					for (unsigned i2 = start[n2]; i2 <= iend[n2]; i2 += iend[n2] - start[n2]) // do 340
					{
						c[n2] = i2;	l2 -= 2; m2 += 1;
						for (unsigned i = start[n]; i <= iend[n]; i++) // do 340
						{
							c[n] = i; dn1 = 0.0; dn2 = 0.0;
							ndx = ndxrb[(c[2] - 1)*JKMulDim + (c[1] - 1)*dim1 + c[0] - 1];
							unsigned curBase1 = (c[2] - 1)*ptsD3Sz + (c[1] - 1)*ptsD2Sz + c[0] - 1;   // do 310
							unsigned curBase2 = l1 * d[n1][2] * ptsD3Sz + l1 * d[n1][1] * ptsD2Sz + l1 * d[n1][0];
							unsigned curBase3 = l2 * d[n2][2] * ptsD3Sz + l2 * d[n2][1] * ptsD2Sz + l2 * d[n2][0];
							dr[0][n1][0] = (ptsXX[curBase1 + curBase2] - ptsXX[curBase1]) *l1;
							dn1 += dr[0][n1][0] * dr[0][n1][0];
							dr[0][n1][1] = (ptsYY[curBase1 + curBase2] - ptsYY[curBase1]) *l1;
							dn1 += dr[0][n1][1] * dr[0][n1][1];
							dr[0][n1][2] = (ptsZZ[curBase1 + curBase2] - ptsZZ[curBase1]) *l1;
							dn1 += dr[0][n1][2] * dr[0][n1][2];

							dr[0][n2][0] = (ptsXX[curBase1 + curBase3] - ptsXX[curBase1])*l2;
							dn2 += dr[0][n2][0] * dr[0][n2][0];
							dr[0][n2][1] = (ptsYY[curBase1 + curBase3] - ptsYY[curBase1])*l2;
							dn2 += dr[0][n2][1] * dr[0][n2][1];
							dr[0][n2][2] = (ptsZZ[curBase1 + curBase3] - ptsZZ[curBase1])*l2;
							dn2 += dr[0][n2][2] * dr[0][n2][2];

							if (start[n] != iend[n]) 	fac = (c[n] - start[n])*ftn;
							else						fac = 0.0;
							id[n] = 0; id[n1] = 1;  id[n2] = 0;
							if (blend[n1] != 222)
							{
								if (dn1 == 0.0) cout << "By C version---dn1=it1..:dn1 " << dn1 << std::endl;
							} //删除部分代码
							unsigned base3 = (ndx - 1) * 24 + id[2] * 12 + id[1] * 6 + id[0] * 3;
							for (unsigned ri = 0; ri < 3; ri++) // do 311
							{
								rb[base3 + rs[ri]] = dn1 * dr[0][n1][rs[ri]]; //rb( 3 , 0:1 , 0:1 , 0:1 , dmrb )
							}// ends do 311
							id[n] = 0;  id[n1] = 0;  id[n2] = 1;
							if (blend[n2] != 222)
								if (dn2 == 0.0) cout << "By C version---dn2=it2..:dn2 " << dn2 << std::endl;
							for (unsigned ri = 0; ri < nrs; ri++) // do 321
							{
								rb[base3 + rs[ri]] = dn2 * dr[0][n2][rs[ri]]; //rb( 3 , 0:1 , 0:1 , 0:1 , dmrb )
							}// ends do 321		
							id[n] = 0;  id[n1] = 1;  id[n2] = 1;
							for (unsigned ri = 0; ri < 3; ri++) // do 331
							{
								rb[base3 + ri] = 0.0; //rb( 3 , 0:1 , 0:1 , 0:1 , dmrb )
							}// ends do 331	
						}//ends do 340
					} //ends do 340
				} //ends do 340
			} //ends if 
		}//ends do 350
	}//ends if
	/////////////////////////////////  Above Fotran File Line 760  /////////////////////////////////
	if (dimold == 3 && prof[0] * prof[1] * prof[2] != 0) //JavaScipt code中是0
	{
		for (unsigned c1 = start[0]; c1 <= iend[0]; c1 += iend[0] - start[0]) {  //do 410
			for (unsigned c2 = start[1]; c2 <= iend[1]; c2 += iend[1] - start[1]) {
				for (unsigned c3 = start[2]; c3 <= iend[2]; c3 += iend[2] - start[2]) {
					ndx = ndxrb[(c3 - 1)*JKMulDim + (c2 - 1)*dim1 + c1 - 1];
					for (unsigned ri = 0; ri < nrs; ri++) //do 410					
						rb[(ndx - 1) * 24 + 21 + rs[ri]] = 0.0;
				}
			}
		}	//ends do 410
	}//ends if 

	for (unsigned nn = 0; nn < dimold; nn++) //do 121
	{
		unsigned tmpN = proj[nn];
		if (start[tmpN] == iend[tmpN]) cout << "By C version---di:n " << tmpN << " " << start[tmpN] << std::endl;
		double di = 1.0 / float(iend[tmpN] - start[tmpN]);
		for (unsigned i = start[tmpN]; i <= iend[tmpN]; i++) //do 121
			fi[tmpN*dimax + i - 1] = (i - start[tmpN])*di;   //fi( dimax , 3 )		
	} //ends do 121
	/////////////////////////////////  Above Fotran File Line 787  /////////////////////////////////
	double fac1, ofac1, fac2, ofac2;
	for (unsigned nn = 0; nn < dimold; nn++) //do 125
	{
		unsigned n = proj[nn];
		unsigned n1 = cy[n][0];
		unsigned n2 = cy[n][1];
		unsigned ndir = 0;
		if (start[n1] != iend[n1])  ++ndir;
		if (start[n2] != iend[n2])  ++ndir;
		if (blend[n] == 222)
		{
			if (ndir == 2)
			{
				for (unsigned i1 = start[n1]; i1 <= iend[n1]; i1++) //do 123
				{
					c[n1] = i1;
					fac1 = fi[n1*dimax + i1 - 1];   ofac1 = 1.0 - fac1;
					for (unsigned i2 = start[n2]; i2 <= iend[n2]; i2++) //do 123
					{
						c[n2] = i2;
						fac2 = fi[n2*dimax + i2 - 1];   ofac2 = 1.0 - fac2;
						if (n == 0) {
							unsigned fBase = (c[2] - 1)*JKMulDim * 3 + (c[1] - 1)*dim1 * 3;
							unsigned aBase1 = dimax * 12;
							unsigned aBase2 = (i1 - 1) * 12;
							unsigned aBase3 = (i2 - 1) * 12;
							for (unsigned i = start[n] - 1; i < iend[n]; i++) //do 1231
							{  //ars( 3 ,   2, 0:1 , dimax , dimax );  arc(3,0:1,0:1,dimax)								
								facin[fBase + i * 3] = ofac1 * ars[i*aBase1 + aBase3] + fac1 * ars[i*aBase1 + aBase3 + 6]
									+ ofac2 * ars[i*aBase1 + aBase2 + 3] + fac2 * ars[i*aBase1 + aBase2 + 9]
									- ofac1 * ofac2*arc[i * 12] - ofac1 * fac2*arc[i * 12 + 6]
									- fac1 * ofac2*arc[i * 12 + 3] - fac1 * fac2*arc[i * 12 + 9];
							}//ends do 1231
						}// if n==0
						else if (n == 1) {
							unsigned fBase = (c[2] - 1)*JKMulDim * 3 + (c[0] - 1) * 3 + n;
							unsigned aBase1 = dimax * 12;
							unsigned aBase2 = (i1 - 1) * 12;
							unsigned aBase3 = (i2 - 1) * 12;
							for (unsigned i = start[n] - 1; i < iend[n]; i++) //do 1232
							{//ars( 3 ,   2, 0:1 , dimax , dimax )							
								facin[fBase + i * dim1 * 3] = ofac1 * ars[i*aBase1 + aBase3 + 1] + fac1 * ars[i*aBase1 + aBase3 + 7]
									+ ofac2 * ars[i*aBase1 + aBase2 + 4] + fac2 * ars[i*aBase1 + aBase2 + 10]
									- ofac1 * ofac2*arc[i * 12 + 1] - ofac1 * fac2*arc[i * 12 + 7]
									- fac1 * ofac2*arc[i * 12 + 4] - fac1 * fac2*arc[i * 12 + 10];
							}//ends do 1232
						}// if n==1
						else if (n == 2) {
							unsigned fBase = (c[1] - 1)*dim1 * 3 + (c[0] - 1) * 3 + n;
							unsigned aBase1 = dimax * 12;
							unsigned aBase2 = (i1 - 1) * 12;
							unsigned aBase3 = (i2 - 1) * 12;
							for (unsigned i = start[n] - 1; i < iend[n]; i++) //do 1233
							{
								facin[i*JKMulDim * 3 + fBase] = ofac1 * ars[i*aBase1 + aBase3 + 2] + fac1 * ars[i*aBase1 + aBase3 + 8]
									+ ofac2 * ars[i*aBase1 + aBase2 + 5] + fac2 * ars[i*aBase1 + aBase2 + 11]
									- ofac1 * ofac2*arc[i * 12 + 2] - ofac1 * fac2*arc[i * 12 + 8]
									- fac1 * ofac2*arc[i * 12 + 5] - fac1 * fac2*arc[i * 12 + 11];
							}//ends do 1233
						}// if n==3
					}//ends do 123
				}//ends do 123
			}//ends do 123
			else
			{
				for (unsigned i1 = start[n1]; i1 <= iend[n1]; i1++) //do 223
				{
					c[n1] = i1;
					fac1 = fi[n1*dimax + c[n1] - 1];   ofac1 = 1.0 - fac1;
					for (unsigned i2 = start[n2]; i2 <= iend[n2]; i2++) //do 223
					{
						c[n2] = i2;
						fac2 = fi[n2*dimax + c[n2] - 1];   ofac2 = 1.0 - fac2;
						if (n == 0) {
							unsigned fBase = (c[2] - 1)*JKMulDim * 3 + (c[1] - 1)*dim1 * 3;
							for (unsigned i = start[n] - 1; i < iend[n]; i++) //do 2231
							{
								facin[fBase + i * 3] = ofac1 * ofac2*arc[i * 12] + ofac1 * fac2*arc[i * 12 + 6]
									+ fac1 * ofac2*arc[i * 12 + 3] + fac1 * fac2*arc[i * 12 + 9];
							}//ends do 2231
						}// if n==0
						else if (n == 1) {
							unsigned fBase = (c[2] - 1)*JKMulDim * 3 + (c[0] - 1) * 3 + n;
							for (unsigned i = start[n] - 1; i < iend[n]; i++) //do 2232
							{
								facin[fBase + i * dim1 * 3] = ofac1 * ofac2*arc[i * 12 + 1] + ofac1 * fac2*arc[i * 12 + 7]
									+ fac1 * ofac2*arc[i * 12 + 4] + fac1 * fac2*arc[i * 12 + 10];
							}//ends do 2232
						}// if n==1
						else if (n == 2) {
							unsigned fBase = (c[1] - 1)*dim1 * 3 + (c[0] - 1) * 3 + n;
							for (unsigned i = start[n] - 1; i < iend[n]; i++) //do 2233
							{
								facin[i*JKMulDim * 3 + fBase] = ofac1 * ofac2*arc[i * 12 + 2] + ofac1 * fac2*arc[i * 12 + 8]
									+ fac1 * ofac2*arc[i * 12 + 5] + fac1 * fac2*arc[i * 12 + 11];
							}//ends do 2233
						}// if n==2
					}//ends do 223
				}//ends do 223
			}//ends do 223			
		}//if(blend[n]==222)
		else if (blend[n] == 143190)
		{
			for (unsigned i1 = start[n1]; i1 <= iend[n1]; i1++) //do 124
			{
				c[n1] = i1;
				for (unsigned i2 = start[n2]; i2 <= iend[n2]; i2++) //do 124
				{
					c[n2] = i2;
					if (n == 0) {
						unsigned fBase = (c[2] - 1)*JKMulDim * 3 + (c[1] - 1)*dim1 * 3;
						for (unsigned i = start[n] - 1; i < iend[n]; i++) { //1241
							facin[fBase + i * 3] = fi[n*dimax + i];
						}//1241
					}
					else if (n == 1) {
						unsigned fBase = (c[2] - 1)*JKMulDim * 3 + (c[0] - 1) * 3 + n;
						for (unsigned i = start[n] - 1; i < iend[n]; i++) { //1242
							facin[fBase + i * dim1 * 3] = fi[n*dimax + i];
						}//1242
					}
					else if (n == 2) {
						unsigned fBase = (c[1] - 1)*dim1 * 3 + (c[0] - 1) * 3 + n;
						for (unsigned i = start[n] - 1; i < iend[n]; i++) { //1243
							facin[i*JKMulDim * 3 + fBase] = fi[n*dimax + i];
						}//1243
					} //ends if
				} //ends do 124
			}//ends do 124
		}//if( blend[n] == 143190)
	}//ends do 125
	/////////////////////////////////  Above Fotran File Line 887  /////////////////////////////////
	//////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////
	for (unsigned nn = 0; nn < dimold; nn++) //do 239
	{
		unsigned n = proj[nn];
		if (prof[n] == 0) //类似JavaScipt code中的0
		{
			double fac = 0.0;
			for (unsigned c3 = start[2] - 1; c3 < iend[2]; c3++) //do 231
			{
				unsigned base1 = c3 * JKMulDim * 3;
				for (unsigned c2 = start[1] - 1; c2 < iend[1]; c2++)
				{
					unsigned base2 = base1 + c2 * dim1 * 3;
					for (unsigned c1 = start[0] - 1; c1 < iend[0]; c1++) {  //blen(3,0:1,0:1,dim1,dim2,dim3)
						fac = facin[base2 + c1 * 3 + n];
						blen[base2 * 4 + c1 * 12 + n] = 1.0 - fac;
						blen[base2 * 4 + c1 * 12 + n + 6] = fac;
					}//ends do 231
				}
			}
		}
		else if (prof[n] == 1)
		{
			double fac = 0.0, facs = 0.0;
			for (unsigned c3 = start[2] - 1; c3 < iend[2]; c3++) //do 232
			{
				unsigned base1 = c3 * JKMulDim * 3;
				for (unsigned c2 = start[1] - 1; c2 < iend[1]; c2++)
				{
					unsigned base2 = base1 + c2 * dim1 * 3;
					for (unsigned c1 = start[0] - 1; c1 < iend[0]; c1++) {  //blen(3,0:1,0:1,dim1,dim2,dim3)
						fac = facin[base2 + c1 * 3 + n];
						facs = fac * fac;
						blen[4 * base2 + 12 * c1 + n] = 1.0 - facs;
						blen[4 * base2 + 12 * c1 + n + 6] = facs;
						blen[4 * base2 + 12 * c1 + n + 3] = fac - facs;
					}//ends do 232
				}
			}
		}
		else if (prof[n] == 2)
		{
			double fac = 0.0, facs = 0.0, fact = 0.0;
			for (unsigned c3 = start[2] - 1; c3 < iend[2]; c3++) //do 233
			{
				unsigned base1 = c3 * JKMulDim * 3;
				for (unsigned c2 = start[1] - 1; c2 < iend[1]; c2++)
				{
					unsigned base2 = base1 + c2 * dim1 * 3;
					for (unsigned c1 = start[0] - 1; c1 < iend[0]; c1++) {  //blen(3,0:1,0:1,dim1,dim2,dim3)
						fac = facin[base2 + c1 * 3 + n];
						facs = fac * fac;
						fact = fac * 2 - facs;
						blen[4 * base2 + 12 * c1 + n] = 1.0 - fact;
						blen[4 * base2 + 12 * c1 + n + 6] = fact;
						blen[4 * base2 + 12 * c1 + n + 9] = facs - fac;
					}//ends do 233
				}
			}
		}
		else if (prof[n] == 3)
		{
			double fac = 0.0, facs = 0.0, facc = 0.0, fact = 0.0;
			for (unsigned c3 = start[2] - 1; c3 < iend[2]; c3++) //do 234
			{
				unsigned base1 = c3 * JKMulDim * 3;
				for (unsigned c2 = start[1] - 1; c2 < iend[1]; c2++)
				{
					unsigned base2 = base1 + c2 * dim1 * 3;
					for (unsigned c1 = start[0] - 1; c1 < iend[0]; c1++) {  //blen(3,0:1,0:1,dim1,dim2,dim3)
						fac = facin[base2 + c1 * 3 + n];
						facs = fac * fac;
						facc = fac * facs;
						fact = facs * 3 - facc * 2;
						blen[4 * base2 + 12 * c1 + n] = 1.0 - fact;
						blen[4 * base2 + 12 * c1 + n + 6] = fact;
						blen[4 * base2 + 12 * c1 + n + 3] = fac - facs * 2 + facc;
						blen[4 * base2 + 12 * c1 + n + 9] = facc - facs;
					}//ends do 234
				}
			}
		}
	} //ends do 239

	//for(unsigned ri=0;ri<nrs;ri++)
	for (unsigned c3 = start[2]; c3 <= iend[2]; c3++) //do 500
		for (unsigned c2 = start[1]; c2 <= iend[1]; c2++)
		{
			unsigned cBase1 = (c3 - 1)*KDim*JDim + (c2 - 1)*JDim - 1;
			unsigned cBase2 = (c3 - 1)*ptsD3Sz + (c2 - 1)*ptsD2Sz;
			for (unsigned c1 = start[0]; c1 <= iend[0]; c1++) {
				if (itype[cBase1 + c1] != lfix)
				{
					ptsXX[cBase2 + (c1 - 1)] = 0.0;
					ptsYY[cBase2 + (c1 - 1)] = 0.0;
					ptsZZ[cBase2 + (c1 - 1)] = 0.0;
				}
			} //ends do 500
		}

	if (proj[2] != 99) {  //99为错误标志
		C1_Terp1(0, blen, rb, nrs, ndxrb, start, iend, rs, ptsXX, ptsYY, ptsZZ, itype, dim1, dim2, dim3, protyp, proj, prof, blend);
		C1_Terp1(1, blen, rb, nrs, ndxrb, start, iend, rs, ptsXX, ptsYY, ptsZZ, itype, dim1, dim2, dim3, protyp, proj, prof, blend);
		C1_Terp1(2, blen, rb, nrs, ndxrb, start, iend, rs, ptsXX, ptsYY, ptsZZ, itype, dim1, dim2, dim3, protyp, proj, prof, blend);
		C2_Terp2(0, 1, blen, rb, nrs, ndxrb, start, iend, rs, ptsXX, ptsYY, ptsZZ, itype, dim1, dim2, dim3, protyp, proj, prof, blend);
		C2_Terp2(0, 2, blen, rb, nrs, ndxrb, start, iend, rs, ptsXX, ptsYY, ptsZZ, itype, dim1, dim2, dim3, protyp, proj, prof, blend);
		C2_Terp2(1, 2, blen, rb, nrs, ndxrb, start, iend, rs, ptsXX, ptsYY, ptsZZ, itype, dim1, dim2, dim3, protyp, proj, prof, blend);
		C3_Terp3(blen, rb, nrs, ndxrb, start, iend, rs, ptsXX, ptsYY, ptsZZ, itype, dim1, dim2, dim3, protyp, proj, prof, blend);
	}
	else if (proj[1] != 99) {  //99为错误标志
		C1_Terp1(proj[0], blen, rb, nrs, ndxrb, start, iend, rs, ptsXX, ptsYY, ptsZZ, itype, dim1, dim2, dim3, protyp, proj, prof, blend);
		C1_Terp1(proj[1], blen, rb, nrs, ndxrb, start, iend, rs, ptsXX, ptsYY, ptsZZ, itype, dim1, dim2, dim3, protyp, proj, prof, blend);
		if (iend[n1] - start[n1] > iend[n2] - start[n2])
			C2_Terp2(proj[0], proj[1], blen, rb, nrs, ndxrb, start, iend, rs, ptsXX, ptsYY, ptsZZ, itype, dim1, dim2, dim3, protyp, proj, prof, blend);
		else
			C2_Terp2(proj[1], proj[0], blen, rb, nrs, ndxrb, start, iend, rs, ptsXX, ptsYY, ptsZZ, itype, dim1, dim2, dim3, protyp, proj, prof, blend);
	}
	else
		C1_Terp1(proj[0], blen, rb, nrs, ndxrb, start, iend, rs, ptsXX, ptsYY, ptsZZ, itype, dim1, dim2, dim3, protyp, proj, prof, blend);

	free(facin);	free(blen);		free(ndxrb);
	free(itype);	free(arc);	free(ars);      free(dar);
	free(tars);		free(arsi); free(rb);  free(fi);
}

void C1_Terp1(unsigned n, double* blen, double * rb, unsigned nrs, unsigned* ndxrb, int* start, int* iend, unsigned* rs,
	double* ptsXX, double* ptsYY, double* ptsZZ, int* itype, unsigned dim10, unsigned dim20, unsigned dim30,
	unsigned protyp, unsigned* proj, unsigned* prof, unsigned* blend)
{
	//f( 3 , 0:dim10 ,0:dim20 ,0:dim30 );
	unsigned lfix = 555;
	unsigned cy[3][2] = { 1, 2, 2, 0, 0, 1 }; //cy( 2,3 )
	unsigned lf[4] = { 0, 1, 1, 1 };   //lf( 0:3 )
	unsigned le[4][2][2] = { 0, 1, 0, 0, 0, 1, 0, 0,
		0, 1, 1, 1, 0, 1, 0, 1 };    //le( 2 , 0:1 , 0:3 )
	unsigned nf[3], ne[3];
	unsigned cc[3][2]; //cc( 0:1 , 3 )
	unsigned c[3]; //c(3)

	unsigned n1 = cy[n][0];
	unsigned n2 = cy[n][1];
	nf[n1] = 0;	ne[n1] = 0;
	nf[n2] = 0;	ne[n2] = 0;
	cc[n][0] = start[n];		cc[n][1] = iend[n];
	unsigned npro = prof[n];
	unsigned ndx;

	unsigned ptsD3Sz = dim10 * dim20; //points( 3 , 0:dim1 , 0:dim2 , 0:dim3 )
	unsigned ptsD2Sz = dim10;

	for (unsigned nf1 = 0; nf1 <= lf[npro]; nf1++) //do 100
	{
		nf[n] = nf1;
		for (unsigned ne1 = le[npro][nf1][0]; ne1 <= le[npro][nf1][1]; ne1++)//do 100
		{
			ne[n] = ne1;
			//for(unsigned ri=0;ri<nrs;ri++) //do 100
			{
				for (unsigned i1 = start[n1]; i1 <= iend[n1]; i1++) //do 100
				{
					c[n1] = i1;
					cc[n1][0] = i1;		cc[n1][1] = i1;
					for (unsigned i2 = start[n2]; i2 <= iend[n2]; i2++) //do 100
					{
						c[n2] = i2;
						cc[n2][0] = i2;		cc[n2][1] = i2;
						unsigned ut3 = cc[2][ne[2]] - 1, ut2 = cc[1][ne[1]] - 1, ut1 = cc[0][ne[0]] - 1;
						ndx = ndxrb[ut3*ptsD3Sz + ut2 * dim10 + ut1];

						double rbVX = rb[(ndx - 1) * 24 + nf[2] * 12 + nf[1] * 6 + nf[0] * 3];
						double rbVY = rb[(ndx - 1) * 24 + nf[2] * 12 + nf[1] * 6 + nf[0] * 3 + 1];
						double rbVZ = rb[(ndx - 1) * 24 + nf[2] * 12 + nf[1] * 6 + nf[0] * 3 + 2];
						if (n == 0)
						{
							unsigned cBase1 = (c[2] - 1)*ptsD3Sz * 12 + (c[1] - 1)*dim10 * 12 + ne1 * 6 + nf1 * 3;
							unsigned cBase2 = (c[2] - 1)*ptsD3Sz + (c[1] - 1)*ptsD2Sz;
							unsigned cBase3 = (c[2] - 1)*ptsD3Sz + (c[1] - 1)*dim10;
							for (unsigned i = start[n] - 1; i < iend[n]; i++) //do 101
								if (itype[cBase3 + i] != lfix)
								{
									double blenV = blen[cBase1 + i * 12];
									ptsXX[cBase2 + i] += blenV * rbVX;
									ptsYY[cBase2 + i] += blenV * rbVY;
									ptsZZ[cBase2 + i] += blenV * rbVZ;
								}
						} //if n==0
						else if (n == 1)
						{
							unsigned cBase1 = (c[2] - 1)*ptsD3Sz * 12 + (c[0] - 1) * 12 + ne1 * 6 + nf1 * 3 + 1;
							unsigned cBase2 = (c[2] - 1)*ptsD3Sz + c[0] - 1;
							unsigned cBase3 = (c[2] - 1)*ptsD3Sz + c[0] - 1;
							for (unsigned i = start[n] - 1; i < iend[n]; i++) //do 102
								if (itype[cBase3 + i * dim10] != lfix)
								{
									double blenV = blen[cBase1 + i * dim10 * 12];
									ptsXX[cBase2 + i * ptsD2Sz] += blenV * rbVX;
									ptsYY[cBase2 + i * ptsD2Sz] += blenV * rbVY;
									ptsZZ[cBase2 + i * ptsD2Sz] += blenV * rbVZ;
								}
						} //if n==1
						else if (n == 2)
						{
							unsigned cBase1 = (c[1] - 1)*dim10 * 12 + (c[0] - 1) * 12 + ne1 * 6 + nf1 * 3 + 2;
							unsigned cBase2 = (c[1] - 1)*ptsD2Sz + c[0] - 1;
							unsigned cBase3 = (c[1] - 1)*dim10 + c[0] - 1;
							for (unsigned i = start[n] - 1; i < iend[n]; i++) //do 103
								if (itype[i*ptsD3Sz + cBase3] != lfix)
								{
									double blenV = blen[cBase1 + i * ptsD3Sz * 12];
									ptsXX[i*ptsD3Sz + cBase2] += blenV * rbVX;
									ptsYY[i*ptsD3Sz + cBase2] += blenV * rbVY;
									ptsZZ[i*ptsD3Sz + cBase2] += blenV * rbVZ;
								}
						} //if n==2
					}//do 100
				}//do 100
			}//do 100
		}
	}//do 100
}
//c2方向超限插值
void C2_Terp2(unsigned n1, unsigned n2, double* blen, double * rb, unsigned nrs, unsigned* ndxrb, int* start, int* iend, unsigned* rs,
	double* ptsXX, double* ptsYY, double* ptsZZ, int* itype, unsigned dim10, unsigned dim20, unsigned dim30,
	unsigned protyp, unsigned* proj, unsigned* prof, unsigned* blend)
{
	unsigned lfix = 555;
	unsigned acy[3][3] = { 99, 2, 1, 2, 99, 0, 1, 0, 99 }; //acy( 3,3 )
	unsigned lf[4] = { 0, 1, 1, 1 };   //lf( 0:3 )
	unsigned le[4][2][2] = { 0, 1, 0, 0, 0, 1, 0, 0,
		0, 1, 1, 1, 0, 1, 0, 1 };    //le( 2 , 0:1 , 0:3 )

	unsigned nf[3], ne[3];
	unsigned cc[3][2]; //cc( 0:1 , 3 )

	unsigned c[3]; //c(3)
	unsigned n = acy[n2][n1]; //如果为-1；则出错；
	if (n == 99) cout << "Maybe error in C2_Terp2. The n is: ", n;
	cc[n1][0] = start[n1];		cc[n1][1] = iend[n1];
	cc[n2][0] = start[n2];		cc[n2][1] = iend[n2];
	nf[n] = 0;	ne[n] = 0;
	unsigned npro1 = prof[n1], npro2 = prof[n2];

	unsigned ptsD3Sz = dim10 * dim20; //points( 3 , 0:dim1 , 0:dim2 , 0:dim3 )
	unsigned ptsD2Sz = dim10;

	unsigned ndx;

	for (unsigned nf1 = 0; nf1 <= lf[npro1]; nf1++) //do 200
	{
		nf[n1] = nf1;
		for (unsigned ne1 = le[npro1][nf1][0]; ne1 <= le[npro1][nf1][1]; ne1++)//do 200
		{
			ne[n1] = ne1;
			for (unsigned nf2 = 0; nf2 <= lf[npro2]; nf2++) //do 200
			{
				nf[n2] = nf2;
				for (unsigned ne2 = le[npro2][nf2][0]; ne2 <= le[npro2][nf2][1]; ne2++)//do 200
				{
					ne[n2] = ne2;
					//for(unsigned ri=0;ri<nrs;ri++) //do 200
					{
						for (unsigned i = start[n]; i <= iend[n]; i++)
						{
							c[n] = i;  cc[n][0] = i;		cc[n][1] = i;
							unsigned ut3 = cc[2][ne[2]] - 1, ut2 = cc[1][ne[1]] - 1, ut1 = cc[0][ne[0]] - 1;
							ndx = ndxrb[ut3*ptsD3Sz + ut2 * dim10 + ut1];
							for (unsigned i2 = start[n2]; i2 <= iend[n2]; i2++) //do 200
							{
								c[n2] = i2;
								double rbVX = rb[(ndx - 1) * 24 + nf[2] * 12 + nf[1] * 6 + nf[0] * 3];
								double rbVY = rb[(ndx - 1) * 24 + nf[2] * 12 + nf[1] * 6 + nf[0] * 3 + 1];
								double rbVZ = rb[(ndx - 1) * 24 + nf[2] * 12 + nf[1] * 6 + nf[0] * 3 + 2];
								if (n1 == 0)
								{
									unsigned base1 = (c[2] - 1)*ptsD3Sz + (c[1] - 1)*dim10;
									unsigned base2 = (c[2] - 1)*ptsD3Sz * 12 + (c[1] - 1)*dim10 * 12 + ne1 * 6 + nf1 * 3 + n1;
									unsigned base3 = (c[2] - 1)*ptsD3Sz * 12 + (c[1] - 1)*dim10 * 12 + ne2 * 6 + nf2 * 3 + n2;
									unsigned base4 = (c[2] - 1)*ptsD3Sz + (c[1] - 1)*ptsD2Sz;
									for (unsigned i1 = start[n1] - 1; i1 < iend[n1]; i1++) //do 201
										if (itype[base1 + i1] != lfix)
										{
											double blenV1 = blen[base2 + i1 * 12];
											double blenV2 = blen[base3 + i1 * 12];
											ptsXX[base4 + i1] -= blenV1 * blenV2 * rbVX;
											ptsYY[base4 + i1] -= blenV1 * blenV2 * rbVY;
											ptsZZ[base4 + i1] -= blenV1 * blenV2 * rbVZ;
										}
								} //if n==0
								else if (n1 == 1)
								{
									unsigned base1 = (c[2] - 1)*ptsD3Sz + c[0] - 1;
									unsigned base2 = (c[2] - 1)*ptsD3Sz * 12 + (c[0] - 1) * 12 + ne1 * 6 + nf1 * 3 + n1;
									unsigned base3 = (c[2] - 1)*ptsD3Sz * 12 + (c[0] - 1) * 12 + ne2 * 6 + nf2 * 3 + n2;
									unsigned base4 = (c[2] - 1)*ptsD3Sz + c[0] - 1;
									for (unsigned i1 = start[n1] - 1; i1 < iend[n1]; i1++) //do 202
										if (itype[base1 + i1 * dim10] != lfix)
										{
											double blenV1 = blen[base2 + i1 * dim10 * 12];
											double blenV2 = blen[base3 + i1 * dim10 * 12];
											ptsXX[base4 + i1 * ptsD2Sz] -= blenV1 * blenV2 * rbVX;
											ptsYY[base4 + i1 * ptsD2Sz] -= blenV1 * blenV2 * rbVY;
											ptsZZ[base4 + i1 * ptsD2Sz] -= blenV1 * blenV2 * rbVZ;
										}
								} //if n==1
								else if (n == 2)
								{
									unsigned base1 = (c[1] - 1)*dim10 + c[0] - 1;
									unsigned base2 = (c[1] - 1)*dim10 * 12 + (c[0] - 1) * 12 + ne1 * 6 + nf1 * 3 + n1;
									unsigned base3 = (c[1] - 1)*dim10 * 12 + (c[0] - 1) * 12 + ne2 * 6 + nf2 * 3 + n2;
									unsigned base4 = (c[1] - 1)*ptsD2Sz + c[0] - 1;
									for (unsigned i1 = start[n1] - 1; i1 < iend[n1]; i1++) //do 203
										if (itype[base1 + i1 * ptsD3Sz] != lfix)
										{
											double blenV1 = blen[i1*ptsD3Sz * 12 + base2];
											double blenV2 = blen[i1*ptsD3Sz * 12 + base3];
											ptsXX[base4 + i1 * ptsD3Sz] -= blenV1 * blenV2 * rbVX;
											ptsYY[base4 + i1 * ptsD3Sz] -= blenV1 * blenV2 * rbVY;
											ptsZZ[base4 + i1 * ptsD3Sz] -= blenV1 * blenV2 * rbVZ;
										}
								} //if n==2
							}//do 100
						}//do 100
					}
				}
			}//do 200
		}
	}//do 200
}
//c3方向超限插值
void C3_Terp3(double* blen, double * rb, unsigned nrs, unsigned* ndxrb, int* start, int* iend, unsigned* rs,
	double* ptsXX, double* ptsYY, double* ptsZZ, int* itype, unsigned dim10, unsigned dim20, unsigned dim30,
	unsigned protyp, unsigned* proj, unsigned* prof, unsigned* blend)
{
	unsigned lfix = 555;
	unsigned lf[4] = { 0, 1, 1, 1 };   //lf( 0:3 )
	unsigned le[4][2][2] = { 0, 1, 0, 0, 0, 1, 0, 0,
		0, 1, 1, 1, 0, 1, 0, 1 };    //le( 2 , 0:1 , 0:3 )
	unsigned cc[3][2]; //cc( 0:1 , 3 )
	cc[0][0] = start[0];		cc[0][1] = iend[0];
	cc[1][0] = start[1];		cc[1][1] = iend[1];
	cc[2][0] = start[2];		cc[2][1] = iend[2];
	unsigned npro1 = prof[0], npro2 = prof[1], npro3 = prof[2];
	unsigned ptsD3Sz = dim10 * dim20; 
	unsigned ptsD2Sz = dim10;
	unsigned ndx;
	for (unsigned nf1 = 0; nf1 <= lf[npro1]; nf1++) //do 300
		for (unsigned ne1 = le[npro1][nf1][0]; ne1 <= le[npro1][nf1][1]; ne1++)
			for (unsigned nf2 = 0; nf2 <= lf[npro2]; nf2++)
				for (unsigned ne2 = le[npro2][nf2][0]; ne2 <= le[npro2][nf2][1]; ne2++)
					for (unsigned nf3 = 0; nf3 <= lf[npro3]; nf3++)
						for (unsigned ne3 = le[npro3][nf3][0]; ne3 <= le[npro3][nf3][1]; ne3++)
						{
							unsigned ut3 = cc[2][ne3] - 1, ut2 = cc[1][ne2] - 1, ut1 = cc[0][ne1] - 1;
							ndx = ndxrb[ut3*ptsD3Sz + ut2 * dim10 + ut1];
							//for(unsigned ri=0;ri<nrs;ri++) //do 200
							{
								double rbVX = rb[(ndx - 1) * 24 + nf3 * 12 + nf2 * 6 + nf1 * 3];
								double rbVY = rb[(ndx - 1) * 24 + nf3 * 12 + nf2 * 6 + nf1 * 3 + 1];
								double rbVZ = rb[(ndx - 1) * 24 + nf3 * 12 + nf2 * 6 + nf1 * 3 + 2];
								for (unsigned c3 = start[2]; c3 <= iend[2]; c3++)
								{
									unsigned base1 = (c3 - 1)*ptsD3Sz;
									for (unsigned c2 = start[1]; c2 <= iend[1]; c2++)
									{
										unsigned base21 = base1 + (c2 - 1)*dim10;
										unsigned base31 = 12 * base21 + ne1 * 6 + nf1 * 3;
										unsigned base32 = 12 * base21 + ne2 * 6 + nf2 * 3;
										unsigned base33 = 12 * base21 + ne3 * 6 + nf3 * 3;
										unsigned base41 = (c3 - 1)*ptsD3Sz + (c2 - 1)*ptsD2Sz;
										for (unsigned c1 = start[0] - 1; c1 < iend[0]; c1++)
										{
											if (itype[base21 + c1] != lfix)
											{
												double blenV1 = blen[base31 + c1 * 12];
												double blenV2 = blen[base32 + c1 * 12 + 1];
												double blenV3 = blen[base33 + c1 * 12 + 2];
												ptsXX[base41 + c1] += blenV1 * blenV2*blenV3 * rbVX;
												ptsYY[base41 + c1] += blenV1 * blenV2*blenV3 * rbVY;
												ptsZZ[base41 + c1] += blenV1 * blenV2*blenV3 * rbVZ;
											}
										}
									}
								}
							}
						}
}

OpenMesh::VPropHandleT<MeshCube::Vertex3D> MeshCube::u_param;
OpenMesh::VPropHandleT<MeshCube::Vertex3D> MeshCube::v_param;

void MeshCube::MeshCoonsSurface(
	int udim, int vdim,
	std::vector<Vertex3D> &us,
	std::vector<Vertex3D> &ue,
	std::vector<Vertex3D> &vs,
	std::vector<Vertex3D> &ve,
	Vertex3D * dompts)
{
	double* ptsxx = nullptr;
	double* ptsyy = nullptr;
	double* ptszz = nullptr;
	int uv = udim * vdim;
	ptsxx = new double[uv] {};
	ptsyy = new double[uv] {};
	ptszz = new double[uv] {};

	int cubase = (vdim - 1)*udim;
	for (int iu = 0; iu < udim; ++iu)
	{
		ptsxx[iu] = us[iu][0];
		ptsyy[iu] = us[iu][1];
		ptszz[iu] = us[iu][2];
		ptsxx[cubase + iu] = ue[iu][0];
		ptsyy[cubase + iu] = ue[iu][1];
		ptszz[cubase + iu] = ue[iu][2];
	}
	int cvbase = udim - 1;
	for (int iv = 0; iv < vdim; ++iv)
	{
		ptsxx[iv*udim] = vs[iv][0];
		ptsyy[iv*udim] = vs[iv][1];
		ptszz[iv*udim] = vs[iv][2];

		ptsxx[cvbase + iv * udim] = ve[iv][0];
		ptsyy[cvbase + iv * udim] = ve[iv][1];
		ptszz[cvbase + iv * udim] = ve[iv][2];
	}
	int start[] = { 1, 1, 1 };
	int iend[] = { udim, vdim, 1 };
	int uvmax = udim < vdim ? vdim : udim;
	MeshCube meshf;
	meshf(ptsxx, ptsyy, ptszz, start, iend, udim, vdim, 1, uvmax);
	int u = udim;
	int v = vdim;
	for (int iv = 1; iv < v - 1; ++iv)
	{
		for (int iu = 1; iu < u - 1; ++iu)
		{
			int idx = iv * u + iu;
			dompts[idx][0] = ptsxx[idx];
			dompts[idx][1] = ptsyy[idx];
			dompts[idx][2] = ptszz[idx];
		}
	}
	for (int iu = 0; iu < u; ++iu)
	{
		float xval = us[iu][0];
		float yval = us[iu][1];
		float zval = us[iu][2];
		dompts[iu] = Vertex3D{ xval, yval, zval };
		dompts[iu][0] = xval;
		dompts[iu][1] = yval;
		dompts[iu][2] = zval;
		xval = ue[iu][0];
		yval = ue[iu][1];
		zval = ue[iu][2];
		dompts[(v - 1)*u + iu] = Vertex3D{ xval, yval, zval };
	}
	for (int iv = 0; iv < v; ++iv)
	{
		float xval = vs[iv][0];
		float yval = vs[iv][1];
		float zval = vs[iv][2];
		dompts[iv*u] = Vertex3D{ xval, yval, zval };

		xval = ve[iv][0];
		yval = ve[iv][1];
		zval = ve[iv][2];
		dompts[iv*u + cvbase] = Vertex3D{ xval, yval, zval };
	}
	delete[] ptsxx;
	delete[] ptsyy;
	delete[] ptszz;

}



//! 生成四边形网格
/*!
  \param  udim u方向维数
  \param  vdim v方向维数
  \param  us,ue u方向边界点数
  \param  vs,ve v方向边界点数
  \return
*/
void MeshCube::MeshPolyMesh(Polygon_mesh& msh,
	int udim, int vdim,
	std::vector<Vertex3D> &us,
	std::vector<Vertex3D> &ue,
	std::vector<Vertex3D> &vs,
	std::vector<Vertex3D> &ve, bool isappran)
{
	Vertex3D * dompts = new Vertex3D[udim*vdim]{};
	MeshCoonsSurface(udim, vdim, us, ue, vs, ve, dompts);
	vector<OpenMesh::VertexHandle> handle;
	handle.reserve(udim*vdim);
	if (isappran)
		msh.clean();
	//生成纹理坐标
	if (!msh.has_vertex_texcoords2D())
		msh.request_vertex_texcoords2D();
	for (int i = 0; i < udim*vdim; ++i)
	{
		auto vh = msh.add_vertex(dompts[i]);
		handle.push_back(vh);
		int iu = i % udim;
		int iv = i / udim;
		float u = static_cast<float>(iu) / (static_cast<float>(udim - 1));
		float v = static_cast<float>(iv) / (static_cast<float>(vdim - 1));
		OpenMesh::Vec2f vp{ u,v };
		msh.set_texcoord2D(vh, vp);
	}
	for (int iv = 0; iv < vdim - 1; ++iv)
	{
		for (int iu = 0; iu < udim - 1; ++iu)
		{
			int one = iv * udim + iu;
			int two = iv * udim + iu + 1;
			int the = iv * udim + iu + 1 + udim;
			int fou = iv * udim + iu + udim;
			vector<OpenMesh::VertexHandle> face;
			face.push_back(handle[one]);
			face.push_back(handle[fou]);
			face.push_back(handle[the]);
			face.push_back(handle[two]);
			msh.add_face(face);
		}
	}
	delete[] dompts;
}

//! 生成三角形网格
/*!
  \param  udim u方向维数
  \param  vdim v方向维数
  \param  us,ue u方向边界点数
  \param  vs,ve v方向边界点数
  \return
*/
void MeshCube::MeshTriangMesh(Triangle_mesh& msh,
	int udim, int vdim,
	std::vector<Vertex3D> &us,
	std::vector<Vertex3D> &ue,
	std::vector<Vertex3D> &vs,
	std::vector<Vertex3D> &ve,
	bool isappran, bool dir)
{
	Vertex3D * dompts = new Vertex3D[udim*vdim]{};
	MeshCoonsSurface(udim, vdim, us, ue, vs, ve, dompts);
	vector<OpenMesh::VertexHandle> handle;
	handle.reserve(udim*vdim);
	if (isappran)
		msh.clean();
	//生成纹理坐标
	if (!msh.has_vertex_texcoords2D())
		msh.request_vertex_texcoords2D();
	for (int i = 0; i < udim*vdim; ++i)
	{
		auto vh = msh.add_vertex(dompts[i]);
		handle.push_back(vh);
		int iu = i % udim;
		int iv = i / udim;
		float u = static_cast<float>(iu) / (static_cast<float>(udim - 1));
		float v = static_cast<float>(iv) / (static_cast<float>(vdim - 1));
		OpenMesh::Vec2f vp{ u,v };
		msh.set_texcoord2D(vh, vp);
	}
	for (int iv = 0; iv < vdim - 1; ++iv)
	{
		for (int iu = 0; iu < udim - 1; ++iu)
		{
			int one = iv * udim + iu;
			int two = iv * udim + iu + 1;
			int the = iv * udim + iu + 1 + udim;
			int fou = iv * udim + iu + udim;
			vector<OpenMesh::VertexHandle> face;
			if (dir)
			{
				face.push_back(handle[one]);
				face.push_back(handle[two]);
				face.push_back(handle[the]);
			}
			else
			{
				face.push_back(handle[one]);
				face.push_back(handle[the]);
				face.push_back(handle[two]);
			}
			msh.add_face(face);
			vector<OpenMesh::VertexHandle> sface;
			if (dir)
			{
				sface.push_back(handle[one]);
				sface.push_back(handle[the]);
				sface.push_back(handle[fou]);
			}
			else
			{
				sface.push_back(handle[one]);
				sface.push_back(handle[fou]);
				sface.push_back(handle[the]);
			}
			msh.add_face(sface);
		}
	}
	delete[] dompts;
}


void MeshCube::MeshTriangleParam(
	Triangle_mesh& msh, int udim, int vdim,
	std::vector<Vertex3D> &us,
	std::vector<Vertex3D> &ue,
	std::vector<Vertex3D> &vs,
	std::vector<Vertex3D> &ve,
	bool isappran , bool dir ,
	float umin , float umax,
	float vmin, float vmax ,
	int umaxp , int vmaxp)
{
	Vertex3D * dompts = new Vertex3D[udim*vdim]{};
	MeshCoonsSurface(udim, vdim, us, ue, vs, ve, dompts);
	vector<OpenMesh::VertexHandle> handle;
	handle.reserve(udim*vdim);
	if (isappran)
		msh.clean();
	//生成纹理坐标
	if (!u_param.is_valid() || !v_param.is_valid()) {
		msh.add_property(u_param);
		msh.add_property(v_param);
	}
	if (!msh.has_vertex_texcoords2D())
		msh.request_vertex_texcoords2D();
	for (int i = 0; i < udim*vdim; ++i)
	{
		auto vh = msh.add_vertex(dompts[i]);
		handle.push_back(vh);
		int iu = i % udim;
		int iv = i / udim;
		float u = static_cast<float>(iu) / (static_cast<float>(10));
		float v = static_cast<float>(iv) / (static_cast<float>(12));
		OpenMesh::Vec2f vp;
		vp[0] = u;
		vp[1] = v;
		//vp[0] = umin+u*(umax-umin);
		//vp[1] = vmin+v*(vmax-vmin);
		msh.set_texcoord2D(vh, vp);
	}
	for (int iv = 0; iv < vdim ; ++iv) {
		for (int iu = 0; iu < udim; ++iu) {
			VertexHandle vh(handle[iv*udim + iu]);
			int v1, v2;
			if (iv < vdim - 1)
				v1 = iv + 1;
			else
				v1 = iv;
			if (iu < udim - 1)
				v2 = iu + 1;
			else
				v2 = iu;
			auto vh0 = handle[(v1 - 1)*udim + v2 - 1];
			auto vh1 = handle[(v1 - 1)*udim + v2];
			auto vh2= handle[(v1)*udim + v2 - 1];
			msh.property(u_param, vh) = (msh.point(vh1) - msh.point(vh0)).normalize();
			msh.property(v_param, vh) = (msh.point(vh2) - msh.point(vh0)).normalize();
		}
	}
	for (int iv = 0; iv < vdim - 1; ++iv)
	{
		for (int iu = 0; iu < udim - 1; ++iu)
		{
			int one = iv * udim + iu;
			int two = iv * udim + iu + 1;
			int the = iv * udim + iu + 1 + udim;
			int fou = iv * udim + iu + udim;
			vector<OpenMesh::VertexHandle> face;
			vector<OpenMesh::VertexHandle> sface;
			if (dir)
			{
				face.push_back(handle[one]);
				face.push_back(handle[two]);
				face.push_back(handle[the]);
				sface.push_back(handle[one]);
				sface.push_back(handle[the]);
				sface.push_back(handle[fou]);
			}
			else
			{
				face.push_back(handle[one]);
				face.push_back(handle[the]);
				face.push_back(handle[two]);
				sface.push_back(handle[one]);
				sface.push_back(handle[fou]);
				sface.push_back(handle[the]);
			}
			msh.add_face(face);
			msh.add_face(sface);
		}
	}
	delete[] dompts;
}


void MeshCube::MeshTriangleMeshTeeth(Triangle_mesh& msh,
	int udim, int vdim,
	std::vector<Vertex3D> &us,
	std::vector<Vertex3D> &ue,
	std::vector<Vertex3D> &vs,
	std::vector<Vertex3D> &ve, bool dir ,
	float tus , float tue , float tvs, float tve)
{
	Vertex3D * dompts = new Vertex3D[udim*vdim]{};
	MeshCoonsSurface(udim, vdim, us, ue, vs, ve, dompts);
	vector<OpenMesh::VertexHandle> handle;
	handle.reserve(udim*vdim);
	//生成纹理坐标
	if (!msh.has_vertex_texcoords2D())
		msh.request_vertex_texcoords2D();
	for (int i = 0; i < udim*vdim; ++i)
	{
		auto vh = msh.add_vertex(dompts[i]);
		handle.push_back(vh);
		float iu = static_cast<float>(i % udim);
		float iv = static_cast<float>(i / udim);
		float u = iu/(static_cast<float>(udim - 1));
		float v = iv/(static_cast<float>(vdim - 1));
		OpenMesh::Vec2f vp;
		vp[0] = tus + u * (tue - tus);
		vp[1] = tvs + v * (tve - tvs);
		msh.set_texcoord2D(vh, vp);
	}
	for (int iv = 0; iv < vdim - 1; ++iv)
	{
		for (int iu = 0; iu < udim - 1; ++iu)
		{
			int one = iv * udim + iu;
			int two = iv * udim + iu + 1;
			int the = iv * udim + iu + 1 + udim;
			int fou = iv * udim + iu + udim;
			vector<OpenMesh::VertexHandle> face;
			if (dir)
			{
				face.push_back(handle[one]);
				face.push_back(handle[two]);
				face.push_back(handle[the]);
			}
			else
			{
				face.push_back(handle[one]);
				face.push_back(handle[the]);
				face.push_back(handle[two]);
			}
			msh.add_face(face);
			vector<OpenMesh::VertexHandle> sface;
			if (dir)
			{
				sface.push_back(handle[one]);
				sface.push_back(handle[the]);
				sface.push_back(handle[fou]);
			}
			else
			{
				sface.push_back(handle[one]);
				sface.push_back(handle[fou]);
				sface.push_back(handle[the]);
			}
			msh.add_face(sface);
		}
	}
	delete[] dompts;
}


void MeshCube::GetUVValueCir(int uvdim, double r,double z,
	std::vector<Vertex3D> &us,
	std::vector<Vertex3D> &ue,
	std::vector<Vertex3D> &vs,
	std::vector<Vertex3D> &ve)
{
	us.reserve(uvdim);
	ue.reserve(uvdim);
	vs.reserve(uvdim);
	ve.reserve(uvdim);
	for (int ik = 0; ik < uvdim; ++ik)
	{
		double ic = 1.0*ik / (double)(uvdim - 1);
		double the = M_PI / 2.0 * ic;
		double us_angle = 1.5*M_PI - the;
		double vs_angle = 1.5*M_PI + the;
		double ue_angle = the;
		double ve_angle = M_PI - the;

		us.push_back(r*Vertex3D(cos(us_angle), sin(us_angle), z));
		ue.push_back(r*Vertex3D(cos(ue_angle), sin(ue_angle), z));
		vs.push_back(r*Vertex3D(cos(vs_angle), sin(vs_angle), z));
		ve.push_back(r*Vertex3D(cos(ve_angle), sin(ve_angle), z));
	}
}
#pragma   warning(pop)