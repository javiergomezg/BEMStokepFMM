#include <math.h>
#include <stdio.h>
#include <iostream>
#include <omp.h>
#define REAL double

void calc_aux_cy(REAL *q , int qSize, 
                    REAL *xq0 , int xq0Size,
                    REAL *xq1 , int xq1Size,
                    REAL *xq2 , int xq2Size,
                    REAL *xi , int xiSize,
                    REAL *yi , int yiSize,
                    REAL *zi , int ziSize,
                    REAL *normal0 , int normal0Size, 
                    REAL *normal1 , int normal1Size, 
                    REAL *normal2 , int normal2Size,
                    int stype,
                    REAL *aux , int auxSize,
                    REAL E);

void calc_aux_cy(REAL *q , int qSize, 
                    REAL *xq0 , int xq0Size,
                    REAL *xq1 , int xq1Size,
                    REAL *xq2 , int xq2Size,
                    REAL *xi , int xiSize,
                    REAL *yi , int yiSize,
                    REAL *zi , int ziSize,
                    REAL *normal0 , int normal0Size, 
                    REAL *normal1 , int normal1Size, 
                    REAL *normal2 , int normal2Size,
                    int stype,
                    REAL *aux , int auxSize,
                    REAL E)
{
    REAL dx_pq[xiSize], dy_pq[yiSize], dz_pq[ziSize], R_pq[xiSize];

    for(int i=0; i<qSize; i++)
    {
        for (int j = 0; j < xiSize; j++)
        {
            dx_pq[j] = xi[j] - xq0[i];
            dy_pq[j] = yi[j] - xq1[i];
            dz_pq[j] = zi[j] - xq2[i];
        }

        for (int j = 0; j < xiSize; j++)
        {
            R_pq[j] = sqrt(dx_pq[j] * dx_pq[j] + dy_pq[j] * dy_pq[j] + dz_pq[j] * dz_pq[j]);
        }


        if (stype == 1)
        {
            for (int j = 0; j < auxSize; j++)
            {
                aux[j] = aux[j] - ( q[i] / ( R_pq[j] * R_pq[j] * R_pq[j] ) * (dx_pq[j] * normal0[j] + dy_pq[j] * normal1[j] + dz_pq[j] * normal2[j]) );
            }
        } else 
        {
            for (int j = 0; j < auxSize; j++)
            {
                aux[j] = aux[j] + q[i] / (E * R_pq[j]) ;
            }
        }
    }
};

