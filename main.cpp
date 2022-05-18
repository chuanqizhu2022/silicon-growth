#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <time.h>
#include <math.h>
#include <omp.h>
#include "CImg.h" //CImg ライブラリ（描画用）使用のためのヘッダ

using namespace cimg_library;
#define RND(x) ((double)(x) / RAND_MAX * rand()) //乱数の関数設定

#define NDX 128 //差分計算における計算領域一辺の分割数
#define NDY 128
#define NDZ 1

#define N 2
#define PI 3.141592

int ndmx = NDX - 1;
int ndmy = NDY - 1; //計算領域の一辺の差分分割数(差分ブロック数), ND-1を定義
int ndmz = NDZ - 1;
int nm = N - 1;
int mid = NDX / 4;
int th_num = 8;
int rows = NDX / th_num;

CImg<unsigned char> phi_fldxy(NDX, NDY, 1, 3);
char outFilePhi_xy[64];

int nstep, pstep;
double dtime, L, dx;
double gamma0, delta, mobi;
double astre, astrem;
double M0, W0, A0, S0, F0;
double Tm, dH;
double temp0, Tg, Tv, Tr;
int ni, nj;
int ix, iy, iz;
double r0, r, th0;
int xx0, yy0, zz0;
int intpos, frapass, curpos, curst, prepos, prest;
double sumplane;
int hasS, allS;
double intvel, inttemp;

double ****phi, ****phi2;
double ***alpha, ***beta;
double ***temp, ***temp2;
int ***phiNum, ****phiIdx;
double **aij, **wij, **mij;
double **anij, **thij, **vpij, **etaij;

//******* メインプログラム ******************************************
int main(int argc, char *argv[])
{
    phi = new double ***[N];
    phi2 = new double ***[N];
    for (ni = 0; ni <= nm; ni++)
    {
        phi[ni] = new double **[NDX];
        phi2[ni] = new double **[NDX];
        for (ix = 0; ix <= ndmx; ix++)
        {
            phi[ni][ix] = new double *[NDY];
            phi2[ni][ix] = new double *[NDY];
            for (iy = 0; iy <= ndmy; iy++)
            {
                phi[ni][ix][iy] = new double[NDZ];
                phi2[ni][ix][iy] = new double[NDZ];
            }
        }
    }

    phiIdx = new int ***[N + 1];
    for (ni = 0; ni <= N; ni++)
    {
        phiIdx[ni] = new int **[NDX];
        for (ix = 0; ix <= ndmx; ix++)
        {
            phiIdx[ni][ix] = new int *[NDY];
            for (iy = 0; iy <= ndmy; iy++)
            {
                phiIdx[ni][ix][iy] = new int[NDZ];
            }
        }
    }

    temp = new double **[NDX];
    temp2 = new double **[NDX];
    alpha = new double **[NDX];
    beta = new double **[NDX];
    for (ix = 0; ix <= ndmx; ix++)
    {
        temp[ix] = new double *[NDY];
        temp2[ix] = new double *[NDY];
        alpha[ix] = new double *[NDY];
        beta[ix] = new double *[NDY];
        for (iy = 0; iy <= ndmy; iy++)
        {
            temp[ix][iy] = new double[NDZ];
            temp2[ix][iy] = new double[NDZ];
            alpha[ix][iy] = new double[NDZ];
            beta[ix][iy] = new double[NDZ];
        }
    }

    phiNum = new int **[NDX];
    for (ix = 0; ix <= ndmx; ix++)
    {
        phiNum[ix] = new int *[NDY];

        for (iy = 0; iy <= ndmy; iy++)
        {
            phiNum[ix][iy] = new int[NDZ];
        }
    }

    aij = new double *[N];
    wij = new double *[N];
    mij = new double *[N];
    anij = new double *[N];
    thij = new double *[N];
    vpij = new double *[N];
    etaij = new double *[N];
    for (ni = 0; ni <= nm; ni++)
    {
        aij[ni] = new double[N];
        wij[ni] = new double[N];
        mij[ni] = new double[N];
        anij[ni] = new double[N];
        thij[ni] = new double[N];
        vpij[ni] = new double[N];
        etaij[ni] = new double[N];
    }

    nstep = 100001;
    pstep = 500;
    dx = 1.0e-5;
    dtime = 1.0e-3;
    delta = 5.0 * dx;
    mobi = 1.68e-9;
    gamma0 = 0.5;
    astre = 0.00;
    astrem = 0.0;
    Tm = 1687.0;
    // dH = 0.7e7;
    dH = 7.0e7;
    temp0 = 1686.0;
    Tg = 0000;
    Tv = 0.00005;
    Tr = Tg * Tv;

    A0 = 8.0 * delta * gamma0 / PI / PI;
    W0 = 4.0 * gamma0 / delta;
    M0 = mobi * PI * PI / (8.0 * delta);

    for (ni = 0; ni <= nm; ni++)
    {
        for (nj = 0; nj <= nm; nj++)
        {
            wij[ni][nj] = W0;
            aij[ni][nj] = A0;
            mij[ni][nj] = M0;
            anij[ni][nj] = 0;
            thij[ni][nj] = 0.0;
            vpij[ni][nj] = 0.0;
            etaij[ni][nj] = 0.0;
            if ((ni == 0) || (nj == 0))
            {
                anij[ni][nj] = 1;
            }
            if (ni == nj)
            {
                wij[ni][nj] = 0.0;
                aij[ni][nj] = 0.0;
                mij[ni][nj] = 0.0;
                anij[ni][nj] = 0;
            }
            if (ni != 0 && nj != 0)
            {
                mij[ni][nj] = 0.1 * mij[ni][nj];
            }
        }
    }

    for (ix = 0; ix <= ndmx; ix++)
    {
        for (iy = 0; iy <= ndmy; iy++)
        {
            for (iz = 0; iz <= ndmz; iz++)
            {
                if ((ix - NDX / 2) * (ix - NDX / 2) + (iy - NDY / 2) * (iy - NDY / 2) < (NDX / 4) * (NDX / 4))
                {
                    phi[1][ix][iy][iz] = 1.0;
                    phi[0][ix][iy][iz] = 0.0;
                    thij[1][0] = 0.0;
                    thij[0][1] = 0.0;
                }
                else
                {
                    phi[0][ix][iy][iz] = 1.0;
                }
            }
        }
    }

    for (ix = 0; ix <= ndmx; ix++)
    {
        for (iy = 0; iy <= ndmy; iy++)
        {
            for (iz = 0; iz <= ndmz; iz++)
            {
                temp[ix][iy][iz] = temp0 - NDX / 16 * dx * Tg + ix * dx * Tg;
            }
        }
    }
    // }

    std::cout << "the stability number is" << mobi * gamma0 * dtime / dx / dx << std::endl;

#pragma omp parallel num_threads(th_num)
    {
        int th_id, offset, start, end;

        th_id = omp_get_thread_num();
        offset = th_id * rows;
        start = offset;
        end = offset + rows - 1;

        int istep = 0;
        int i, j, k, ii, jj, kk;    //整数
        int ip, im, jp, jm, kp, km; //整数
        int n1, n2, n3, phinum;     //整数

        double th, vp, eta;
        double thii, vpii, etaii;
        double thetax, thetay;
        double epsilon0;
        double termiikk, termjjkk;
        double miijj;

        double phidx, phidy, phidz;
        double phidxx, phidyy, phidzz;
        double phidxy, phidxz, phidyz;
        double phiabs, phiabsii;
        double sum1, pddtt;

    start:;

        if (istep % pstep == 0 && th_id == 0)
        {
            cimg_forXY(phi_fldxy, x, y)
            {
                phi_fldxy(x, y, 0) = 255. * alpha[x][y][NDZ / 2] * (1.0 - beta[x][y][NDZ / 2]);         // red
                phi_fldxy(x, y, 1) = 255. * (1.0 - alpha[x][y][NDZ / 2]) * (1.0 - beta[x][y][NDZ / 2]); // green
                phi_fldxy(x, y, 2) = 128. * beta[x][y][NDZ / 2];                                        // blue
            }
            sprintf(outFilePhi_xy, "figures/phi/2dxy%d.png", istep);
            phi_fldxy.save_jpeg(outFilePhi_xy);
            std::cout << "-----------" << std::endl;
            std::cout << "the interface position is " << intpos << std::endl;
            std::cout << "the interface velocity is " << intvel << std::endl;
            std::cout << "the interface temperature is " << temp[intpos][0][0] << std::endl;
        }

        for (i = start; i <= end; i++)
        {
            for (j = 0; j <= ndmy; j++)
            {
                for (k = 0; k <= ndmz; k++)
                {
                    ip = i + 1;
                    im = i - 1;
                    jp = j + 1;
                    jm = j - 1;
                    kp = k + 1;
                    km = k - 1;
                    if (i == ndmx)
                    {
                        ip = ndmx;
                    }
                    if (i == 0)
                    {
                        im = 0;
                    }
                    if (j == ndmy)
                    {
                        jp = 0;
                    }
                    if (j == 0)
                    {
                        jm = ndmy;
                    }
                    if (k == ndmz)
                    {
                        kp = 0;
                    }
                    if (k == 0)
                    {
                        km = ndmz;
                    }
                    phinum = 0;
                    for (ii = 0; ii <= nm; ii++)
                    {
                        if ((phi[ii][i][j][k] > 0.0) ||
                            ((phi[ii][i][j][k] == 0.0) && (phi[ii][ip][j][k] > 0.0) ||
                             (phi[ii][im][j][k] > 0.0) ||
                             (phi[ii][i][jp][k] > 0.0) ||
                             (phi[ii][i][jm][k] > 0.0) ||
                             (phi[ii][i][j][kp] > 0.0) ||
                             (phi[ii][i][j][km] > 0.0)))
                        {
                            phinum++;
                            phiIdx[phinum][i][j][k] = ii;
                        }
                    }
                    phiNum[i][j][k] = phinum;
                }
            }
        }
#pragma omp barrier

        // Evolution Equations
        for (i = start; i <= end; i++)
        {
            for (j = 0; j <= ndmy; j++)
            {
                for (k = 0; k <= ndmz; k++)
                {
                    ip = i + 1;
                    im = i - 1;
                    jp = j + 1;
                    jm = j - 1;
                    kp = k + 1;
                    km = k - 1;
                    if (i == ndmx)
                    {
                        ip = ndmx;
                    }
                    if (i == 0)
                    {
                        im = 0;
                    }
                    if (j == ndmy)
                    {
                        jp = 0;
                    }
                    if (j == 0)
                    {
                        jm = ndmy;
                    }
                    if (k == ndmz)
                    {
                        kp = 0;
                    }
                    if (k == 0)
                    {
                        km = ndmz;
                    }

                    for (n1 = 1; n1 <= phiNum[i][j][k]; n1++)
                    {
                        ii = phiIdx[n1][i][j][k];
                        pddtt = 0.0;
                        for (n2 = 1; n2 <= phiNum[i][j][k]; n2++)
                        {
                            jj = phiIdx[n2][i][j][k];
                            sum1 = 0.0;
                            for (n3 = 1; n3 <= phiNum[i][j][k]; n3++)
                            {
                                kk = phiIdx[n3][i][j][k];

                                phidxx = (phi[kk][ip][j][k] + phi[kk][im][j][k] - 2.0 * phi[kk][i][j][k]) / dx / dx;
                                phidyy = (phi[kk][i][jp][k] + phi[kk][i][jm][k] - 2.0 * phi[kk][i][j][k]) / dx / dx;
                                phidzz = (phi[kk][i][j][kp] + phi[kk][i][j][km] - 2.0 * phi[kk][i][j][k]) / dx / dx;

                                termiikk = aij[ii][kk] * (phidxx + phidyy + phidzz);
                                termjjkk = aij[jj][kk] * (phidxx + phidyy + phidzz);

                                sum1 += 0.5 * (termiikk - termjjkk) + (wij[ii][kk] - wij[jj][kk]) * phi[kk][i][j][k];
                            }
                            if (ii > 0 && jj == 0)
                            {
                                F0 = -(temp[i][j][k] - Tm) * dH / Tm;
                            }
                            else if (ii == 0 && jj > 0)
                            {
                                F0 = (temp[i][j][k] - Tm) * dH / Tm;
                            }
                            else
                            {
                                F0 = 0.0;
                            }
                            pddtt += -2.0 * mij[ii][jj] / (double)phiNum[i][j][k] * (sum1 - 8.0 / PI * F0 * sqrt(phi[ii][i][j][k] * phi[jj][i][j][k]));
                        }
                        phi2[ii][i][j][k] = phi[ii][i][j][k] + pddtt * dtime;
                        if (phi2[ii][i][j][k] >= 1.0)
                        {
                            phi2[ii][i][j][k] = 1.0;
                        } //フェーズフィールドの変域補正
                        if (phi2[ii][i][j][k] <= 0.0)
                        {
                            phi2[ii][i][j][k] = 0.0;
                        }
                    }
                } // j
            }     // i
        }

        for (i = start; i <= end; i++)
        {
            for (j = 0; j <= ndmy; j++)
            {
                for (k = 0; k <= ndmz; k++)
                {
                    for (ii = 0; ii <= nm; ii++)
                    {
                        phi[ii][i][j][k] = phi2[ii][i][j][k];
                    }
                    temp[i][j][k] -= Tr * dtime;
                }
            }
        }

        //
        for (i = start; i <= end; i++)
        {
            for (j = 0; j <= ndmy; j++)
            {
                for (k = 0; k <= ndmz; k++)
                {
                    sum1 = 0.0;
                    for (ii = 0; ii <= nm; ii++)
                    {
                        sum1 += phi[ii][i][j][k];
                    }
                    for (ii = 0; ii <= nm; ii++)
                    {
                        phi[ii][i][j][k] = phi[ii][i][j][k] / sum1;
                    }
                }
            }
        }

        for (i = start; i <= end; i++)
        {
            for (j = 0; j <= ndmy; j++)
            {
                for (k = 0; k <= ndmz; k++)
                {
                    sum1 = 0.0;
                    for (ii = 1; ii <= nm; ii++)
                    {
                        sum1 += abs(thij[ii][0] * 4 / PI) * phi[ii][i][j][k];
                    }
                    alpha[i][j][k] = sum1;
                    beta[i][j][k] = phi[0][i][j][k];
                }
            }
        }

#pragma omp barrier

        istep = istep + 1;
        if (istep < nstep)
        {
            goto start;
        }

    end:;
    }
    return 0;
}
