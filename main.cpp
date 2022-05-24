#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <time.h>
#include <math.h>
#include <omp.h>
#include "CImg.h" //CImg ライブラリ（描画用）使用のためのヘッダ

using namespace cimg_library;
#define RND(x) ((double)(x) / RAND_MAX * rand()) //乱数の関数設定

#define NDX 64 //差分計算における計算領域一辺の分割数
#define NDY 64
#define NDZ 64

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
double astre, astrem, al0;
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
    astre = 0.2;
    al0 = 35.0 / 180.0 * PI;
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

    // thij[0][1] = PI / 4.0;
    // thij[1][0] = PI / 4.0;

    for (ix = 0; ix <= ndmx; ix++)
    {
        for (iy = 0; iy <= ndmy; iy++)
        {
            for (iz = 0; iz <= ndmz; iz++)
            {
                if ((ix - NDX / 2) * (ix - NDX / 2) + (iy - NDY / 2) * (iy - NDY / 2) + (iz - NDZ / 2) * (iz - NDZ / 2) < (NDX / 4) * (NDX / 4))
                {
                    phi[1][ix][iy][iz] = 1.0;
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

        double xxp, xyp, xzp, yxp, yyp, yzp, zxp, zyp, zzp;
        double I1, I2, I3;
        double cosap, ep;
        double dcosapdphix, dcosapdphiy, dcosapdphiz;
        double dphiabs2dx, dphiabs2dy, dphiabs2dz;
        double ddcosapdphixdx, ddcosapdphiydy, ddcosapdphizdz;
        double dcosapdx, dcosapdy, dcosapdz;

        double phidx, phidy, phidz;
        double phidxp, phidyp, phidzp;
        double phidxx, phidyy, phidzz;
        double phidxpx, phidypx, phidzpx;
        double phidxpy, phidypy, phidzpy;
        double phidxpz, phidypz, phidzpz;
        double phidxy, phidxz, phidyz;
        double phiabs2, phiabsii;
        double termx0, termx1, termx2, termx3;
        double termy0, termy1, termy2, termy3;
        double termz0, termz1, termz2, termz3;
        double termx, termy, termz;
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

            FILE *stream;
            char buffer[30];
            sprintf(buffer, "data/phi/3d%d.vtk", istep);
            stream = fopen(buffer, "a");

            fprintf(stream, "# vtk DataFile Version 1.0\n");
            fprintf(stream, "phi_%d.vtk\n", istep);
            fprintf(stream, "ASCII\n");
            fprintf(stream, "DATASET STRUCTURED_POINTS\n");
            fprintf(stream, "DIMENSIONS %d %d %d\n", NDX, NDY, NDZ);
            fprintf(stream, "ORIGIN 0.0 0.0 0.0\n");
            fprintf(stream, "ASPECT_RATIO 1.0 1.0 1.0\n");
            fprintf(stream, "\n");
            fprintf(stream, "POINT_DATA %d\n", NDX * NDY * NDZ);
            fprintf(stream, "SCALARS scalars float\n");
            fprintf(stream, "LOOKUP_TABLE default\n");

            for (k = 0; k <= ndmz; k++)
            {
                for (j = 0; j <= ndmy; j++)
                {
                    for (i = 0; i <= ndmx; i++)
                    {
                        fprintf(stream, "%e\n", phi[0][i][j][k]);
                    }
                }
            }
            fclose(stream);

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

                                phidx = (phi[kk][ip][j][k] - phi[kk][im][j][k]) / 2.0 / dx;
                                phidy = (phi[kk][i][jp][k] - phi[kk][i][jm][k]) / 2.0 / dx;
                                phidz = (phi[kk][i][j][kp] - phi[kk][i][j][km]) / 2.0 / dx;

                                phidxx = (phi[kk][ip][j][k] + phi[kk][im][j][k] - 2.0 * phi[kk][i][j][k]) / dx / dx;
                                phidyy = (phi[kk][i][jp][k] + phi[kk][i][jm][k] - 2.0 * phi[kk][i][j][k]) / dx / dx;
                                phidzz = (phi[kk][i][j][kp] + phi[kk][i][j][km] - 2.0 * phi[kk][i][j][k]) / dx / dx;

                                phidxy = (phi[kk][ip][jp][k] + phi[kk][im][jm][k] - phi[kk][im][jp][k] - phi[kk][ip][jm][k]) / 4.0 / dx / dx;
                                phidxz = (phi[kk][ip][j][kp] + phi[kk][im][j][km] - phi[kk][im][j][kp] - phi[kk][ip][j][km]) / 4.0 / dx / dx;
                                phidyz = (phi[kk][i][jp][kp] + phi[kk][i][jm][km] - phi[kk][i][jm][kp] - phi[kk][i][jp][km]) / 4.0 / dx / dx;

                                phiabs2 = phidx * phidx + phidy * phidy + phidz * phidz;

                                if (anij[ii][kk] == 1)
                                {
                                    epsilon0 = sqrt(aij[ii][kk]);

                                    th = thij[ii][kk];
                                    vp = vpij[ii][kk];
                                    eta = etaij[ii][kk];

                                    xxp = cos(th) * cos(vp);
                                    yxp = sin(th) * cos(vp);
                                    zxp = sin(vp);
                                    xyp = -sin(th) * cos(eta) - cos(th) * sin(vp) * sin(eta);
                                    yyp = cos(th) * cos(eta) - sin(th) * sin(vp) * sin(eta);
                                    zyp = cos(vp) * sin(eta);
                                    xzp = sin(eta) * sin(th) - cos(eta) * cos(th) * sin(vp);
                                    yzp = -sin(eta) * cos(th) - cos(eta) * sin(th) * sin(vp);
                                    zzp = cos(eta) * cos(vp);

                                    I1 = xxp + yxp + zxp;
                                    I2 = xyp + yyp + zyp;
                                    I3 = xzp + yzp + zzp;

                                    phidxp = phidx * xxp + phidy * yxp + phidz * zxp;
                                    phidyp = phidx * xyp + phidy * yyp + phidz * zyp;
                                    phidzp = phidx * xzp + phidy * yzp + phidz * zzp;

                                    phidxpx = phidxx * xxp + phidxy * yxp + phidxz * zxp;
                                    phidypx = phidxx * xyp + phidxy * yyp + phidxz * zyp;
                                    phidzpx = phidxx * xzp + phidxy * yzp + phidxz * zzp;

                                    phidxpy = phidxy * xxp + phidyy * yxp + phidyz * zxp;
                                    phidypy = phidxy * xyp + phidyy * yyp + phidyz * zyp;
                                    phidzpy = phidxy * xzp + phidyy * yzp + phidyz * zzp;

                                    phidxpz = phidxz * xxp + phidyz * yxp + phidzz * zxp;
                                    phidypz = phidxz * xyp + phidyz * yyp + phidzz * zyp;
                                    phidzpz = phidxz * xzp + phidyz * yzp + phidzz * zzp;

                                    cosap = sqrt(abs(pow((I1 * abs(phidxp) + I2 * abs(phidyp) + I3 * abs(phidzp)) / sqrt(3.0 * phiabs2), 2.0) - 1.0e-8));

                                    if (cosap > cos(al0))
                                    {
                                        ep = epsilon0 * (1.0 + astre * (abs(cosap) + tan(al0) * sqrt(1.0 - cosap * cosap) - 1.0 / cos(al0)));

                                        dcosapdphix = (I1 * phidxp / abs(phidxp) * xxp - sqrt(3.0) * phidx / sqrt(phiabs2) * cosap) / sqrt(3.0 * phiabs2);
                                        dcosapdphiy = (I2 * phidyp / abs(phidyp) * yyp - sqrt(3.0) * phidy / sqrt(phiabs2) * cosap) / sqrt(3.0 * phiabs2);
                                        dcosapdphiz = (I3 * phidzp / abs(phidzp) * zzp - sqrt(3.0) * phidz / sqrt(phiabs2) * cosap) / sqrt(3.0 * phiabs2);

                                        dcosapdx = (I1 * phidxp / abs(phidxp) * phidxpx + I2 * phidyp / abs(phidyp) * phidypx + I3 * phidzp / abs(phidzp) * phidzpx - sqrt(3.0) * (phidx * phidxx + phidy * phidxy + phidz * phidxz) / sqrt(phiabs2) * cosap) / sqrt(3.0 * phiabs2);
                                        dcosapdy = (I1 * phidxp / abs(phidxp) * phidxpy + I2 * phidyp / abs(phidyp) * phidypy + I3 * phidzp / abs(phidzp) * phidzpy - sqrt(3.0) * (phidx * phidxy + phidy * phidyy + phidz * phidyz) / sqrt(phiabs2) * cosap) / sqrt(3.0 * phiabs2);
                                        dcosapdz = (I1 * phidxp / abs(phidxp) * phidxpz + I2 * phidyp / abs(phidyp) * phidypz + I3 * phidzp / abs(phidzp) * phidzpz - sqrt(3.0) * (phidx * phidxz + phidy * phidyz + phidz * phidzz) / sqrt(phiabs2) * cosap) / sqrt(3.0 * phiabs2);

                                        dphiabs2dx = 2.0 * phidx * phidxx + 2.0 * phidy * phidxy + 2.0 * phidz * phidxz;
                                        dphiabs2dy = 2.0 * phidx * phidxy + 2.0 * phidy * phidyy + 2.0 * phidz * phidyz;
                                        dphiabs2dz = 2.0 * phidx * phidxz + 2.0 * phidy * phidyz + 2.0 * phidz * phidzz;

                                        ddcosapdphixdx = I1 * xxp / (3.0 * phiabs2 * phidxp * phidxp) * (phidxpx * sqrt(3.0 * phiabs2 * phidxp * phidxp) - 1.0 / sqrt(3.0 * phiabs2 * phidxp * phidxp) * (1.5 * dphiabs2dx * phidxp * phidxp + 3.0 * phidxp * phidxpx * phiabs2)) - sqrt(3.0) / pow(phiabs2, 2.0) * ((phidxx * cosap + phidx * dcosapdx) * phiabs2 - dphiabs2dx * phidx * cosap);
                                        ddcosapdphiydy = I2 * yyp / (3.0 * phiabs2 * phidyp * phidyp) * (phidypy * sqrt(3.0 * phiabs2 * phidyp * phidyp) - 1.0 / sqrt(3.0 * phiabs2 * phidyp * phidyp) * (1.5 * dphiabs2dy * phidyp * phidyp + 3.0 * phidyp * phidypy * phiabs2)) - sqrt(3.0) / pow(phiabs2, 2.0) * ((phidyy * cosap + phidy * dcosapdy) * phiabs2 - dphiabs2dy * phidy * cosap);
                                        ddcosapdphizdz = I3 * zzp / (3.0 * phiabs2 * phidzp * phidzp) * (phidzpz * sqrt(3.0 * phiabs2 * phidzp * phidzp) - 1.0 / sqrt(3.0 * phiabs2 * phidzp * phidzp) * (1.5 * dphiabs2dz * phidzp * phidzp + 3.0 * phidzp * phidzpz * phiabs2)) - sqrt(3.0) / pow(phiabs2, 2.0) * ((phidzz * cosap + phidz * dcosapdz) * phiabs2 - dphiabs2dz * phidz * cosap);

                                        termx0 = 2.0 * ep * dcosapdx * epsilon0 * astre * (1.0 - tan(al0) * abs(cosap) / sqrt(1.0 - cosap * cosap)) * phidx + phidxx * ep * ep;
                                        termy0 = 2.0 * ep * dcosapdy * epsilon0 * astre * (1.0 - tan(al0) * abs(cosap) / sqrt(1.0 - cosap * cosap)) * phidy + phidyy * ep * ep;
                                        termz0 = 2.0 * ep * dcosapdz * epsilon0 * astre * (1.0 - tan(al0) * abs(cosap) / sqrt(1.0 - cosap * cosap)) * phidz + phidzz * ep * ep;

                                        termx1 = dcosapdx * dcosapdphix * pow(epsilon0 * astre * (1.0 - tan(al0) * abs(cosap) / sqrt(1.0 - cosap * cosap)), 2.0) * phiabs2;
                                        termy1 = dcosapdy * dcosapdphiy * pow(epsilon0 * astre * (1.0 - tan(al0) * abs(cosap) / sqrt(1.0 - cosap * cosap)), 2.0) * phiabs2;
                                        termz1 = dcosapdz * dcosapdphiz * pow(epsilon0 * astre * (1.0 - tan(al0) * abs(cosap) / sqrt(1.0 - cosap * cosap)), 2.0) * phiabs2;

                                        termx2 = ep * epsilon0 * astre * dcosapdphix * (-tan(al0)) / pow(1.0 - cosap * cosap, 1.5) * dcosapdx * phiabs2 + ep * epsilon0 * astre * ddcosapdphixdx * (1.0 - tan(al0) * cosap / sqrt(1.0 - cosap * cosap)) * phiabs2;
                                        termy2 = ep * epsilon0 * astre * dcosapdphiy * (-tan(al0)) / pow(1.0 - cosap * cosap, 1.5) * dcosapdy * phiabs2 + ep * epsilon0 * astre * ddcosapdphiydy * (1.0 - tan(al0) * cosap / sqrt(1.0 - cosap * cosap)) * phiabs2;
                                        termz2 = ep * epsilon0 * astre * dcosapdphiz * (-tan(al0)) / pow(1.0 - cosap * cosap, 1.5) * dcosapdz * phiabs2 + ep * epsilon0 * astre * ddcosapdphizdz * (1.0 - tan(al0) * cosap / sqrt(1.0 - cosap * cosap)) * phiabs2;

                                        termx3 = ep * epsilon0 * astre * dcosapdphix * (1.0 - tan(al0) * abs(cosap) / sqrt(1.0 - cosap * cosap)) * dphiabs2dx;
                                        termy3 = ep * epsilon0 * astre * dcosapdphiy * (1.0 - tan(al0) * abs(cosap) / sqrt(1.0 - cosap * cosap)) * dphiabs2dy;
                                        termz3 = ep * epsilon0 * astre * dcosapdphiz * (1.0 - tan(al0) * abs(cosap) / sqrt(1.0 - cosap * cosap)) * dphiabs2dz;

                                        termx = termx0 + termx1 + termx2 + termx3;
                                        termy = termy0 + termy1 + termy2 + termy3;
                                        termz = termz0 + termz1 + termz2 + termz3;

                                        termiikk = termx + termy + termz;
                                    }
                                    else
                                    {
                                        termiikk = aij[ii][kk] * (phidxx + phidyy + phidzz);
                                    }
                                }
                                else
                                {
                                    termiikk = aij[ii][kk] * (phidxx + phidyy + phidzz);
                                }

                                if (anij[jj][kk] == 1)
                                {
                                    epsilon0 = sqrt(aij[jj][kk]);

                                    th = thij[jj][kk];
                                    vp = vpij[jj][kk];
                                    eta = etaij[jj][kk];

                                    xxp = cos(th) * cos(vp);
                                    yxp = sin(th) * cos(vp);
                                    zxp = sin(vp);
                                    xyp = -sin(th) * cos(eta) - cos(th) * sin(vp) * sin(eta);
                                    yyp = cos(th) * cos(eta) - sin(th) * sin(vp) * sin(eta);
                                    zyp = cos(vp) * sin(eta);
                                    xzp = sin(eta) * sin(th) - cos(eta) * cos(th) * sin(vp);
                                    yzp = -sin(eta) * cos(th) - cos(eta) * sin(th) * sin(vp);
                                    zzp = cos(eta) * cos(vp);

                                    I1 = xxp + yxp + zxp;
                                    I2 = xyp + yyp + zyp;
                                    I3 = xzp + yzp + zzp;

                                    phidxp = phidx * xxp + phidy * yxp + phidz * zxp;
                                    phidyp = phidx * xyp + phidy * yyp + phidz * zyp;
                                    phidzp = phidx * xzp + phidy * yzp + phidz * zzp;

                                    phidxpx = phidxx * xxp + phidxy * yxp + phidxz * zxp;
                                    phidypx = phidxx * xyp + phidxy * yyp + phidxz * zyp;
                                    phidzpx = phidxx * xzp + phidxy * yzp + phidxz * zzp;

                                    phidxpy = phidxy * xxp + phidyy * yxp + phidyz * zxp;
                                    phidypy = phidxy * xyp + phidyy * yyp + phidyz * zyp;
                                    phidzpy = phidxy * xzp + phidyy * yzp + phidyz * zzp;

                                    phidxpz = phidxz * xxp + phidyz * yxp + phidzz * zxp;
                                    phidypz = phidxz * xyp + phidyz * yyp + phidzz * zyp;
                                    phidzpz = phidxz * xzp + phidyz * yzp + phidzz * zzp;

                                    cosap = sqrt(abs(pow((I1 * abs(phidxp) + I2 * abs(phidyp) + I3 * abs(phidzp)) / sqrt(3.0 * phiabs2), 2.0) - 1.0e-8));

                                    if (cosap > cos(al0))
                                    {
                                        ep = epsilon0 * (1.0 + astre * (abs(cosap) + tan(al0) * sqrt(1.0 - cosap * cosap) - 1.0 / cos(al0)));

                                        dcosapdphix = (I1 * phidxp / abs(phidxp) * xxp - sqrt(3.0) * phidx / sqrt(phiabs2) * cosap) / sqrt(3.0 * phiabs2);
                                        dcosapdphiy = (I2 * phidyp / abs(phidyp) * yyp - sqrt(3.0) * phidy / sqrt(phiabs2) * cosap) / sqrt(3.0 * phiabs2);
                                        dcosapdphiz = (I3 * phidzp / abs(phidzp) * zzp - sqrt(3.0) * phidz / sqrt(phiabs2) * cosap) / sqrt(3.0 * phiabs2);

                                        dcosapdx = (I1 * phidxp / abs(phidxp) * phidxpx + I2 * phidyp / abs(phidyp) * phidypx + I3 * phidzp / abs(phidzp) * phidzpx - sqrt(3.0) * (phidx * phidxx + phidy * phidxy + phidz * phidxz) / sqrt(phiabs2) * cosap) / sqrt(3.0 * phiabs2);
                                        dcosapdy = (I1 * phidxp / abs(phidxp) * phidxpy + I2 * phidyp / abs(phidyp) * phidypy + I3 * phidzp / abs(phidzp) * phidzpy - sqrt(3.0) * (phidx * phidxy + phidy * phidyy + phidz * phidyz) / sqrt(phiabs2) * cosap) / sqrt(3.0 * phiabs2);
                                        dcosapdz = (I1 * phidxp / abs(phidxp) * phidxpz + I2 * phidyp / abs(phidyp) * phidypz + I3 * phidzp / abs(phidzp) * phidzpz - sqrt(3.0) * (phidx * phidxz + phidy * phidyz + phidz * phidzz) / sqrt(phiabs2) * cosap) / sqrt(3.0 * phiabs2);

                                        dphiabs2dx = 2.0 * phidx * phidxx + 2.0 * phidy * phidxy + 2.0 * phidz * phidxz;
                                        dphiabs2dy = 2.0 * phidx * phidxy + 2.0 * phidy * phidyy + 2.0 * phidz * phidyz;
                                        dphiabs2dz = 2.0 * phidx * phidxz + 2.0 * phidy * phidyz + 2.0 * phidz * phidzz;

                                        ddcosapdphixdx = I1 * xxp / (3.0 * phiabs2 * phidxp * phidxp) * (phidxpx * sqrt(3.0 * phiabs2 * phidxp * phidxp) - 1.0 / sqrt(3.0 * phiabs2 * phidxp * phidxp) * (1.5 * dphiabs2dx * phidxp * phidxp + 3.0 * phidxp * phidxpx * phiabs2)) - sqrt(3.0) / pow(phiabs2, 2.0) * ((phidxx * cosap + phidx * dcosapdx) * phiabs2 - dphiabs2dx * phidx * cosap);
                                        ddcosapdphiydy = I2 * yyp / (3.0 * phiabs2 * phidyp * phidyp) * (phidypy * sqrt(3.0 * phiabs2 * phidyp * phidyp) - 1.0 / sqrt(3.0 * phiabs2 * phidyp * phidyp) * (1.5 * dphiabs2dy * phidyp * phidyp + 3.0 * phidyp * phidypy * phiabs2)) - sqrt(3.0) / pow(phiabs2, 2.0) * ((phidyy * cosap + phidy * dcosapdy) * phiabs2 - dphiabs2dy * phidy * cosap);
                                        ddcosapdphizdz = I3 * zzp / (3.0 * phiabs2 * phidzp * phidzp) * (phidzpz * sqrt(3.0 * phiabs2 * phidzp * phidzp) - 1.0 / sqrt(3.0 * phiabs2 * phidzp * phidzp) * (1.5 * dphiabs2dz * phidzp * phidzp + 3.0 * phidzp * phidzpz * phiabs2)) - sqrt(3.0) / pow(phiabs2, 2.0) * ((phidzz * cosap + phidz * dcosapdz) * phiabs2 - dphiabs2dz * phidz * cosap);

                                        termx0 = 2.0 * ep * dcosapdx * epsilon0 * astre * (1.0 - tan(al0) * abs(cosap) / sqrt(1.0 - cosap * cosap)) * phidx + phidxx * ep * ep;
                                        termy0 = 2.0 * ep * dcosapdy * epsilon0 * astre * (1.0 - tan(al0) * abs(cosap) / sqrt(1.0 - cosap * cosap)) * phidy + phidyy * ep * ep;
                                        termz0 = 2.0 * ep * dcosapdz * epsilon0 * astre * (1.0 - tan(al0) * abs(cosap) / sqrt(1.0 - cosap * cosap)) * phidz + phidzz * ep * ep;

                                        termx1 = dcosapdx * dcosapdphix * pow(epsilon0 * astre * (1.0 - tan(al0) * abs(cosap) / sqrt(1.0 - cosap * cosap)), 2.0) * phiabs2;
                                        termy1 = dcosapdy * dcosapdphiy * pow(epsilon0 * astre * (1.0 - tan(al0) * abs(cosap) / sqrt(1.0 - cosap * cosap)), 2.0) * phiabs2;
                                        termz1 = dcosapdz * dcosapdphiz * pow(epsilon0 * astre * (1.0 - tan(al0) * abs(cosap) / sqrt(1.0 - cosap * cosap)), 2.0) * phiabs2;

                                        termx2 = ep * epsilon0 * astre * dcosapdphix * (-tan(al0)) / pow(1.0 - cosap * cosap, 1.5) * dcosapdx * phiabs2 + ep * epsilon0 * astre * ddcosapdphixdx * (1.0 - tan(al0) * cosap / sqrt(1.0 - cosap * cosap)) * phiabs2;
                                        termy2 = ep * epsilon0 * astre * dcosapdphiy * (-tan(al0)) / pow(1.0 - cosap * cosap, 1.5) * dcosapdy * phiabs2 + ep * epsilon0 * astre * ddcosapdphiydy * (1.0 - tan(al0) * cosap / sqrt(1.0 - cosap * cosap)) * phiabs2;
                                        termz2 = ep * epsilon0 * astre * dcosapdphiz * (-tan(al0)) / pow(1.0 - cosap * cosap, 1.5) * dcosapdz * phiabs2 + ep * epsilon0 * astre * ddcosapdphizdz * (1.0 - tan(al0) * cosap / sqrt(1.0 - cosap * cosap)) * phiabs2;

                                        termx3 = ep * epsilon0 * astre * dcosapdphix * (1.0 - tan(al0) * abs(cosap) / sqrt(1.0 - cosap * cosap)) * dphiabs2dx;
                                        termy3 = ep * epsilon0 * astre * dcosapdphiy * (1.0 - tan(al0) * abs(cosap) / sqrt(1.0 - cosap * cosap)) * dphiabs2dy;
                                        termz3 = ep * epsilon0 * astre * dcosapdphiz * (1.0 - tan(al0) * abs(cosap) / sqrt(1.0 - cosap * cosap)) * dphiabs2dz;

                                        termx = termx0 + termx1 + termx2 + termx3;
                                        termy = termy0 + termy1 + termy2 + termy3;
                                        termz = termz0 + termz1 + termz2 + termz3;

                                        termjjkk = termx + termy + termz;
                                    }
                                    else
                                    {
                                        termjjkk = aij[jj][kk] * (phidxx + phidyy + phidzz);
                                    }
                                }
                                else
                                {
                                    termjjkk = aij[jj][kk] * (phidxx + phidyy + phidzz);
                                }

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
