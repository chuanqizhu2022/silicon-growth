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
double ee;
double M0, W0, A0, S0, F0;
double Tm, dH;
double temp0, Tg, Tv, Tr;
int ni, nj;
int ix, iy, iz;
double r0, r, th0;
int xx0, yy0, zz0;
int intpos, frapass, curpos, curst, prepos, prest;
double sumplane;
double sump;
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
    astre = 1.0;
    ee = 0.0;
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

        double al, al0m, As, P, COS, SIN;
        double alii;

        double phidxii, phidyii, phidzii, phiabs2ii, nxii, nyii, nzii;

        double xxp, xyp, xzp, yxp, yyp, yzp, zxp, zyp, zzp;
        double I1, I2, I3;
        double cosap, ep;
        double dcosapdphix, dcosapdphiy, dcosapdphiz;
        double dphiabs2dx, dphiabs2dy, dphiabs2dz;
        double ddcosapdphixdx, ddcosapdphiydy, ddcosapdphizdz;
        double dcosapdx, dcosapdy, dcosapdz;

        double nxp, nyp, nzp;
        double nxpx, nxpy, nxpz;
        double nypx, nypy, nypz;
        double nzpx, nzpy, nzpz;
        double nx, ny, nz, nxx, nxy, nxz, nyx, nyy, nyz, nzx, nzy, nzz;

        double dPdx, dPdphix, ddPdphixdx;
        double dAsdx, dAsdphix, ddAsdphixdx;
        double nxphix, nyphix, nzphix;
        double nxphixdx, nyphixdx, nzphixdx;
        double ddPdphixdx_xy, ddPdphixdx_xz, ddPdphixdx_yz;

        double dPdy, dPdphiy, ddPdphiydy;
        double dAsdy, dAsdphiy, ddAsdphiydy;
        double nxphiy, nyphiy, nzphiy;
        double nxphiydy, nyphiydy, nzphiydy;
        double ddPdphiydy_xy, ddPdphiydy_xz, ddPdphiydy_yz;

        double dPdz, dPdphiz, ddPdphizdz;
        double dAsdz, dAsdphiz, ddAsdphizdz;
        double nxphiz, nyphiz, nzphiz;
        double nxphizdz, nyphizdz, nzphizdz;
        double ddPdphizdz_xy, ddPdphizdz_xz, ddPdphizdz_yz;

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
                        // sump = 0.0;
                        // for (kk = 0; kk <= nm; kk++)
                        // {
                        //     sump += phi[kk][i][j][k] * phi[kk][i][j][k];
                        // }
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

                        phidxii = (phi[ii][ip][j][k] - phi[ii][im][j][k]) / 2.0 / dx;
                        phidyii = (phi[ii][i][jp][k] - phi[ii][i][jm][k]) / 2.0 / dx;
                        phidzii = (phi[ii][i][j][kp] - phi[ii][i][j][km]) / 2.0 / dx;
                        phiabs2ii = phidxii * phidxii + phidyii * phidyii + phidzii * phidzii;
                        nxii = phidxii / sqrt(phiabs2ii);
                        nyii = phidyii / sqrt(phiabs2ii);
                        nzii = phidzii / sqrt(phiabs2ii);

                        alii = acos((abs(nxii) + abs(nyii) + abs(nzii)) / sqrt(3.0));

                        pddtt = 0.0;
                        for (n2 = 1; n2 <= phiNum[i][j][k]; n2++)
                        {
                            jj = phiIdx[n2][i][j][k];
                            sum1 = 0.0;
                            for (n3 = 1; n3 <= phiNum[i][j][k]; n3++)
                            {
                                kk = phiIdx[n3][i][j][k];

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

                                phidx = (phi[kk][ip][j][k] - phi[kk][im][j][k]) / 2.0 / dx;
                                phidy = (phi[kk][i][jp][k] - phi[kk][i][jm][k]) / 2.0 / dx;
                                phidz = (phi[kk][i][j][kp] - phi[kk][i][j][km]) / 2.0 / dx;

                                phidxp = phidx * xxp + phidy * yxp + phidz * zxp;
                                phidyp = phidx * xyp + phidy * yyp + phidz * zyp;
                                phidzp = phidx * xzp + phidy * yzp + phidz * zzp;

                                phidxx = (phi[kk][ip][j][k] + phi[kk][im][j][k] - 2.0 * phi[kk][i][j][k]) / dx / dx;
                                phidyy = (phi[kk][i][jp][k] + phi[kk][i][jm][k] - 2.0 * phi[kk][i][j][k]) / dx / dx;
                                phidzz = (phi[kk][i][j][kp] + phi[kk][i][j][km] - 2.0 * phi[kk][i][j][k]) / dx / dx;

                                phidxy = (phi[kk][ip][jp][k] + phi[kk][im][jm][k] - phi[kk][im][jp][k] - phi[kk][ip][jm][k]) / 4.0 / dx / dx;
                                phidxz = (phi[kk][ip][j][kp] + phi[kk][im][j][km] - phi[kk][im][j][kp] - phi[kk][ip][j][km]) / 4.0 / dx / dx;
                                phidyz = (phi[kk][i][jp][kp] + phi[kk][i][jm][km] - phi[kk][i][jm][kp] - phi[kk][i][jp][km]) / 4.0 / dx / dx;

                                phiabs2 = phidx * phidx + phidy * phidy + phidz * phidz;

                                dphiabs2dx = 2.0 * phidx * phidxx + 2.0 * phidy * phidxy + 2.0 * phidz * phidxz;
                                dphiabs2dy = 2.0 * phidx * phidxy + 2.0 * phidy * phidyy + 2.0 * phidz * phidyz;
                                dphiabs2dz = 2.0 * phidx * phidxz + 2.0 * phidy * phidyz + 2.0 * phidz * phidzz;

                                nx = phidx / sqrt(phiabs2);
                                ny = phidy / sqrt(phiabs2);
                                nz = phidz / sqrt(phiabs2);

                                nxp = phidxp / sqrt(phiabs2);
                                nyp = phidyp / sqrt(phiabs2);
                                nzp = phidzp / sqrt(phiabs2);

                                nxphix = -phidx * phidx / pow(phiabs2, 1.5) + 1.0 / sqrt(phiabs2);
                                nyphix = -phidx * phidy / pow(phiabs2, 1.5);
                                nzphix = -phidx * phidz / pow(phiabs2, 1.5);

                                // nxpphix = nxphix * xxp + nyphix * yxp + nzphix * zxp;
                                // nypphix = nxphix * xyp + nyphix * yyp + nzphix * zyp;
                                // nzpphix = nxphix * xzp + nyphix * yzp + nzphix * zzp;

                                nxphixdx = -1.0 / pow(phiabs2, 3.0) * (2.0 * phidx * phidxx * pow(phiabs2, 1.5) - 1.5 * sqrt(phiabs2) * dphiabs2dx * phidx * phidx) - 0.5 * pow(phiabs2, -1.5) * dphiabs2dx;
                                nyphixdx = -1.0 / pow(phiabs2, 3.0) * ((phidxx * phidy + phidx * phidxy) * pow(phiabs2, 1.5) - 1.5 * sqrt(phiabs2) * dphiabs2dx * phidx * phidy);
                                nzphixdx = -1.0 / pow(phiabs2, 3.0) * ((phidxx * phidz + phidx * phidxz) * pow(phiabs2, 1.5) - 1.5 * sqrt(phiabs2) * dphiabs2dx * phidx * phidz);

                                // nxpphixdx = nxphixdx * xxp + nyphixdx * yxp + nzphixdx * zxp;
                                // nypphixdx = nxphixdx * xyp + nyphixdx * yyp + nzphixdx * zyp;
                                // nzpphixdx = nxphixdx * xzp + nyphixdx * yzp + nzphixdx * zzp;

                                nxphiy = -phidx * phidy / pow(phiabs2, 1.5);
                                nyphiy = -phidy * phidy / pow(phiabs2, 1.5) + 1.0 / sqrt(phiabs2);
                                nzphiy = -phidz * phidy / pow(phiabs2, 1.5);

                                // nxpphiy = nxphiy * xxp + nyphiy * yxp + nzphiy * zxp;
                                // nypphiy = nxphiy * xyp + nyphiy * yyp + nzphiy * zyp;
                                // nzpphiy = nxphiy * xzp + nyphiy * yzp + nzphiy * zzp;

                                nxphiydy = -1.0 / pow(phiabs2, 3.0) * ((phidxy * phidy + phidx * phidyy) * pow(phiabs2, 1.5) - 1.5 * sqrt(phiabs2) * dphiabs2dy * phidx * phidy);
                                nyphiydy = -1.0 / pow(phiabs2, 3.0) * (2.0 * phidy * phidyy * pow(phiabs2, 1.5) - 1.5 * sqrt(phiabs2) * dphiabs2dy * phidy * phidy) - 0.5 * pow(phiabs2, -1.5) * dphiabs2dy;
                                nzphiydy = -1.0 / pow(phiabs2, 3.0) * ((phidyz * phidy + phidz * phidyy) * pow(phiabs2, 1.5) - 1.5 * sqrt(phiabs2) * dphiabs2dy * phidz * phidy);

                                // nxpphiydy = nxphiydy * xxp + nyphiydy * yxp + nzphiydy * zxp;
                                // nypphiydy = nxphiydy * xyp + nyphiydy * yyp + nzphiydy * zyp;
                                // nzpphiydy = nxphiydy * xzp + nyphiydy * yzp + nzphiydy * zzp;

                                nxphiz = -phidx * phidz / pow(phiabs2, 1.5);
                                nyphiz = -phidy * phidz / pow(phiabs2, 1.5);
                                nzphiz = -phidz * phidz / pow(phiabs2, 1.5) + 1.0 / sqrt(phiabs2);

                                // nxpphiz = nxphiz * xxp + nyphiz * yxp + nzphiz * zxp;
                                // nypphiz = nxphiz * xyp + nyphiz * yyp + nzphiz * zyp;
                                // nzpphiz = nxphiz * xzp + nyphiz * yzp + nzphiz * zzp;

                                nxphizdz = -1.0 / pow(phiabs2, 3.0) * ((phidxz * phidz + phidx * phidzz) * pow(phiabs2, 1.5) - 1.5 * sqrt(phiabs2) * dphiabs2dz * phidx * phidz);
                                nyphizdz = -1.0 / pow(phiabs2, 3.0) * ((phidyz * phidz + phidy * phidzz) * pow(phiabs2, 1.5) - 1.5 * sqrt(phiabs2) * dphiabs2dz * phidy * phidz);
                                nzphizdz = -1.0 / pow(phiabs2, 3.0) * (2.0 * phidz * phidzz * pow(phiabs2, 1.5) - 1.5 * sqrt(phiabs2) * dphiabs2dz * phidz * phidz) - 0.5 * pow(phiabs2, -1.5) * dphiabs2dz;

                                // nxpphizdz = nxphizdz * xxp + nyphizdz * yxp + nzphizdz * zxp;
                                // nypphizdz = nxphizdz * xyp + nyphizdz * yyp + nzphizdz * zyp;
                                // nzpphizdz = nxphizdz * xzp + nyphizdz * yzp + nzphizdz * zzp;

                                nxx = phidxx / sqrt(phiabs2) - phidx / (2.0 * pow(phiabs2, 1.5)) * dphiabs2dx;
                                nyx = phidxy / sqrt(phiabs2) - phidy / (2.0 * pow(phiabs2, 1.5)) * dphiabs2dx;
                                nzx = phidxz / sqrt(phiabs2) - phidz / (2.0 * pow(phiabs2, 1.5)) * dphiabs2dx;

                                nxy = phidxy / sqrt(phiabs2) - phidx / (2.0 * pow(phiabs2, 1.5)) * dphiabs2dy;
                                nyy = phidyy / sqrt(phiabs2) - phidy / (2.0 * pow(phiabs2, 1.5)) * dphiabs2dy;
                                nzy = phidyz / sqrt(phiabs2) - phidz / (2.0 * pow(phiabs2, 1.5)) * dphiabs2dy;

                                nxz = phidxz / sqrt(phiabs2) - phidx / (2.0 * pow(phiabs2, 1.5)) * dphiabs2dz;
                                nyz = phidyz / sqrt(phiabs2) - phidy / (2.0 * pow(phiabs2, 1.5)) * dphiabs2dz;
                                nzz = phidzz / sqrt(phiabs2) - phidz / (2.0 * pow(phiabs2, 1.5)) * dphiabs2dz;

                                nxpx = nxx * xxp + nyx * yxp + nzx * zxp;
                                nypx = nxx * xyp + nyx * yyp + nzx * zyp;
                                nzpx = nxx * xzp + nyx * yzp + nzx * zzp;

                                nxpy = nxy * xxp + nyy * yxp + nzy * zxp;
                                nypy = nxy * xyp + nyy * yyp + nzy * zyp;
                                nzpy = nxy * xzp + nyy * yzp + nzy * zzp;

                                nxpz = nxz * xxp + nyz * yxp + nzz * zxp;
                                nypz = nxz * xyp + nyz * yyp + nzz * zyp;
                                nzpz = nxz * xzp + nyz * yzp + nzz * zzp;

                                al = acos((abs(nxp) * I1 + abs(nyp) * I2 + abs(nzp) * I3) / sqrt(3.0));
                                al0m = asin(sqrt((ee * ee * (tan(al0) * tan(al0) - 1.0) + tan(al0) * tan(al0)) / (1.0 + tan(al0) * tan(al0))));
                                P = abs(nxp * nyp) * I1 * I2 + abs(nxp * nzp) * I1 * I3 + abs(nyp * nzp) * I2 * I3 + nxp * nxp * I1 * I1 / 2.0 + nyp * nyp * I2 * I2 / 2.0 + nzp * nzp * I3 * I3 / 2.0;
                                COS = sqrt((2.0 * P + 3.0 * ee * ee) / 3.0);
                                SIN = sqrt((3.0 - 2.0 * P + 3.0 * ee * ee) / 3.0);
                                As = 1.0 + astre * (COS + SIN * tan(al0));

                                if (anij[ii][kk] == 1)
                                {
                                    if (al < al0m && al > 1.0e-2)
                                    {
                                        epsilon0 = sqrt(aij[ii][kk]);

                                        dPdx = nxp * nyp / abs(nxp * nyp) * I1 * I2 * (nxpx * nyp + nxp * nypx) + nxp * nzp / abs(nxp * nzp) * I1 * I3 * (nxpx * nzp + nxp * nzpx) + nyp * nzp / abs(nyp * nzp) * I2 * I3 * (nypx * nzp + nyp * nzpx) + nxp * nxpx * I1 * I1 + nyp * nypx * I2 * I2 + nzp * nzpx * I3 * I3;
                                        dPdy = nxp * nyp / abs(nxp * nyp) * I1 * I2 * (nxpy * nyp + nxp * nypy) + nxp * nzp / abs(nxp * nzp) * I1 * I3 * (nxpy * nzp + nxp * nzpy) + nyp * nzp / abs(nyp * nzp) * I2 * I3 * (nypy * nzp + nyp * nzpy) + nxp * nxpy * I1 * I1 + nyp * nypy * I2 * I2 + nzp * nzpy * I3 * I3;
                                        dPdz = nxp * nyp / abs(nxp * nyp) * I1 * I2 * (nxpz * nyp + nxp * nypz) + nxp * nzp / abs(nxp * nzp) * I1 * I3 * (nxpz * nzp + nxp * nzpz) + nyp * nzp / abs(nyp * nzp) * I2 * I3 * (nypz * nzp + nyp * nzpz) + nxp * nxpz * I1 * I1 + nyp * nypz * I2 * I2 + nzp * nzpz * I3 * I3;

                                        dPdphix = 0.5 / abs(nx * ny) * I1 * I2 * (2.0 * nx * nxphix * ny * ny + nx * nx * 2.0 * ny * nyphix) + 0.5 / abs(nx * nz) * I1 * I3 * (2.0 * nx * nxphix * nz * nz + nx * nx * 2.0 * nz * nzphix) + 0.5 / abs(ny * nz) * I2 * I3 * (2.0 * ny * nyphix * nz * nz + ny * ny * 2.0 * nz * nzphix) + nx * nxphix * I1 * I1 + ny * nyphix * I2 * I2 + nz * nzphix * I3 * I3;
                                        dPdphiy = 0.5 / abs(nx * ny) * I1 * I2 * (2.0 * nx * nxphiy * ny * ny + nx * nx * 2.0 * ny * nyphiy) + 0.5 / abs(nx * nz) * I1 * I3 * (2.0 * nx * nxphiy * nz * nz + nx * nx * 2.0 * nz * nzphiy) + 0.5 / abs(ny * nz) * I2 * I3 * (2.0 * ny * nyphiy * nz * nz + ny * ny * 2.0 * nz * nzphiy) + nx * nxphiy * I1 * I1 + ny * nyphiy * I2 * I2 + nz * nzphiy * I3 * I3;
                                        dPdphiz = 0.5 / abs(nx * ny) * I1 * I2 * (2.0 * nx * nxphiz * ny * ny + nx * nx * 2.0 * ny * nyphiz) + 0.5 / abs(nx * nz) * I1 * I3 * (2.0 * nx * nxphiz * nz * nz + nx * nx * 2.0 * nz * nzphiz) + 0.5 / abs(ny * nz) * I2 * I3 * (2.0 * ny * nyphiz * nz * nz + ny * ny * 2.0 * nz * nzphiz) + nx * nxphiz * I1 * I1 + ny * nyphiz * I2 * I2 + nz * nzphiz * I3 * I3;

                                        ddPdphixdx_xy = I1 * I2 * (1.0 / pow(nx * ny, 2.0) * ((nxx * ny + nx * nyx) * abs(nx * ny) - 1.0 / abs(nx * ny) * (nx * nxx * ny * ny + nx * nx * ny * nyx) * nx * ny) * (ny * nxphix + nx * nyphix) + nx * ny / abs(nx * ny) * (nyx * nxphix + ny * nxphixdx + nxx * nyphix + nx * nyphixdx)) + I1 * I1 * (nxx * nxphix + nx * nxphixdx);
                                        ddPdphixdx_yz = I2 * I3 * (1.0 / pow(ny * nz, 2.0) * ((nyx * nz + ny * nzx) * abs(ny * nz) - 1.0 / abs(ny * nz) * (ny * nyx * nz * nz + ny * ny * nz * nzx) * ny * nz) * (nz * nyphix + ny * nzphix) + ny * nz / abs(ny * nz) * (nyx * nzphix + ny * nzphixdx + nyx * nzphix + ny * nzphixdx)) + I2 * I2 * (nyx * nyphix + ny * nyphixdx);
                                        ddPdphixdx_xz = I1 * I3 * (1.0 / pow(nx * nz, 2.0) * ((nxx * nz + nx * nzx) * abs(nx * nz) - 1.0 / abs(nx * nz) * (nx * nxx * nz * nz + nx * nx * nz * nzx) * nx * nz) * (nz * nxphix + nx * nzphix) + nx * nz / abs(nx * nz) * (nzx * nxphix + nz * nxphixdx + nxx * nzphix + nx * nzphixdx)) + I3 * I3 * (nzx * nzphix + nz * nzphixdx);

                                        ddPdphiydy_xy = I1 * I2 * (1.0 / pow(nx * ny, 2.0) * ((nxy * ny + nx * nyy) * abs(nx * ny) - 1.0 / abs(nx * ny) * (nx * nxy * ny * ny + nx * nx * ny * nyy) * nx * ny) * (ny * nxphiy + nx * nyphiy) + nx * ny / abs(nx * ny) * (nyy * nxphiy + ny * nxphiydy + nxy * nyphiy + nx * nyphiydy)) + I1 * I1 * (nxy * nxphiy + nx * nxphiydy);
                                        ddPdphiydy_yz = I2 * I3 * (1.0 / pow(ny * nz, 2.0) * ((nyy * nz + ny * nzy) * abs(ny * nz) - 1.0 / abs(ny * nz) * (ny * nyy * nz * nz + ny * ny * nz * nzy) * ny * nz) * (nz * nyphiy + ny * nzphiy) + ny * nz / abs(ny * nz) * (nyy * nzphiy + ny * nzphiydy + nyy * nzphiy + ny * nzphiydy)) + I2 * I2 * (nyy * nyphiy + ny * nyphiydy);
                                        ddPdphiydy_xz = I1 * I3 * (1.0 / pow(nx * nz, 2.0) * ((nxy * nz + nx * nzy) * abs(nx * nz) - 1.0 / abs(nx * nz) * (nx * nxy * nz * nz + nx * nx * nz * nzy) * nx * nz) * (nz * nxphiy + nx * nzphiy) + nx * nz / abs(nx * nz) * (nzy * nxphiy + nz * nxphiydy + nxy * nzphiy + nx * nzphiydy)) + I3 * I3 * (nzy * nzphiy + nz * nzphiydy);

                                        ddPdphizdz_xy = I1 * I2 * (1.0 / pow(nx * ny, 2.0) * ((nxz * ny + nx * nyz) * abs(nx * ny) - 1.0 / abs(nx * ny) * (nx * nxz * ny * ny + nx * nx * ny * nyz) * nx * ny) * (ny * nxphiz + nx * nyphiz) + nx * ny / abs(nx * ny) * (nyz * nxphiz + ny * nxphizdz + nxz * nyphiz + nx * nyphizdz)) + I1 * I1 * (nxz * nxphiz + nx * nxphizdz);
                                        ddPdphizdz_yz = I2 * I3 * (1.0 / pow(ny * nz, 2.0) * ((nyz * nz + ny * nzz) * abs(ny * nz) - 1.0 / abs(ny * nz) * (ny * nyz * nz * nz + ny * ny * nz * nzz) * ny * nz) * (nz * nyphiz + ny * nzphiz) + ny * nz / abs(ny * nz) * (nyz * nzphiz + ny * nzphizdz + nyz * nzphiz + ny * nzphizdz)) + I2 * I2 * (nyz * nyphiz + ny * nyphizdz);
                                        ddPdphizdz_xz = I1 * I3 * (1.0 / pow(nx * nz, 2.0) * ((nxz * nz + nx * nzz) * abs(nx * nz) - 1.0 / abs(nx * nz) * (nx * nxz * nz * nz + nx * nx * nz * nzz) * nx * nz) * (nz * nxphiz + nx * nzphiz) + nx * nz / abs(nx * nz) * (nzz * nxphiz + nz * nxphizdz + nxz * nzphiz + nx * nzphizdz)) + I3 * I3 * (nzz * nzphiz + nz * nzphizdz);

                                        ddPdphixdx = ddPdphixdx_xy + ddPdphixdx_xz + ddPdphixdx_yz;
                                        ddPdphiydy = ddPdphiydy_xy + ddPdphiydy_xz + ddPdphiydy_yz;
                                        ddPdphizdz = ddPdphizdz_xy + ddPdphizdz_xz + ddPdphizdz_yz;

                                        dAsdx = astre / 3.0 * (1.0 / COS - tan(al0) / SIN) * dPdx;
                                        dAsdy = astre / 3.0 * (1.0 / COS - tan(al0) / SIN) * dPdy;
                                        dAsdz = astre / 3.0 * (1.0 / COS - tan(al0) / SIN) * dPdz;

                                        dAsdphix = astre / 3.0 * (1.0 / COS - tan(al0m) / SIN) * dPdphix;
                                        dAsdphiy = astre / 3.0 * (1.0 / COS - tan(al0m) / SIN) * dPdphiy;
                                        dAsdphiz = astre / 3.0 * (1.0 / COS - tan(al0m) / SIN) * dPdphiz;

                                        ddAsdphixdx = astre * (-1.0 / 9.0) * (1.0 / pow(COS, 3.0) + tan(al0) / pow(SIN, 3.0)) * dPdx * dPdphix + astre * (1.0 / 3.0 / COS - tan(al0) / 3.0 / SIN) * ddPdphixdx;
                                        ddAsdphiydy = astre * (-1.0 / 9.0) * (1.0 / pow(COS, 3.0) + tan(al0) / pow(SIN, 3.0)) * dPdy * dPdphiy + astre * (1.0 / 3.0 / COS - tan(al0) / 3.0 / SIN) * ddPdphiydy;
                                        ddAsdphizdz = astre * (-1.0 / 9.0) * (1.0 / pow(COS, 3.0) + tan(al0) / pow(SIN, 3.0)) * dPdz * dPdphiz + astre * (1.0 / 3.0 / COS - tan(al0) / 3.0 / SIN) * ddPdphizdz;

                                        termx0 = 2.0 * epsilon0 * As * epsilon0 * dAsdx * phidx + phidxx * pow(epsilon0 * As, 2.0);
                                        termy0 = 2.0 * epsilon0 * As * epsilon0 * dAsdy * phidy + phidyy * pow(epsilon0 * As, 2.0);
                                        termz0 = 2.0 * epsilon0 * As * epsilon0 * dAsdz * phidz + phidzz * pow(epsilon0 * As, 2.0);

                                        termx1 = epsilon0 * dAsdx * epsilon0 * dAsdphix * phiabs2;
                                        termy1 = epsilon0 * dAsdy * epsilon0 * dAsdphiy * phiabs2;
                                        termz1 = epsilon0 * dAsdz * epsilon0 * dAsdphiz * phiabs2;

                                        termx2 = epsilon0 * As * epsilon0 * ddAsdphixdx * phiabs2;
                                        termy2 = epsilon0 * As * epsilon0 * ddAsdphiydy * phiabs2;
                                        termz2 = epsilon0 * As * epsilon0 * ddAsdphizdz * phiabs2;

                                        termx3 = epsilon0 * As * epsilon0 * dAsdphix * dphiabs2dx;
                                        termy3 = epsilon0 * As * epsilon0 * dAsdphiy * dphiabs2dy;
                                        termz3 = epsilon0 * As * epsilon0 * dAsdphiz * dphiabs2dz;

                                        termx = termx0 + termx1 + termx2 + termx3;
                                        termy = termy0 + termy1 + termy2 + termy3;
                                        termz = termz0 + termz1 + termz2 + termz3;

                                        termiikk = termx + termy + termz;
                                    }
                                    else if (al <= 1.0e-2)
                                    {
                                        termiikk = aij[ii][kk] * pow(1.0 + astre * (sqrt(cos(al) * cos(al) + ee * ee) + tan(al0) * sqrt(sin(al) * sin(al) + ee * ee)), 2.0) * (phidxx + phidyy + phidzz);
                                    }
                                    else
                                    {
                                        termiikk = aij[ii][kk] * pow(1.0 + astre * sqrt((1.0 + 2.0 * ee * ee) * (1.0 + tan(al0) * tan(al0))), 2.0) * (phidxx + phidyy + phidzz);
                                    }
                                }
                                else
                                {
                                    termiikk = aij[ii][kk] * (phidxx + phidyy + phidzz);
                                }

                                if (anij[jj][kk] == 1)
                                {
                                    if (al < al0m && al > 1.0e-2)
                                    {
                                        epsilon0 = sqrt(aij[jj][kk]);

                                        dPdx = nxp * nyp / abs(nxp * nyp) * I1 * I2 * (nxpx * nyp + nxp * nypx) + nxp * nzp / abs(nxp * nzp) * I1 * I3 * (nxpx * nzp + nxp * nzpx) + nyp * nzp / abs(nyp * nzp) * I2 * I3 * (nypx * nzp + nyp * nzpx) + nxp * nxpx * I1 * I1 + nyp * nypx * I2 * I2 + nzp * nzpx * I3 * I3;
                                        dPdy = nxp * nyp / abs(nxp * nyp) * I1 * I2 * (nxpy * nyp + nxp * nypy) + nxp * nzp / abs(nxp * nzp) * I1 * I3 * (nxpy * nzp + nxp * nzpy) + nyp * nzp / abs(nyp * nzp) * I2 * I3 * (nypy * nzp + nyp * nzpy) + nxp * nxpy * I1 * I1 + nyp * nypy * I2 * I2 + nzp * nzpy * I3 * I3;
                                        dPdz = nxp * nyp / abs(nxp * nyp) * I1 * I2 * (nxpz * nyp + nxp * nypz) + nxp * nzp / abs(nxp * nzp) * I1 * I3 * (nxpz * nzp + nxp * nzpz) + nyp * nzp / abs(nyp * nzp) * I2 * I3 * (nypz * nzp + nyp * nzpz) + nxp * nxpz * I1 * I1 + nyp * nypz * I2 * I2 + nzp * nzpz * I3 * I3;

                                        dPdphix = 0.5 / abs(nx * ny) * I1 * I2 * (2.0 * nx * nxphix * ny * ny + nx * nx * 2.0 * ny * nyphix) + 0.5 / abs(nx * nz) * I1 * I3 * (2.0 * nx * nxphix * nz * nz + nx * nx * 2.0 * nz * nzphix) + 0.5 / abs(ny * nz) * I2 * I3 * (2.0 * ny * nyphix * nz * nz + ny * ny * 2.0 * nz * nzphix) + nx * nxphix * I1 * I1 + ny * nyphix * I2 * I2 + nz * nzphix * I3 * I3;
                                        dPdphiy = 0.5 / abs(nx * ny) * I1 * I2 * (2.0 * nx * nxphiy * ny * ny + nx * nx * 2.0 * ny * nyphiy) + 0.5 / abs(nx * nz) * I1 * I3 * (2.0 * nx * nxphiy * nz * nz + nx * nx * 2.0 * nz * nzphiy) + 0.5 / abs(ny * nz) * I2 * I3 * (2.0 * ny * nyphiy * nz * nz + ny * ny * 2.0 * nz * nzphiy) + nx * nxphiy * I1 * I1 + ny * nyphiy * I2 * I2 + nz * nzphiy * I3 * I3;
                                        dPdphiz = 0.5 / abs(nx * ny) * I1 * I2 * (2.0 * nx * nxphiz * ny * ny + nx * nx * 2.0 * ny * nyphiz) + 0.5 / abs(nx * nz) * I1 * I3 * (2.0 * nx * nxphiz * nz * nz + nx * nx * 2.0 * nz * nzphiz) + 0.5 / abs(ny * nz) * I2 * I3 * (2.0 * ny * nyphiz * nz * nz + ny * ny * 2.0 * nz * nzphiz) + nx * nxphiz * I1 * I1 + ny * nyphiz * I2 * I2 + nz * nzphiz * I3 * I3;

                                        ddPdphixdx_xy = I1 * I2 * (1.0 / pow(nx * ny, 2.0) * ((nxx * ny + nx * nyx) * abs(nx * ny) - 1.0 / abs(nx * ny) * (nx * nxx * ny * ny + nx * nx * ny * nyx) * nx * ny) * (ny * nxphix + nx * nyphix) + nx * ny / abs(nx * ny) * (nyx * nxphix + ny * nxphixdx + nxx * nyphix + nx * nyphixdx)) + I1 * I1 * (nxx * nxphix + nx * nxphixdx);
                                        ddPdphixdx_yz = I2 * I3 * (1.0 / pow(ny * nz, 2.0) * ((nyx * nz + ny * nzx) * abs(ny * nz) - 1.0 / abs(ny * nz) * (ny * nyx * nz * nz + ny * ny * nz * nzx) * ny * nz) * (nz * nyphix + ny * nzphix) + ny * nz / abs(ny * nz) * (nyx * nzphix + ny * nzphixdx + nyx * nzphix + ny * nzphixdx)) + I2 * I2 * (nyx * nyphix + ny * nyphixdx);
                                        ddPdphixdx_xz = I1 * I3 * (1.0 / pow(nx * nz, 2.0) * ((nxx * nz + nx * nzx) * abs(nx * nz) - 1.0 / abs(nx * nz) * (nx * nxx * nz * nz + nx * nx * nz * nzx) * nx * nz) * (nz * nxphix + nx * nzphix) + nx * nz / abs(nx * nz) * (nzx * nxphix + nz * nxphixdx + nxx * nzphix + nx * nzphixdx)) + I3 * I3 * (nzx * nzphix + nz * nzphixdx);

                                        ddPdphiydy_xy = I1 * I2 * (1.0 / pow(nx * ny, 2.0) * ((nxy * ny + nx * nyy) * abs(nx * ny) - 1.0 / abs(nx * ny) * (nx * nxy * ny * ny + nx * nx * ny * nyy) * nx * ny) * (ny * nxphiy + nx * nyphiy) + nx * ny / abs(nx * ny) * (nyy * nxphiy + ny * nxphiydy + nxy * nyphiy + nx * nyphiydy)) + I1 * I1 * (nxy * nxphiy + nx * nxphiydy);
                                        ddPdphiydy_yz = I2 * I3 * (1.0 / pow(ny * nz, 2.0) * ((nyy * nz + ny * nzy) * abs(ny * nz) - 1.0 / abs(ny * nz) * (ny * nyy * nz * nz + ny * ny * nz * nzy) * ny * nz) * (nz * nyphiy + ny * nzphiy) + ny * nz / abs(ny * nz) * (nyy * nzphiy + ny * nzphiydy + nyy * nzphiy + ny * nzphiydy)) + I2 * I2 * (nyy * nyphiy + ny * nyphiydy);
                                        ddPdphiydy_xz = I1 * I3 * (1.0 / pow(nx * nz, 2.0) * ((nxy * nz + nx * nzy) * abs(nx * nz) - 1.0 / abs(nx * nz) * (nx * nxy * nz * nz + nx * nx * nz * nzy) * nx * nz) * (nz * nxphiy + nx * nzphiy) + nx * nz / abs(nx * nz) * (nzy * nxphiy + nz * nxphiydy + nxy * nzphiy + nx * nzphiydy)) + I3 * I3 * (nzy * nzphiy + nz * nzphiydy);

                                        ddPdphizdz_xy = I1 * I2 * (1.0 / pow(nx * ny, 2.0) * ((nxz * ny + nx * nyz) * abs(nx * ny) - 1.0 / abs(nx * ny) * (nx * nxz * ny * ny + nx * nx * ny * nyz) * nx * ny) * (ny * nxphiz + nx * nyphiz) + nx * ny / abs(nx * ny) * (nyz * nxphiz + ny * nxphizdz + nxz * nyphiz + nx * nyphizdz)) + I1 * I1 * (nxz * nxphiz + nx * nxphizdz);
                                        ddPdphizdz_yz = I2 * I3 * (1.0 / pow(ny * nz, 2.0) * ((nyz * nz + ny * nzz) * abs(ny * nz) - 1.0 / abs(ny * nz) * (ny * nyz * nz * nz + ny * ny * nz * nzz) * ny * nz) * (nz * nyphiz + ny * nzphiz) + ny * nz / abs(ny * nz) * (nyz * nzphiz + ny * nzphizdz + nyz * nzphiz + ny * nzphizdz)) + I2 * I2 * (nyz * nyphiz + ny * nyphizdz);
                                        ddPdphizdz_xz = I1 * I3 * (1.0 / pow(nx * nz, 2.0) * ((nxz * nz + nx * nzz) * abs(nx * nz) - 1.0 / abs(nx * nz) * (nx * nxz * nz * nz + nx * nx * nz * nzz) * nx * nz) * (nz * nxphiz + nx * nzphiz) + nx * nz / abs(nx * nz) * (nzz * nxphiz + nz * nxphizdz + nxz * nzphiz + nx * nzphizdz)) + I3 * I3 * (nzz * nzphiz + nz * nzphizdz);

                                        ddPdphixdx = ddPdphixdx_xy + ddPdphixdx_xz + ddPdphixdx_yz;
                                        ddPdphiydy = ddPdphiydy_xy + ddPdphiydy_xz + ddPdphiydy_yz;
                                        ddPdphizdz = ddPdphizdz_xy + ddPdphizdz_xz + ddPdphizdz_yz;

                                        dAsdx = astre / 3.0 * (1.0 / COS - tan(al0) / SIN) * dPdx;
                                        dAsdy = astre / 3.0 * (1.0 / COS - tan(al0) / SIN) * dPdy;
                                        dAsdz = astre / 3.0 * (1.0 / COS - tan(al0) / SIN) * dPdz;

                                        dAsdphix = astre / 3.0 * (1.0 / COS - tan(al0m) / SIN) * dPdphix;
                                        dAsdphiy = astre / 3.0 * (1.0 / COS - tan(al0m) / SIN) * dPdphiy;
                                        dAsdphiz = astre / 3.0 * (1.0 / COS - tan(al0m) / SIN) * dPdphiz;

                                        ddAsdphixdx = astre * (-1.0 / 9.0) * (1.0 / pow(COS, 3.0) + tan(al0) / pow(SIN, 3.0)) * dPdx * dPdphix + astre * (1.0 / 3.0 / COS - tan(al0) / 3.0 / SIN) * ddPdphixdx;
                                        ddAsdphiydy = astre * (-1.0 / 9.0) * (1.0 / pow(COS, 3.0) + tan(al0) / pow(SIN, 3.0)) * dPdy * dPdphiy + astre * (1.0 / 3.0 / COS - tan(al0) / 3.0 / SIN) * ddPdphiydy;
                                        ddAsdphizdz = astre * (-1.0 / 9.0) * (1.0 / pow(COS, 3.0) + tan(al0) / pow(SIN, 3.0)) * dPdz * dPdphiz + astre * (1.0 / 3.0 / COS - tan(al0) / 3.0 / SIN) * ddPdphizdz;

                                        termx0 = 2.0 * epsilon0 * As * epsilon0 * dAsdx * phidx + phidxx * pow(epsilon0 * As, 2.0);
                                        termy0 = 2.0 * epsilon0 * As * epsilon0 * dAsdy * phidy + phidyy * pow(epsilon0 * As, 2.0);
                                        termz0 = 2.0 * epsilon0 * As * epsilon0 * dAsdz * phidz + phidzz * pow(epsilon0 * As, 2.0);

                                        termx1 = epsilon0 * dAsdx * epsilon0 * dAsdphix * phiabs2;
                                        termy1 = epsilon0 * dAsdy * epsilon0 * dAsdphiy * phiabs2;
                                        termz1 = epsilon0 * dAsdz * epsilon0 * dAsdphiz * phiabs2;

                                        termx2 = epsilon0 * As * epsilon0 * ddAsdphixdx * phiabs2;
                                        termy2 = epsilon0 * As * epsilon0 * ddAsdphiydy * phiabs2;
                                        termz2 = epsilon0 * As * epsilon0 * ddAsdphizdz * phiabs2;

                                        termx3 = epsilon0 * As * epsilon0 * dAsdphix * dphiabs2dx;
                                        termy3 = epsilon0 * As * epsilon0 * dAsdphiy * dphiabs2dy;
                                        termz3 = epsilon0 * As * epsilon0 * dAsdphiz * dphiabs2dz;

                                        termx = termx0 + termx1 + termx2 + termx3;
                                        termy = termy0 + termy1 + termy2 + termy3;
                                        termz = termz0 + termz1 + termz2 + termz3;

                                        termjjkk = termx + termy + termz;
                                    }
                                    else if (al <= 1.0e-2)
                                    {
                                        termjjkk = aij[jj][kk] * pow(1.0 + astre * (sqrt(cos(al) * cos(al) + ee * ee) + tan(al0) * sqrt(sin(al) * sin(al) + ee * ee)), 2.0) * (phidxx + phidyy + phidzz);
                                    }
                                    else
                                    {
                                        termjjkk = aij[jj][kk] * pow(1.0 + astre * sqrt((1.0 + 2.0 * ee * ee) * (1.0 + tan(al0) * tan(al0))), 2.0) * (phidxx + phidyy + phidzz);
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
                            if (alii < (2.0 / 180.0 * PI))
                            {
                                miijj = mij[ii][jj] / 3.0;
                            }
                            else if ((alii >= (2.0 / 180.0 * PI)) && (alii < (4.0 / 180.0 * PI)))
                            {
                                miijj = (sin(90.0 * (alii - 3.0 / 180.0 * PI)) + 2.0) / 3.0 * mij[ii][jj];
                            }
                            else
                            {
                                miijj = mij[ii][jj];
                            }
                            pddtt += -2.0 * miijj / (double)phiNum[i][j][k] * (sum1 - 8.0 / PI * F0 * sqrt(phi[ii][i][j][k] * phi[jj][i][j][k]));
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
