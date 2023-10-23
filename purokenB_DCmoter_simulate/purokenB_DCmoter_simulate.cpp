/********************************************************************
Title	Project_B_No1

【更新履歴】
2021.10.13
2021.10.15 ルンゲクッタ法記述
2022.10.31 ASR記述, 2022パワエレ第6回の課題に合わせて微修正
*********************************************************************/

#pragma warning(disable : 4996) // エラー消し

#define _CRT_SECURE_NO_WARNINGS

//#include <iostream>
//#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

// /*キャリア周波数*/
// #define FS 5000.0

/*数学*/
#define pi 3.141529265 /*円周率*/

/***時間系***/
#define zt_end 20        /*シミュレーション終了時刻*/
#define zh 0.0000001     /*シミュレーション周期 1[us]*/
#define zh_CNTL 0.000100 /*制御周期 100[μs]*/
// #define zh_CNTL (1/Fs)
#define zh_disp 1.00    /*結果表示記録周期 1.0[s]*/
#define zh_rec 0.001000 /*結果記録周期 1[ms]*/

// モデル諸元//
#define R 1.0   /*巻線抵抗*/
#define L 0.04    /*自己インダクタンス*/
#define KPhif 1.0 /*起電力係数*/
#define Jm 0.01   /*慣性負荷*/
#define Tl 0.10   /*負荷トルク*/

/*データ記録用*/
#define data1 "zt, wre_ref, wre, wre_refF, wm_ref, wm, wm_refF, va, ia, dwre, zdwre, e_induced\n" /*変数名*/
#define data2 "%lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf\n"                      /*枠*/
#define data3 zt, wre_ref, wre, wre_refF, wm_ref, wm, wm_refF, va, ia, dwre, zdwre, e_induced     /*変数*/

// 各定数定義//
double ia = 0.0;        /*電機子電流*/
double wre = 0.0;       /*角速度 [rad/s]*/
double wm = 0.0;        /*角速度 [rpm]*/
double va = 0.0;        /*電圧*/
double e_induced = 0.0; // 誘起電圧

/*ASR*/
double KpASR = 0.0;
double KiASR = 0.0;
double wre_ref = 0.0;  /*角速度指令値 [rad/s]*/
double wre_refF = 0.0; /*角速度指令値の1次遅れ [rad/s]*/
double wm_ref = 0.0;   /*角速度指令値 [rpm]*/
double wm_refF = 0.0;  /*角速度指令値の1次遅れ [rpm]*/
double dwre = 0.0;     /*角速度の偏差*/
double zdwre = 0.0;    /*角速度の偏差積分値*/
double ia_ref = 0.0;   /*ASRから出力される電機子電流の指令値*/
double wc_ASR = 0.0;   /*ASRのカットオフ周波数[rad/s]*/

double wc = 0.0; /*速度応答のカットオフ周波数[rad/s]*/
double Td = 0.0; /*速度応答の目標時定数[s]*/

/*ACR*/
double Td_ACR = 0.0; /*ACRの時定数*/
double KpACR = 0.0;  // 5.0;
double KiACR = 0.0;  // 125.0;
double dia = 0.0;    /*電流の偏差*/
double zdia = 0.0;   /*電流の偏差積分値*/
double va_PI = 0.0;  /*ACRから出力される電機子電圧*/
double va_FF = 0.0;  /*誘起電圧分の電機子電圧を補償するフィードフォワード項*/

/* 時間変数 */
double zt = 0.0;      /*時刻*/
double zt_CNTL = 0.0; /*制御時刻*/
double zt_disp = 0.0; /*結果表示用時刻*/
double zt_rec = 0.0;  /*結果記録用時刻*/

int result; /* 結果を記録するためのフラグ */

// ルンゲクッタ用
double k1ia, k2ia, k3ia, k4ia;
double k1wre, k2wre, k3wre, k4wre;

/***1次遅れフィルター****/
double funcdely(double zhh, double Tm, double uF_old, double u)
{
    double aF_new;
    aF_new = (u * zhh + uF_old * Tm) / (Tm + zhh);
    return (aF_new);
}
/*********************/

// 電流の微分方程式
double funcdiffia(double aia, double va, double awre)
{
    double diffaia;
    diffaia = (va - R * aia - KPhif * awre) / L; /*電機子の回路方程式*/
    return (diffaia);
}
/*****************/

// 速度の微分方程式
double funcdiffwre(double aia)
{
    double diffawre;
    diffawre = (KPhif * aia - Tl) / Jm; /*回転計の方程式*/
    return (diffawre);
}
/*****************/

// rpm -> rad/s 単位変換
double rpm_rads(double rpm)
{
    double rads;
    rads = rpm * 2 * pi / 60;
    return (rads);
}
/****************/

// rad/s -> rpm 単位変換
double rads_rpm(double rads)
{
    double rpm;
    rpm = rads * 60 / 2 / pi;
    return (rpm);
}
/****************/

int main()
{
    /**************************************************************/
    wc = 25.0;      /*カットオフ周波数25[rad/s]*/
    Td = 1.0 / wc;  /*時定数*/
    Td_ACR = 0.008; /*電流制御系時定数8[ms]*/
    wc_ASR = 2.5;   /*速度制御系カットオフ周波数5[rad/s]*/
    /***************************************************************/

    KpACR = (L / Td_ACR) * 1.0;
    KiACR = (R / Td_ACR) * 1.0;

    KpASR = Jm * wc / KPhif;
    KiASR = KpASR * wc_ASR;

    FILE* fp;
    time_t start, end;
    struct tm* ltstart;
    char filename[35];

    /*ファイル作成*/
    time(&start);
    ltstart = localtime(&start);

    filename[0] = 'P';
    filename[1] = 'o';
    filename[2] = 'w';
    filename[3] = 'e';
    filename[4] = 'r';
    filename[5] = 'E';
    filename[6] = 'l';
    filename[7] = 'e';
    filename[8] = '-';
    filename[9] = 'N';
    filename[10] = 'O';
    filename[11] = '6';
    filename[12] = '_';
    filename[13] = '0' + (ltstart->tm_year + 1900) / 1000;
    filename[14] = '0' + ((ltstart->tm_year + 1900) % 1000) / 100;
    filename[15] = '0' + ((ltstart->tm_year + 1900) % 100) / 10;
    filename[16] = '0' + ((ltstart->tm_year + 1900) % 10);
    filename[17] = ',';
    filename[18] = '0' + (ltstart->tm_mon + 1) / 10;
    filename[19] = '0' + (ltstart->tm_mon + 1) % 10;
    filename[20] = '_';
    filename[21] = '0' + ltstart->tm_mday / 10;
    filename[22] = '0' + ltstart->tm_mday % 10;
    filename[23] = ',';
    filename[24] = '0' + ltstart->tm_hour / 10;
    filename[25] = '0' + ltstart->tm_hour % 10;
    filename[26] = '_';
    filename[27] = '0' + ltstart->tm_min / 10;
    filename[28] = '0' + ltstart->tm_min % 10;
    filename[29] = '.';
    filename[30] = 'c';
    filename[31] = 's';
    filename[32] = 'v';
    filename[33] = '\0';

    fp = fopen(filename, "w");

    printf("シミュレーション開始 \n");
    printf("\n");

    /*** システムパラメータ設定 ***/
    // zt_end = 20.0; /* シミュレーション終了時間 */
    /******************************/
    printf("KpACR=%lf, KiACR=%lf, KpASR=%lf, KiASR=%lf\n", KpACR, KiACR, KpASR, KiASR);
    printf("zt wm_ref wm\n");
    fprintf(fp, data1);

    result = 1; /** 結果を記録するためのフラグ **/

    /****************************** MAIN LOOP ******************************/
    for (zt = 0.0; zt < zt_end; zt += zh)
    {

        if (zt >= zt_rec)
        {
            /* 結果を表示 */
            if (result == 1)
            {
                fprintf(fp, data2, data3);
            }

            zt_rec = zt_rec + zh_rec; /* 次のRESULT後まで結果は記録しない */
        }

        if (zt >= zt_disp)
        {
            printf("%lf %lf %lf\n", zt, wm_ref, wm);

            zt_disp = zt_disp + zh_disp; /* 次のRESULT1後まで結果は表示しない */
        }

        // 制御ループ
        if (zt > zt_CNTL)
        {
            /*運転方法*/
            if (zt < 1)
            {
                wm_ref = 0.0;
                wre_ref = rpm_rads(wm_ref);
            }
            if (zt >= 1.0 && zt < 8.0) /**定常状態となるよう十分な時間とる**/
            {
                wm_ref = 600.0;
                wre_ref = rpm_rads(wm_ref);
            }
            if (zt >= 8.0)
            {
                wm_ref = 800.0;
                wre_ref = rpm_rads(wm_ref);
            }

            /**ASR**/
            dwre = wre_ref - wre;             /*偏差*/
            zdwre = zdwre + (dwre * zh_CNTL); /*積分値*/
            ia_ref = KpASR * dwre + KiASR * zdwre;

            /**ACR**/
            dia = ia_ref - ia;
            zdia = zdia + (dia * zh_CNTL);
            va_PI = KpACR * dia + KiACR * zdia;
            va_FF = KPhif * wre;

            va = va_PI + va_FF;
            if (va > 100)
            {
                va = 100;
            }

            zt_CNTL = zt_CNTL + zh_CNTL;
        }

        /********** ルンゲクッタ **********/
        k1ia = zh * funcdiffia(ia, va, wre);
        k1wre = zh * funcdiffwre(ia);

        k2ia = zh * funcdiffia(ia + k1ia / 2.0, va, wre + k1wre / 2.0);
        k2wre = zh * funcdiffwre(ia + k1ia / 2.0);

        k3ia = zh * funcdiffia(ia + k2ia / 2.0, va, wre + k2wre / 2.0);
        k3wre = zh * funcdiffwre(ia + k2ia / 2.0);

        k4ia = zh * funcdiffia(ia + k1ia, va, wre + k1wre);
        k4wre = zh * funcdiffwre(ia + k1ia);

        ia = ia + (k1ia * 2.0 + k2ia * 2.0 + k3ia + k4ia) / 6.0;
        wre = wre + (k1wre * 2.0 + k2wre * 2.0 + k3wre + k4wre) / 6.0;
        /********************************/

        wm = rads_rpm(wre);

        wre_refF = funcdely(zh, Td, wre_refF, wre_ref);
        wm_refF = rads_rpm(wre_refF);

        /*誘起電圧*/
        e_induced = wre * KPhif;
    }
    fclose(fp);
    time(&end);
    putchar('\n');
    printf("結果を %s に出力しました\n", filename);
    printf("シミュレーションに%d秒かかりました\n", end - start);
    // getchar();

    return 0;
}