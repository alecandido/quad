#include <math.h>
#include <float.h>
#include <stdio.h>

int hello() {
    return 42;
}

void modify_double( double *a){
    *a = 10.0;
    printf("%f", a);
}

void looop(){
    int i,j ;
    for (j = 0; j < 100 ; j++){
        i += 1;
        printf("%d",i);
    }
}

// ----------------------------------------------------------------

#define GSL_DBL_MIN        2.2250738585072014e-308
#define GSL_DBL_EPSILON        2.2204460492503131e-16

double rescale_error (double err, const double result_abs, const double result_asc)
{
  err = fabs(err) ;

  if (result_asc != 0 && err != 0)
      {
        double scale = pow((200 * err / result_asc), 1.5) ;

        if (scale < 1)
          {
            err = result_asc * scale ;
          }
        else
          {
            err = result_asc ;
          }
      }
  if (result_abs > GSL_DBL_MIN / (50 * GSL_DBL_EPSILON))
    {
      double min_err = 50 * GSL_DBL_EPSILON * result_abs ;

      if (min_err > err)
        {
          err = min_err ;
        }
    }
  return err ;
}

// --------------------------------------------------------------

static const double xgk[31] =   /* abscissae of the 61-point kronrod rule */
{
  0.999484410050490637571325895705811,
  0.996893484074649540271630050918695,
  0.991630996870404594858628366109486,
  0.983668123279747209970032581605663,
  0.973116322501126268374693868423707,
  0.960021864968307512216871025581798,
  0.944374444748559979415831324037439,
  0.926200047429274325879324277080474,
  0.905573307699907798546522558925958,
  0.882560535792052681543116462530226,
  0.857205233546061098958658510658944,
  0.829565762382768397442898119732502,
  0.799727835821839083013668942322683,
  0.767777432104826194917977340974503,
  0.733790062453226804726171131369528,
  0.697850494793315796932292388026640,
  0.660061064126626961370053668149271,
  0.620526182989242861140477556431189,
  0.579345235826361691756024932172540,
  0.536624148142019899264169793311073,
  0.492480467861778574993693061207709,
  0.447033769538089176780609900322854,
  0.400401254830394392535476211542661,
  0.352704725530878113471037207089374,
  0.304073202273625077372677107199257,
  0.254636926167889846439805129817805,
  0.204525116682309891438957671002025,
  0.153869913608583546963794672743256,
  0.102806937966737030147096751318001,
  0.051471842555317695833025213166723,
  0.000000000000000000000000000000000
};

static const double wgk[31] =   /* weights of the 61-point kronrod rule */
{
  0.001389013698677007624551591226760,
  0.003890461127099884051267201844516,
  0.006630703915931292173319826369750,
  0.009273279659517763428441146892024,
  0.011823015253496341742232898853251,
  0.014369729507045804812451432443580,
  0.016920889189053272627572289420322,
  0.019414141193942381173408951050128,
  0.021828035821609192297167485738339,
  0.024191162078080601365686370725232,
  0.026509954882333101610601709335075,
  0.028754048765041292843978785354334,
  0.030907257562387762472884252943092,
  0.032981447057483726031814191016854,
  0.034979338028060024137499670731468,
  0.036882364651821229223911065617136,
  0.038678945624727592950348651532281,
  0.040374538951535959111995279752468,
  0.041969810215164246147147541285970,
  0.043452539701356069316831728117073,
  0.044814800133162663192355551616723,
  0.046059238271006988116271735559374,
  0.047185546569299153945261478181099,
  0.048185861757087129140779492298305,
  0.049055434555029778887528165367238,
  0.049795683427074206357811569379942,
  0.050405921402782346840893085653585,
  0.050881795898749606492297473049805,
  0.051221547849258772170656282604944,
  0.051426128537459025933862879215781,
  0.051494729429451567558340433647099
};

static const double wg[15] =    /* weights of the 30-point gauss rule */
{
  0.007968192496166605615465883474674,
  0.018466468311090959142302131912047,
  0.028784707883323369349719179611292,
  0.038799192569627049596801936446348,
  0.048402672830594052902938140422808,
  0.057493156217619066481721689402056,
  0.065974229882180495128128515115962,
  0.073755974737705206268243850022191,
  0.080755895229420215354694938460530,
  0.086899787201082979802387530715126,
  0.092122522237786128717632707087619,
  0.096368737174644259639468626351810,
  0.099593420586795267062780282103569,
  0.101762389748405504596428952168554,
  0.102852652893558840341285636705415
};

void
gsl_integration_qk (const int n,double a, double b,double *result, double *abserr,double *resabs, double *resasc)
{
    double fv1[31], fv2[31];


   const double center = 0.5 * (a + b);
   const double half_length = 0.5 * (b - a);
   const double abs_half_length = fabs(half_length);
   const double f_center = cos(center);
   double result_gauss = 0;
   double result_kronrod = f_center * wgk[n - 1];

   double result_abs = fabs(result_kronrod);
   double result_asc = 0;
   double mean = 0, err = 0;

   int j;
   if (n % 2 == 0)
       {
         result_gauss = f_center * wg[n / 2 - 1];
       }

   for (j = 0; j < (n - 1) / 2; j++)
     {
       const int jtw = j * 2 + 1;  /* in original fortran j=1,2,3 jtw=2,4,6 */
       const double abscissa = half_length * xgk[jtw];
       const double fval1 = cos(center - abscissa);
       const double fval2 = cos(center + abscissa);
       const double fsum = fval1 + fval2;
       fv1[jtw] = fval1;
       fv2[jtw] = fval2;
       result_gauss += wg[j] * fsum;
       result_kronrod += wgk[jtw] * fsum;
       result_abs += wgk[jtw] * (fabs (fval1) + fabs (fval2));
     }

   for (j = 0; j < n / 2; j++)
     {
       int jtwm1 = j * 2;
       const double abscissa = half_length * xgk[jtwm1];
       const double fval1 = cos(center - abscissa);
       const double fval2 = cos(center + abscissa);
       fv1[jtwm1] = fval1;
       fv2[jtwm1] = fval2;
       result_kronrod += wgk[jtwm1] * (fval1 + fval2);
       result_abs += wgk[jtwm1] * (fabs (fval1) + fabs (fval2));
     };
   mean = result_kronrod * 0.5;

 result_asc = wgk[n - 1] * fabs (f_center - mean);

 for (j = 0; j < n - 1; j++)
   {
     result_asc += wgk[j] * (fabs (fv1[j] - mean) + fabs (fv2[j] - mean));
   }
 /* scale by the width of the integration region */

 err = (result_kronrod - result_gauss) * half_length;

 result_kronrod *= half_length;
 result_abs *= abs_half_length;
 result_asc *= abs_half_length;

 *result = result_kronrod;
 *resabs = result_abs;
 *resasc = result_asc;
 *abserr = rescale_error (err, result_abs, result_asc);

}
