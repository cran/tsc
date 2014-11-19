#include "model_t_s1.h"
#include <iostream>
#include <math.h>
#include <numeric>
#include <vector>
#include <algorithm>  //std::nth_element
#include "R.h"
#include "Rmath.h"
#define pi 3.1415926
#define delta 0.1


//using namespace std;

//function to get the proposed test-statistics
double function(vector<double> y1, vector<double> y2)
{
   int n1 = y1.size();  //sample size of y1
   int n2 = y2.size();  //sample size of y2
   vector<double> combine=y1;  //combine vector y1 and vector y2
   combine.insert(combine.end(),y2.begin(),y2.end());
   double sum1=std::accumulate(combine.begin(),combine.end(),0.0);
   double mean1=sum1/combine.size();
   double sq_sum=std::inner_product(combine.begin(),combine.end(),combine.begin(),0.0);
   double variance=sq_sum / combine.size() - mean1 * mean1;
    //get the value of m1(1,2,...minimum of round(n1^(1-delta)) and round(n1/2)
    //get the minimum value of round(n1^(1-delta)) and round(n1/2)
   int m1_min=min(round(pow(n1,1-delta)),round(n1/2));
    //now get the vector for m1.
    //get the value of m2(1,2,...minimum of round(n2^(1-delta)) and round(n2/2)
    //get the minimum value of round(n2^(1-delta)) and round(n2/2)
   int m2_min=min(round(pow(n2,1-delta)),round(n2/2));
    //now get the vector for m1.
   double partc_1 =  10000000;
   sort(y1.begin(),y1.end());   //sort y1

   for (int m1=1; m1<=m1_min; m1++) { //get the part C
   // double temp=function1(y1,i,n1);
   int lower1,upper1;
   double test1 = 0;

   for (int j1=1; j1<=n1;j1++)  {
   if (j1<1+m1) lower1=1; //if j-m <1, then j-m should be equal to 1
          else if (j1>=1+m1)  lower1=j1-m1; //if j-m >=1, then j-m should be j-m
    //get the j+m value
   if (j1+m1>n1) upper1=n1;  //if j+m >n(sample size), then j+m should be equal to n(sample size)
          else upper1=j1+m1; //if j+m <=n(sample size), then j+m should be equal to j+m(sample size)
    //y so we can get the (j-m)th order statistics and (j+m)th order statistics

   double y_lower1=y1[lower1-1];  //find the (j-m)th order statistics

   double y_upper1=y1[upper1-1];  //find the (j+m)th order statistics
   double y_dif1=y_upper1-y_lower1;

   if (y_dif1==0) y_dif1=(n1+n2)/2;
   //j is from 1 to n,the propose is to get the production of m/[n*(y_j+m-y_j-m)]
        //get the j-m value
         //the differenc
   double test_1=log(2*((int)m1)/((double)n1*y_dif1)); //each j(j=1...n), we have each m/[j*(y_dif)]
   test1=test1+test_1; //sum all the test1 (we have n test_1)
   }
   if(test1<partc_1)
       partc_1 =test1;

    //now we can have m1 partc, then we get the min of partc
   }

   double partd_1 = 10000000;
   sort(y2.begin(),y2.end());   //sort y2
   for ( int m2=1; m2<=m2_min; m2++) { //get the part d
   // Rprintf("m2=%d\n",m2);
   // double temp=function1(y2,i,n2);
   int lower2,upper2;
   double test2 = 0;
    // Rprintf("n=%d\n",n);
   for ( int j2=1; j2<=n2;j2++)  {
       // Rprintf("j2=%d\n",j2);
    if (j2<1+m2) lower2=1; //if j-m <1, then j-m should be equal to 1
          else  lower2=j2-m2; //if j-m >=1, then j-m should be j-m
         // Rprintf("lower2=%d\n",lower2);
        //get the j+m value
        if (j2+m2>n2) upper2=n2;  //if j+m >n(sample size), then j+m should be equal to n(sample size)
          else upper2=j2+m2; //if j+m <=n(sample size), then j+m should be equal to j+m(sample size)
      //  Rprintf("upper=%lf\n",upper);
   //sort y so we can get the (j-m)th order statistics and (j+m)th order statistics

   double y_lower2=y2[lower2-1];  //find the (j-m)th order statistics

   double y_upper2=y2[upper2-1];  //find the (j+m)th order statistics
   //Rprintf("y_lower=%lf,y_upper=%lf\n",y_lower,y_upper);

   double y_dif2=y_upper2-y_lower2;
   if (y_dif2==0) y_dif2=(n1+n2)/2;
   //Rprintf("y_dif2=%lf\n",y_dif2);
    //j is from 1 to n,the propose is to get the production of m/[n*(y_j+m-y_j-m)]
        //get the j-m value
         //the differenc
   double test_2=log(2*((int)m2)/((double)n2*y_dif2)); //each j(j=1...n), we have each m/[j*(y_dif)]
  // Rprintf("test_1=%lf\n",test_1);
   test2=test2+test_2; //multiply all the test1 (we have n test_1)
   }

    if(test2<partd_1)
        partd_1 = test2;
    //now we can have m2 partd, then we get the min of partd
    }
      //Rprintf("partd_1=%lf",partd_1);

    //now get the test statistics

    double test_stat1=((double)(n1+n2)/2)*log(2*pi*exp(1)*variance)+partc_1+partd_1;
    // double test_stat1=((double)(n1+n2)/2)*log(2*pi*exp(1)*variance)+partc_1+test2;
    //double test_stat=pow(2*pi*exp(1)*variance,((double)(n1)/2))*pow(2*pi*exp(1)*variance,((double)(n2)/2))*partc_1*partd_1;
    //Rprintf("a=%lf,b=%lf\n",2*pi*exp(1)*variance,((double)(n1)/2));
    //Rprintf("A=%lf, B=%lf, C=%lf, D=%lf",pow(2*pi*exp(1)*variance,((double)(n1)/2)),pow(2*pi*exp(1)*variance,((double)(n2)/2)), partc_1, partd_1 );
    //Rprintf("pow1=%lf\n",test_stat);
    //double test_stat1=log(test_stat);
   // Rprintf("pow2=%lf\n",test_stat1);
    return test_stat1;
    //Rprintf("test=%lf\n",test_stat1);

}
//int *ndouble1, int *ndouble2, double* y1, double* y2, double* mc , double* test_statistics
void test (int *ndouble1, int *ndouble2, double* y1, double* y2, double* test_statistics){
    vector<double> v_y1;
    vector<double> v_y2;


   for(int i = 0; i < *ndouble1; i++){
        v_y1.push_back(y1[i]);
    }
    for(int i = 0; i < *ndouble2; i++){
        v_y2.push_back(y2[i]);
    }
    *test_statistics = function(v_y1,v_y2);

}


double foo(vector<double> y1, vector<double> y2,int num)
{
  double p, p1;
  int count;
  count=0;
  //lower1=0;
  //upper1=0;
  double t=0;
  double log_sum=0;
  double t_2_sum=0;
  //double log_t_sum=0;
  double test=function(y1,y2);  //get the test statistics of the data which the user provides to us
  //Monte Carlo simulation to get the 10000 test statistics from standard normal distribution
  vector<double>test11(num);


  for (int i=1; i <=num; i++)  {
        vector<double>x1(y1.size());
        vector<double>x2(y2.size());
        for(unsigned j = 0; j<y1.size(); j++ ){
        x1[j]=rnorm(0,1)
        ;
        }
        for(unsigned j =0; j < y2.size(); j++){
        x2[j]=rnorm(0,1);
        }
   double test1=function(x1,x2);

   test11[i-1]=function(x1,x2);

   double t_2=test1*test1;
   //Rprintf("t_2=%lf\n",t_2);
   t_2_sum=t_2_sum+t_2;//sum(t^2)
   //Rprintf("t_2_sum=%lf\n",t_2_sum);
   if (test1>test) count=count+1; //get the the number of how many test1(10000) is less than test(from the user's data)
        else  count=count;
   t=t+test1;//sum(t)
   //Rprintf("t=%lf\n",t);
//   if (num1>num+100) num1=1;
  //    else num1=num1;
  // Rprintf("num1=%lf\n",num1);
  // double log1=log(num1);
 //   Rprintf("log1=%lf\n",log1);
  //  Rprintf("test11=%lf\n",test11[i-1]);
 //  double log_t=log1*test1;
//    Rprintf("log_t=%lf\n",log_t);
 //   log_t_sum=log_t_sum+log_t; //sum(y*t)
   // log_sum=log_sum+log1; //sum(y)
    // Rprintf("log_t_sum=%lf\n",log_t_sum);
    //Rprintf("log_sum=%lf\n",log_sum);



  //Rprintf("log_t_sum=%lf\n",log_t_sum);
//Rprintf("log1=%lf\n",log1);

// Rprintf("count=%lf\n",count);
  //Rprintf("mc=%lf\n",num);

  }
sort(test11.begin(),test11.end());
 double log_t_sum=0;
 for (int i=1; i <=num; i++)  {
 double num1=i/(num-i+0.0000001);
 if (num1>num+100) num1=1;
      else num1=num1;
    double log1=log(num1);
 //   Rprintf("log1=%lf\n",log1);
  //  Rprintf("test11=%lf\n",test11[i-1]);
   double log_t=log1*test11[i-1];
//    Rprintf("log_t=%lf\n",log_t);
    log_t_sum=log_t_sum+log_t; //sum(y*t)
    log_sum=log_sum+log1; //sum(y)
 }
double c1=(log_t_sum-(log_sum/(num-1))*t)/((t_2_sum)-(t/(num-1))*t);
double c0=log_sum/(num-1)-c1*(t/(num-1));
 // Rprintf("test1=%lf\n",test1);

 //   Rprintf("log_t_sum=%lf\n",log_t_sum);
 //Rprintf("log1=%lf\n",log1);

 // Rprintf("count=%lf\n",count);
 //Rprintf("mc=%lf\n",num);
 // Rprintf("t=%lf\n",t);
  //Rprintf("count=%lf\n",count);
  //Rprintf("lower=%lf\n",lower);
  //Rprintf("upper=%lf\n",upper);
 //  Rprintf("logtest1=%lf\n",logtest1);
 if (count==num) p1=0.000001; //get the the number of how many test1(10000) is less than test(from the user's data)
   else  p1=double(count)/num;   //get the p-value
  if ((p1==0) || (p1==1)) p=1/(1+exp(c0+c1*test));
     else p=p1;
 // Rprintf("num=%d\n",num);
  return p;
  //return test1;
 }

void vexler(int *ndouble1, int *ndouble2, double* y1, double* y2, double* mc , double* test_statistics, double* p_value ){
   vector<double> v_y1;
   vector<double> v_y2;

    for(int i = 0; i < *ndouble1; i++){
        v_y1.push_back(y1[i]);
    }
    for(int i = 0; i < *ndouble2; i++){
        v_y2.push_back(y2[i]);
    }
    *test_statistics=function(v_y1,v_y2);
   // Rprintf("%d",(int)(*mc));
    *p_value = foo(v_y1,v_y2,(int)(*mc));


}


