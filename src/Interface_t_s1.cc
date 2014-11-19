#include "model_t_s1.h"

//void test1(int *ndouble1, int *ndouble2, double* y1, double* y2, double* test_statistics){
//    test(ndouble1, ndouble2, y1, y2, test_statistics);
//}

//void vexler1(int *ndouble1, int *ndouble2, double* y1, double* y2, double* mc, double* test_statistics, double* p_value){
//    vexler(ndouble1, ndouble2, y1, y2, mc, test_statistics, p_value);
//}


extern "C" {
    void CWrapper1(int *ndouble1, int *ndouble2, double* y1, double* y2, double* test_statistics){
        test(ndouble1, ndouble2, y1, y2, test_statistics);
    }
}
extern "C" {
    void CWrapper(int *ndouble1, int *ndouble2, double* y1, double* y2, double* mc, double* test_statistics, double* p_value){
    vexler(ndouble1, ndouble2, y1, y2, mc, test_statistics,p_value);
    }
}

