#include <cppunit/extensions/HelperMacros.h>

#include <gimli.h>
#include <pos.h>
#include <vector.h>
#include <blockmatrix.h>
#include <matrix.h>
#include <sparsematrix.h>
#include <vectortemplates.h>
#include <vector>

#include <stdexcept>

using namespace GIMLI;

class VectorTest : public CppUnit::TestFixture  {
    CPPUNIT_TEST_SUITE(VectorTest);

    CPPUNIT_TEST(testIterator);
    CPPUNIT_TEST(testEquality);
    CPPUNIT_TEST(testFill);
    CPPUNIT_TEST(testSetVal);
    CPPUNIT_TEST(testUnaryOperations);
    CPPUNIT_TEST(testBinaryOperations);
    CPPUNIT_TEST(testExpressionOperators);
    CPPUNIT_TEST(testFunctions);
    CPPUNIT_TEST(testStdVectorTemplates);
    CPPUNIT_TEST(testCVector);
    CPPUNIT_TEST(testRVector3);
    CPPUNIT_TEST(testMatrix);
    CPPUNIT_TEST(testBlockMatrix);
    CPPUNIT_TEST(testSparseMapMatrix);
    CPPUNIT_TEST(testFind);
    CPPUNIT_TEST(testIO);

    CPPUNIT_TEST_SUITE_END();

public:

    void testIterator(){
        RVector v(10); v.fill(x__ + 1.0);
        testIterator_(v);

        R3Vector v1(10);
        testIterator_(v1);
    }

    template < class Vec > void testIterator_(Vec & v){
        GIMLI::Index len = v.size();
        GIMLI::VectorIterator<typename Vec::ValType> iter = v.begin();
        uint i = 0;
        do {
            i++;
//             __MS(i << " " << v[i-1])
            CPPUNIT_ASSERT(v[i-1] == iter.nextVal());
        } while (iter.hasMore());
        CPPUNIT_ASSERT(i == len);

        CPPUNIT_ASSERT(v(0,len).beginPyIter()[0] == v[0]);
        i = 0;
        for (typename Vec::iterator it = v.begin(); it != v.end(); it++, i ++){
//             __MS(v[i] << " " << *it)
            CPPUNIT_ASSERT(v[i] == *it);
        }
    }

    void testEquality(){
        testEquality_< RVector > ();
        testEquality_< CVector > ();
    }

    template < class Vec > void testEquality_(){
        Vec v1 = Vec(10);

        try{ v1.getVal(0); CPPUNIT_ASSERT(1); } catch(...){}
        try{ v1.getVal(10); CPPUNIT_ASSERT(0); } catch(...){}

        Vec v2 = Vec(v1);
        CPPUNIT_ASSERT(v1 == v1);
        CPPUNIT_ASSERT(v1 == v2);

        Vec v;          CPPUNIT_ASSERT(v1 != v);
        v = v1;         CPPUNIT_ASSERT(v1 == v);
        v = v1 + 1.0;   CPPUNIT_ASSERT(v1 != v);
        v = v1 + 1.0;   CPPUNIT_ASSERT(v == v1 + 1.0);
    }

    void testFill(){
        double vec[5]={0.0, 1.0, 2.0, 3.0, 4.0};

        RVector v1(5); v1.fill(vec);
        Placeholder x__;
        RVector v2(5); v2.fill(x__);
        CPPUNIT_ASSERT(v1 == v2);
        CPPUNIT_ASSERT(fliplr(v1)[0] == 4.0);
        CPPUNIT_ASSERT(fliplr(v1)[4] == 0.0);
        CPPUNIT_ASSERT(fliplr(fliplr(v1)) == v1);
    }

    void testUnaryOperations(){
        testUnaryOperations_< double > (3.1415);
        testUnaryOperations_< Complex > (Complex(3.1415, 1.0));
        RVector v(10, 42.42);
        RVector v1(sign(square(v))); CPPUNIT_ASSERT(sum(v1) == 10);
        randn(v);
        CPPUNIT_ASSERT(RVector(abs(v + 1e-6)) == RVector(exp10(log10(abs(v + 1e-6)))));
        CPPUNIT_ASSERT(RVector(abs(v + 1e-6)) == RVector(exp(log(abs(v + 1e-6)))));
        CPPUNIT_ASSERT(RVector(abs(v + 1e-6)) != RVector(exp10(log(abs(v + 1e-6)))));
    }
    template < class ValueType > void testUnaryOperations_(const ValueType & fac){
        typedef Vector < ValueType > Vec;
        Vec v1 = Vec(10);
        Vec v(v1), v2(v1);
        v.fill(fac);
        v2.fill(fac);
        v += fac; CPPUNIT_ASSERT(v == Vec(v1.size(), fac * 2.0));
        v -= fac; CPPUNIT_ASSERT(v == Vec(v1.size(), fac));
        v *= fac; CPPUNIT_ASSERT(v == Vec(v1.size(), fac * fac));
        v /= fac; CPPUNIT_ASSERT(v == Vec(v1.size(), fac));
        v += v2;  CPPUNIT_ASSERT(v == Vec(v1.size(), fac * 2.0));
        v -= v2;  CPPUNIT_ASSERT(v == Vec(v1.size(), fac));
        v *= v2;  CPPUNIT_ASSERT(v == Vec(v1.size(), fac * fac));
        v /= v2;  CPPUNIT_ASSERT(v == Vec(v1.size(), fac));
        v = -v2;  CPPUNIT_ASSERT(v == v2 * -1.);
    }

    void testExpressionOperators(){

        RVector m(10, 2.0);
        RVector t(10, 2.0);

        RVector t1((m * m) / t[0]);
        RVector t2(RVector(m * m) / t[0]);

//         __MS(t1)
//         __MS(t2)
        CPPUNIT_ASSERT(t1 == t2);
//         exit(0);
    }

    void testSetVal(){
        typedef Vector < double > Vec;
        Vec v1(10);
        v1.fill(x__ + 1.0);
        Vec v2(10, 1.0);
        Vec v3(10);

        CPPUNIT_ASSERT_THROW(v1.setVal(1, -1), std::out_of_range);
        CPPUNIT_ASSERT_THROW(v1.setVal(1, 11), std::out_of_range);

        //** setVal(const ValueType & val, const BVector & bv)
        CPPUNIT_ASSERT(sum(v2.setVal(2, BVector(v2.size(), true))) == 20);
        CPPUNIT_ASSERT(sum(v2.setVal(1, BVector(v2.size(), false))) == 20);
        CPPUNIT_ASSERT(sum(v2.setVal(1, BVector(v2.size(), true))) == 10);
        CPPUNIT_ASSERT_THROW(v1.setVal(1, BVector(5)), std::length_error);

        CPPUNIT_ASSERT(sum(v3.setVal(v1, 0, 10)) == 55);
        CPPUNIT_ASSERT(sum(v3.addVal(v2, 5, 10)) == 60);
        CPPUNIT_ASSERT(sum(v3.setVal(1.0)) == 10);

        CPPUNIT_ASSERT(sum(v3.setVal(2.0, find(v1 > 5))) == 15);

        RVector v4(5, 3.0);
        CPPUNIT_ASSERT(sum(v3.setVal(v4, 5, 10)) == 20);
        CPPUNIT_ASSERT(sum(v3.addVal(v4, 5, 10)) == 35);

        //** setVal(const ValueType & val, Index start, SIndex end)
        // set from 0 to end == 1
        CPPUNIT_ASSERT(sum(v1.setVal(1.0, 0, -1)) == 10);
        // set from 5 to end == 0
        CPPUNIT_ASSERT(sum(v1.setVal(0.0, 5, 100)) == 5);
        // set from 5 to 7 == 1
        CPPUNIT_ASSERT(sum(v1.setVal(1.0, 5, 7)) == 7);
        // set from -1 to 7 -> from 7 to 7 == 1
        CPPUNIT_ASSERT(sum(v1.setVal(1.0, -1, 7)) == 7);
        // set from 10 to 7 -> from 7 to 7 == 1
        CPPUNIT_ASSERT(sum(v1.setVal(1.0, 10, 7)) == 7);
        // set from 0 to 17 -> from 0 to 10 == 1
        CPPUNIT_ASSERT(sum(v1.setVal(1.0, 0, 17)) == 10);
        // set from 0 to -2 -> from 0 to 10 == 1
        CPPUNIT_ASSERT(sum(v1.setVal(1.0, 0, -2)) == 10);

        //**  setVal(const Vector < ValueType > & vals, const IndexArray & iArray)

        CPPUNIT_ASSERT(sum(v1.setVal(v4(find(v4 > 0)), find(v4 > 0))) == 20);
        // set v1[1]=2
        CPPUNIT_ASSERT(v1.setVal(RVector(1,2), IndexArray(1,1))[1] == 2);
        CPPUNIT_ASSERT(sum(v1.setVal(RVector(1,2), IndexArray(1,1))) == 19);

        CPPUNIT_ASSERT_THROW(v1.setVal(RVector(1,1), IndexArray(1,10)),
                              std::out_of_range);
        CPPUNIT_ASSERT_THROW(v1.setVal(v4(find(v4 > 0)), find(v4 > 3)),
                             std::length_error);

    }

    void testBinaryOperations(){
        RVector v1(*v1_);
        RVector v2(*v1_);
        double fac = 3.1415;
        v1.fill(fac); v2.fill(fac);
        CPPUNIT_ASSERT(RVector(v1 + v2) == RVector(v1_->size(), 2.0 * fac));
        CPPUNIT_ASSERT(RVector(v1 - v2) == RVector(v1_->size(), 0.0));
        CPPUNIT_ASSERT(RVector(v2 - v1) == RVector(v1_->size(), 0.0));
        CPPUNIT_ASSERT(RVector(v1 * v2) == RVector(v1_->size(), fac * fac));
        CPPUNIT_ASSERT(RVector(v1 / v2) == RVector(v1_->size(), 1.0));
        CPPUNIT_ASSERT(RVector(v2 / v1) == RVector(v1_->size(), 1.0));

        v2 *= 2.0;
        CPPUNIT_ASSERT(RVector(v1 - v2) == RVector(v1_->size(), -fac));
        CPPUNIT_ASSERT(RVector(v2 - v1) == RVector(v1_->size(),  fac));
        CPPUNIT_ASSERT(RVector(v1 - v2) != RVector(v1_->size(),  fac));
        CPPUNIT_ASSERT(RVector(v1 / v2) == RVector(v1_->size(),  0.5));
        CPPUNIT_ASSERT(RVector(v2 / v1) == RVector(v1_->size(),  2.0));
    }

    void testFunctions(){
        typedef Vector < double > Vec;

        Vec v(*v1_);
        v.fill(x__ + 1.0);
        CPPUNIT_ASSERT(sum(v) == (v[0] + v[v1_->size() -1]) *
                                                (::floor(v1_->size() / 2.0) +
                                                    (v1_->size() / 2.0 - ::floor(v1_->size() / 2.0))));
        CPPUNIT_ASSERT(GIMLI::sum(v - v) == 0.0);


        CPPUNIT_ASSERT(GIMLI::sum(v(0, -1) - v(0, v.size()-1)) == 0.0);
        CPPUNIT_ASSERT(min(v) == 1.0);
        CPPUNIT_ASSERT(max(v) == v1_->size());
        RVector vs(10);
        vs[1] = 1; vs[2] = 1; vs[3] = 2; vs[4] = 2; vs[7] = 1;
        CPPUNIT_ASSERT(unique(vs).size() == 6); // 0 1 2 0 1 0
        CPPUNIT_ASSERT(unique(sort(vs)).size() == 3); // 0 1 2

        vs.fill(x__ + 1.0);
        vs *= 1.111;
        vs.round(0.01);
        CPPUNIT_ASSERT(::fabs(vs[0] - 1.11) < TOLERANCE);
        CPPUNIT_ASSERT(::fabs(vs[4] - 5.56) < TOLERANCE);
        CPPUNIT_ASSERT(::fabs(vs[6] - 7.78) < TOLERANCE);
        CPPUNIT_ASSERT(::fabs(vs[9] - 11.11) < TOLERANCE);

        vs *=0;
        CPPUNIT_ASSERT(nonZero(vs) == false);
        CPPUNIT_ASSERT(zero(vs) == true);

        Vec r;
        r = GIMLI::increasingRange(0.1, 1.0, 3);
        r = GIMLI::increasingRange(0.1, 3.3, 5);
    }

    void testStdVectorTemplates(){
        std::vector < int > iVec(10, 0);
        iVec[1] = 1; iVec[2] = 1; iVec[3] = 2; iVec[4] = 2; iVec[7] = 1;
        CPPUNIT_ASSERT(unique(iVec).size() == 6); // 0 1 2 0 1 0
        CPPUNIT_ASSERT(unique(sort(iVec)).size() == 3); // 0 1 2

        std::vector < double > dVec(10, 1.0);
        CPPUNIT_ASSERT(GIMLI::max(dVec) == 1.0);
    }

    void testRVector3(){
        RVector3 p1(1.0, 1.0, 1.0);
        RVector3 p2(-1.0, -1.0, -1.0);
        CPPUNIT_ASSERT(p1 == p1);
        CPPUNIT_ASSERT(p1 != p2);
        CPPUNIT_ASSERT(p1 == -p2);
        CPPUNIT_ASSERT(p1 != p2);
        CPPUNIT_ASSERT(p1 != p1 * -1.0);
        CPPUNIT_ASSERT(p1 == (p2 * -1.0));
        CPPUNIT_ASSERT(p1 == p1 / 1.0);
        CPPUNIT_ASSERT(p1 * -1.0 == p1 / -1.0);
        CPPUNIT_ASSERT(p1 + p2 == p2 + p1);
        p2 *= -1.0;
        CPPUNIT_ASSERT(p1 == p2);
        p2 = p1; p2 += p1;  CPPUNIT_ASSERT(p1 + p1 == p2);
        p2 = p1; p2 += 1.0;  CPPUNIT_ASSERT(p1 + p1 == p2);
        p2 = p1; p2.translate(p1);  CPPUNIT_ASSERT(p1 + p1 == p2);
        p2 = p1; p2.scale(p1);  CPPUNIT_ASSERT(p1 == p2);
    }

    void testCVector(){
        size_t dim = 10;
//         RVector freq(dim);
//         Complex i_unit(0.0 , 1.0);
//         CVector k1(dim) ;
//         CVector g(dim, Complex(1.0, 0.0));
//         CVector z = toComplex(freq) / (k1 * g) * i_unit;

        CVector c1(dim, Complex(1.0, 0.0));
        RVector r1(dim, 1.0);
        CVector c2;
        c2 = c1 * r1;
        CPPUNIT_ASSERT(c1 == c2);
        c2 = r1 * c1;
        CPPUNIT_ASSERT(c1 == c2);
        c2 = toComplex(real(c1), imag(c1));
        CPPUNIT_ASSERT(c1 == c2);
        CPPUNIT_ASSERT(real(c1 + conj(c1)) == RVector(2.0 * real(c1)));
        CPPUNIT_ASSERT(imag(c1 - conj(c1)) == RVector(2.0 * imag(c1)));

        r1 *= 0.0;
        CPPUNIT_ASSERT(imag(c1 * conj(c1)) == r1);
        CPPUNIT_ASSERT(RVector(angle(c1) + angle(conj(c1))) == r1);
        CPPUNIT_ASSERT(abs(c1) == RVector(sqrt(real(c1) * real(c1)
                                                  + imag(c1) * imag(c1))));
        CPPUNIT_ASSERT(abs(c1) == RVector(exp(real(log(c1)))));
        //real, imag, angle, conj, abs

        c2 = toComplex(real(c1), real(c1));
        CVector c3(1./c2);
        CPPUNIT_ASSERT(c3[0] == Complex(0.5, -0.5));
        CPPUNIT_ASSERT(c3[0] == 1./ c2[0]);
    }

    template < class ValueType > void testMatrix_(){
        typedef Matrix < ValueType > Mat;
        typedef Vector < ValueType > Vec;

        Mat A(5, 5);
#if not defined ( __APPLE__ )
        try{ A.row(11); CPPUNIT_ASSERT(0); } catch(...){}
        try{ A.row(-1); CPPUNIT_ASSERT(0); } catch(...){}
        try{ CPPUNIT_ASSERT(0); } catch(...){}
#else
	std::cout << "check! cppunit exception check for macos " << std::endl;
#endif
        CPPUNIT_ASSERT(A.rows() == 5);
        CPPUNIT_ASSERT(A.cols() == 5);

        A.resize(3, 2); CPPUNIT_ASSERT(A.cols() == 2);
        A.resize(8, 9); CPPUNIT_ASSERT(A.cols() == 9);

        A[0][0] = 1.0;
        A[1] = A[0];

        CPPUNIT_ASSERT(A[0] == A[1]);
        CPPUNIT_ASSERT(A.row(2) != A[1]);
        CPPUNIT_ASSERT(A[0][0] == 1.0);
        CPPUNIT_ASSERT(A[1][0] == 1.0);

        CPPUNIT_ASSERT(fliplr(fliplr(A)) == A);

        A.push_back(A[0]); CPPUNIT_ASSERT(A.rows() == 9);
        CPPUNIT_ASSERT(A[A.rows()-1] == A.back());

        Mat B(A);
        CPPUNIT_ASSERT(B.rows() == 9);
        CPPUNIT_ASSERT(B.cols() == 9);
        CPPUNIT_ASSERT(B[0][0] == 1.0);

        B[1][1] = 1.0;
        A = B;
        CPPUNIT_ASSERT(A == B);
        CPPUNIT_ASSERT(A[1][1] == 1.0);

        A.clear();
        CPPUNIT_ASSERT(A.rows() == 0);
        CPPUNIT_ASSERT(A.cols() == 0);

        A.resize(3, 4);
        A *= 0.0;
        A += 1.0;
        CPPUNIT_ASSERT(sum(A[0]) == A.cols());

        Vec a(A.rows()); a.fill(x__ * 2.5);
        Vec b(A.cols()); b.fill((x__ + 1.0) * 4.2);

        B = A;
        Mat C(B);
        scaleMatrix(B, a, b);
        rank1Update(C, a, b);
        CPPUNIT_ASSERT((C-B) == A);
        CPPUNIT_ASSERT(::fabs(sum(A * b) - sum(b) * A.rows()) < TOLERANCE);

        A *= 1.11;
        A.round(0.1);
        CPPUNIT_ASSERT(::fabs(A[0][0] -1.1) < TOLERANCE);

        A *= 0.;
        A += 2.;

        // test col, row access
        A[1].fill(2);
        CPPUNIT_ASSERT(sum(A.row(1)) == A.cols()*2);
        CPPUNIT_ASSERT(sum(A.col(1)) == A.rows()*2);
    }

    void testFind(){
        typedef Vector < double > Vec;
        Vec x(10); x.fill(x__ + 1.0);
        Vec y(10); y.fill(x__ + 2.0);

        x[5] = 10.0;

        BVector b(y > x);
        CPPUNIT_ASSERT((y > x).size() == 10);
        CPPUNIT_ASSERT((y < x).size() == 10);
        CPPUNIT_ASSERT(find(y > x).size() == 9);
        CPPUNIT_ASSERT(find(x > y).size() == 1);
        CPPUNIT_ASSERT(find(x > y)[0] == 5);
        CPPUNIT_ASSERT(x(find(x > y))[0] == 10.0);

        x[2] = 5.0;
        x[3] = 5.0;
	// x = [1, 2, 5, 5, 5, 6, 7, 8, 9, 10]
        CPPUNIT_ASSERT((x == 5).size() == 10);
        CPPUNIT_ASSERT(find(x == 5).size() == 3);
        CPPUNIT_ASSERT(find(~(x == 5)).size() == 7);
        CPPUNIT_ASSERT(find(x == 5)[0] == 2);
        CPPUNIT_ASSERT(find(x <= 5).size() == 5);
        CPPUNIT_ASSERT(find(x > 5).size() == 5);
        CPPUNIT_ASSERT(find(x < 5).size() == 2);
        CPPUNIT_ASSERT(find((x > 5) & (x < 5)).size() == 0);
        CPPUNIT_ASSERT(find((x > 5) | (x < 5)).size() == 7);
        CPPUNIT_ASSERT(find(x+x == 4.).size() == 1);
        CPPUNIT_ASSERT(find(x+x >= 4.).size() == 9);
        CPPUNIT_ASSERT(find(x+x > 4.).size() == 8);
        CPPUNIT_ASSERT(find(abs(x+x) > 4.).size() == 8);
        CPPUNIT_ASSERT(find(~(abs(x+x) > 4.)).size() == 2);
        x[0] = ::log(0);
        CPPUNIT_ASSERT(find(isInf(x)).size() == 1);
        CPPUNIT_ASSERT(find(isInf(x+x)).size() == 1);
        x[1] += ::sqrt(-1) ;
        CPPUNIT_ASSERT(find(isNaN(x)).size() == 1);
        CPPUNIT_ASSERT(find(isInfNaN(x)).size() == 2);
        //CPPUNIT_ASSERT(find(isnan(x+x)).size() == 0);

        GIMLI::CVector c(10);
        CPPUNIT_ASSERT(find(c < Complex(.0, 0.0)).size() == 0);
//         std::less< std::complex<double> > l;
//         std::cout << l(Complex(.0, 0.0), Complex(.0, 0.0)) << std::endl;
    }

    void testMatrix(){
        testMatrix_< double >();
        testMatrixMult();
//        testMatrix_< float >();
    }

    void testMatrixMult(){

        // GIMLI::matMult(A, B, C, -1);
        // CPPUNIT_ASSERT(C[0] == GIMLI::RVector(std::vector< double >{0., 0.}));
        // CPPUNIT_ASSERT(C[1] == GIMLI::RVector(std::vector< double >{0., 0.}));

        Index m = 2;
        Index n = 3;
        Index k = 4;

        double *_A = new double[m * k];
        double *_B = new double[k * n];
        for (Index i = 0; i < m*k; i ++ ){_A[i] = i+1;}
        for (Index i = 0; i < k*n; i ++ ){_B[i] = i+1;}

        GIMLI::RMatrix A(m, k, _A);
        GIMLI::RMatrix B(k, n, _B);
        GIMLI::RMatrix C;
        GIMLI::matMult(A, B, C, 1.0, 0.0);
        CPPUNIT_ASSERT(C.rows() == m);
        CPPUNIT_ASSERT(C.cols() == n);
        CPPUNIT_ASSERT(C[0] ==
                       GIMLI::RVector(std::vector< double >{70., 80., 90}));
        CPPUNIT_ASSERT(C[1] ==
                       GIMLI::RVector(std::vector< double >{158, 184, 210}));

        double *_AT = new double[k*m];

        GIMLI::RMatrix AT(k, m, _A);
        for (Index i = 0; i < k; i ++ ){
            for (Index j = 0; j < m; j ++ ){
                AT[i][j] = A[j][i];
            }
        }

        GIMLI::RMatrix C2;
        GIMLI::matTransMult(AT, B, C2, 1.0, 0.0);
        CPPUNIT_ASSERT(C2.rows() == m);
        CPPUNIT_ASSERT(C2.cols() == n);
        CPPUNIT_ASSERT(C2[0] ==
                       GIMLI::RVector(std::vector< double >{70., 80., 90}));
        CPPUNIT_ASSERT(C2[1] ==
                       GIMLI::RVector(std::vector< double >{158, 184, 210}));
    }

    void testBlockMatrix(){
        GIMLI::BlockMatrix < double > A(false);

        GIMLI::RMatrix B;
        GIMLI::Index m1 = A.addMatrix(&B);

        GIMLI::RSparseMapMatrix C;
        GIMLI::Index m2 = A.addMatrix(&C);

        A.addMatrixEntry(m1, 0, 0);
        A.addMatrixEntry(m1, 3, 3);
        A.addMatrixEntry(m1, 3, 0);
        A.addMatrixEntry(m2, 6, 0, -1.0);

        B.resize(3, 3);
        B += 1.0;

        C.resize(3, 3);
        for (uint i = 0; i < 3; i ++ ){
            for (uint j = 0; j < 3; j ++ ){
                C.addVal(i,j,1);
            }
        }

        CPPUNIT_ASSERT(A.rows() == 9);
        CPPUNIT_ASSERT(A.cols() == 6);

        GIMLI::RVector b(A.cols(), 1);
        CPPUNIT_ASSERT(sum(A.mult(b)) == 18);
        GIMLI::RVector c(A.rows(), 1);
        CPPUNIT_ASSERT(sum(A.transMult(c)) == 18);
    }

    void testSparseMapMatrix(){
        GIMLI::RSparseMapMatrix A(2, 2);
        A.addVal(0, 0, 1.0);
        A.addVal(1, 1, 1.0);
        GIMLI::RSparseMapMatrix B;
        B = A;
        GIMLI::RSparseMapMatrix C(B+A);
        CPPUNIT_ASSERT(C.getVal(0, 0) == 2.0);
        CPPUNIT_ASSERT((C+C)[0][0] == 4.0);
        CPPUNIT_ASSERT((C+C).getVal(0, 0) == 4.0);
        CPPUNIT_ASSERT((C*2.0).getVal(0, 0) == 4.0);
        CPPUNIT_ASSERT(((C+C)*2.0).getVal(1, 1) == 8.0);

        GIMLI::RSparseMapMatrix D(3, 4);
        for (GIMLI::Index i = 0; i < D.rows(); i ++ ){
            for (GIMLI::Index j = 0; j < D.cols(); j ++ ){
                D[i][j] = 1.0;
            }
        }
        CPPUNIT_ASSERT(D.col(2) == GIMLI::RVector(D.rows(), 1.0));
        CPPUNIT_ASSERT(D.row(2) == GIMLI::RVector(D.cols(), 1.0));

        D.cleanRow(1);
        CPPUNIT_ASSERT(D.col(2) == GIMLI::RVector(std::vector< double >{1., 0., 1.}));
    }

    void testIO(){
        RVector v(100);
        randn(v);
        save(v, "test.vec");
        RVector t1("test.vec");
        CPPUNIT_ASSERT(v == t1);

        save(t1, "test", Binary);
        RVector t2("test", Binary);
        CPPUNIT_ASSERT(v == t2);
    }
    void setUp(){
        v1_ = new RVector(10);
        v2_ = new RVector(*v1_);
    }

    void tearDown(){
        delete v1_;
        delete v2_;
    }

private:
    RVector * v1_, * v2_;

};


CPPUNIT_TEST_SUITE_REGISTRATION(VectorTest);
