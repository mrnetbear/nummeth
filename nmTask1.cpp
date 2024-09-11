#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <algorithm>
//#include <math.h>
#include <cmath>
#include <time.h>
#include <random>
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TMatrixD.h"
#include "TDecompLU.h"

int interPoly(){

    //Define the coefficients of the polynomial
    const int n = 7;
    TMatrixD A(n,n);
    TMatrixD f(n,1);
    TMatrixD aj(n,1);
    TMatrixD x(n,1);

    //Set the coefficients

    const double left = 0, right = 3;

    for (int i=0; i<n; i++){
        A(i,0) = 1;
        x(i,0) = right/(double)n * i;

        f(i,0) = sin(x(i,0))*exp(-x(i,0));
        for (int j=1; j<n; j++){
            A(i,j) = pow(x(i,0), j);
        }
    }
    //Print all coefficients
    for (int i=0; i<n; i++){
        for (int j = 0; j < n; j++)
        {
            std::cout << A(i,j) << " ";
        }
        std::cout << "| " << x(i,0) << " | " << f(i,0) << std::endl;
        
    }
    //Perform the transformation
    double det = A.Determinant();
    if (!det) return 1;

    aj = A.Invert(&det)*f;

    //Print inverted coefficients
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++)
        {
            std::cout << A(i,j) << " ";
        }
        std::cout << std::endl;
    }

    //Print found coefficients
    for (int i = 0; i < n; i++){
        std::cout << "aj(" << i << ") = " << aj(i,0) << std::endl;
    }

    double x1 [100];
    double fx1 [100];
    for (int i = 0; i < 100; i++){
        x1[i] = 3.0/100.00 * i;
        fx1[i] = sin(x1[i])*exp(-x1[i]);
    }

    double x2 [100];
    double fx2 [100];
    for (int i = 0; i < 100; i++){
        x2[i] = 3.0/100.00 * i;
        for (int j = 0; j < n; j++){
            fx2[i] += aj(j,0)*pow(x2[i], j);
        }
    }

    TGraph* origin = new TGraph (100, x1, fx1);
    origin->SetMarkerStyle(20);
    origin->SetMarkerColor(kRed);
    origin->SetLineColor(kRed);
    origin->SetLineWidth(2);

    TGraph* approx = new TGraph (100, x2, fx2);
    approx->SetMarkerStyle(20);
    approx->SetMarkerColor(kBlue);
    approx->SetLineColor(kBlue);
    approx->SetLineWidth(2);

    TLegend* legend = new TLegend(0.1,0.7,0.4,0.9);
    origin->Draw();
    approx->Draw("same");
    
    return 0;
}