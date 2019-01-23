#include "driven_cavity.h"

int main(int argc, char* argv[]) {
    
    if (argc < 3) {
        std::cout << "Programm braucht Name des Input und des Output Files";
        return 0;
    }
    
    std::string inputname = argv[1];
    std::string outputname = argv[2];
    
    int imax; int jmax; double xlength; double ylength; double delt; double t_end; double tau; double del_vec; double eps; double omg; double alpha; int itermax = 0; double GX; double GY; double Re; double UI; double VI; double PI;
    read_parameters_from_file(inputname, imax, jmax, xlength, ylength, delt, t_end, tau, del_vec, eps, omg, alpha, itermax, GX, GY, Re, UI, VI, PI);
    
    
    
    std::unique_ptr<std::unique_ptr<double[]>[]> U = std::make_unique<std::unique_ptr<double[]>[]>(imax + 2);
    for (int i = 0; i < imax + 2; i++) {
        U[i] = std::make_unique<double[]>(jmax + 2);
    }
    std::unique_ptr<std::unique_ptr<double[]>[]> V = std::make_unique<std::unique_ptr<double[]>[]>(imax + 2);
    for (int i = 0; i < imax + 2; i++) {
        V[i] = std::make_unique<double[]>(jmax + 2);
    }
    std::unique_ptr<std::unique_ptr<double[]>[]> U2 = std::make_unique<std::unique_ptr<double[]>[]>(imax + 2);
    for (int i = 0; i < imax + 2; i++) {
        U2[i] = std::make_unique<double[]>(jmax + 2);
    }
    std::unique_ptr<std::unique_ptr<double[]>[]> V2 = std::make_unique<std::unique_ptr<double[]>[]>(imax + 2);
    for (int i = 0; i < imax + 2; i++) {
        V2[i] = std::make_unique<double[]>(jmax + 2);
    }
    std::unique_ptr<std::unique_ptr<double[]>[]> P = std::make_unique<std::unique_ptr<double[]>[]>(imax + 2);
    for (int i = 0; i < imax + 2; i++) {
        P[i] = std::make_unique<double[]>(jmax + 2);
    }
    std::unique_ptr<std::unique_ptr<double[]>[]> F = std::make_unique<std::unique_ptr<double[]>[]>(imax + 2);
    for (int i = 0; i < imax + 2; i++) {
        F[i] = std::make_unique<double[]>(jmax + 2);
    }
    std::unique_ptr<std::unique_ptr<double[]>[]> G = std::make_unique<std::unique_ptr<double[]>[]>(imax + 2);
    for (int i = 0; i < imax + 2; i++) {
        G[i] = std::make_unique<double[]>(jmax + 2);
    }
    std::unique_ptr<std::unique_ptr<double[]>[]> RHS = std::make_unique<std::unique_ptr<double[]>[]>(imax + 2);
    for (int i = 0; i < imax + 2; i++) {
        RHS[i] = std::make_unique<double[]>(jmax + 2);
    }
    
    
    
    std::string appendix;
    int number_of_files = 0;
    double del_vec2 = del_vec;
    double res = 100000;
    double delx = xlength / imax;
    double dely = ylength / jmax;
    
    felder_initialisieren(U,imax,jmax,UI);
    felder_initialisieren(V,imax,jmax,VI);
    felder_initialisieren(U2,imax,jmax,0.0);
    felder_initialisieren(V2,imax,jmax,0.0);
    felder_initialisieren(P,imax,jmax,PI);
    felder_initialisieren(F,imax,jmax,0.0);
    felder_initialisieren(G,imax,jmax,0.0);
    felder_initialisieren(RHS,imax,jmax,0.0);
    double t = 0;
    
    while (t < t_end) {
        double delt= zeitschrittweite(tau,Re,delt,delx,dely,U,V,imax+2,jmax+2);
        
        u_randwerte(U,imax,jmax);
        v_randwerte(V,imax,jmax);
        
        for (int i=1;i<imax+1;i++) {
            U[i][jmax+1] = 2.0 - U[i][jmax];
        }  
        berechne_F(F,imax,jmax,delt,delx,dely,Re,U,V,GX,alpha);
        berechne_G(G,imax,jmax,delt,delx,dely,Re,U,V,GY,alpha);
            
        berechne_rhs(delt,delx,dely,imax,jmax,F,G,RHS);
            
        SOR(itermax,imax,jmax,omg,delx,dely,eps,RHS,P,res);
            
        berechne_U(imax,jmax,delt,delx,P,F,U);
        berechne_V(imax,jmax,delt,dely,P,G,V);
        
    
        if (t > del_vec) {
            number_of_files += 1;
            
            if(number_of_files<10) appendix = "_00"+ std::to_string(number_of_files);
            else {
                if(number_of_files<100) appendix = "_0"+ std::to_string(number_of_files);
                else appendix = "_"+std::to_string(number_of_files);
            }
            write_output_into_file(outputname+appendix, xlength, ylength, imax, jmax, U, V, P, U2, V2);
            del_vec = t + del_vec2;
        }
        
        t += delt;
    }
    number_of_files += 1;
    
    if (number_of_files < 10) appendix = "_00" + std::to_string(number_of_files);
    else {
        if (number_of_files < 100) appendix = "_0" + std::to_string(number_of_files);
        else appendix = "_"+std::to_string(number_of_files);
    }
    write_output_into_file(outputname+appendix, xlength, ylength, imax, jmax, U, V, P, U2, V2);
}
