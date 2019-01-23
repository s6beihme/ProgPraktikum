#pragma once
#include <iostream>
#include <memory>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>


template <typename T>
T convert_to(const std::string &str) {
    std::istringstream ss(str);
    T num;
    ss >> num;
    return num;
}

template<typename T>
T use_partition_to_convert_to_T(std::string s) {
    std::size_t found = s.find("=");
    s.erase(0, found + 1);
    return convert_to<T>(s);
}


template<typename T>
void read_parameters_from_file(std::string filename, int& imax, int& jmax, T& xlength, T& ylength, T& delt, T& t_end, T& tau, T& del_vec, T& eps, T& omg, T& alpha, int& itermax, T& GX, T& GY, T& Re, T& UI, T& VI, T& PI) {
    std::ifstream myfile;
    myfile.open(filename);
    
    if (!myfile.is_open()) {
        std::cout << "FILE DIDNT OPEN!\n";
        exit(0);
    }
    std::string temp;
    
    myfile >> temp;
    imax = use_partition_to_convert_to_T<int>(temp);
    myfile >> temp;
    jmax = use_partition_to_convert_to_T<int>(temp);
    myfile >> temp;
    xlength = use_partition_to_convert_to_T<T>(temp);
    myfile >> temp;
    ylength = use_partition_to_convert_to_T<T>(temp);
    myfile >> temp;
    delt = use_partition_to_convert_to_T<T>(temp);
    myfile >> temp;
    t_end = use_partition_to_convert_to_T<T>(temp);
    myfile >> temp;
    tau = use_partition_to_convert_to_T<T>(temp);
    myfile >> temp;
    del_vec = use_partition_to_convert_to_T<T>(temp);
    myfile >> temp;
    eps = use_partition_to_convert_to_T<T>(temp);
    myfile >> temp;
    omg = use_partition_to_convert_to_T<T>(temp);
    myfile >> temp;
    alpha = use_partition_to_convert_to_T<T>(temp);
    myfile >> temp;
    itermax = use_partition_to_convert_to_T<int>(temp);
    myfile >> temp;
    GX = use_partition_to_convert_to_T<T>(temp);
    myfile >> temp;
    GY = use_partition_to_convert_to_T<T>(temp);
    myfile >> temp;
    Re = use_partition_to_convert_to_T<int>(temp);
    myfile >> temp;
    UI = use_partition_to_convert_to_T<T>(temp);
    myfile >> temp;
    VI = use_partition_to_convert_to_T<T>(temp);
    myfile >> temp;
    PI = use_partition_to_convert_to_T<T>(temp);
}


template <typename T>
T max_eintrag(std::unique_ptr<std::unique_ptr<T[]>[]> &matrix,int zeilen,int spalten) {
    T max_ = 0.0;
    for(int i=0; i<zeilen; i++) {
        for(int j=0; j<spalten; j++) {
            if(abs(matrix[i][j]) > max_)   max_ = matrix[i][j];
        }
    }

    return max_;
}

template <typename T>
void felder_initialisieren(std::unique_ptr<std::unique_ptr<T[]>[]> &matrix,int imax,int jmax,T anfangswert) {
    for(int i=0;i<imax+2;i++) {
        for(int j=0;j<jmax+2;j++) {
            matrix[i][j] = anfangswert;
        }   
    }
}

template <typename T>
T zeitschrittweite(T tau,T Re,T delt,T delx,T dely,std::unique_ptr<std::unique_ptr<T[]>[]> &U,std::unique_ptr<std::unique_ptr<T[]>[]> &V,int zeilen,int spalten) {
    if (tau<0) return delt;
	else if (max_eintrag(U,zeilen,spalten) == 0 && max_eintrag(V,zeilen,spalten) != 0) return tau * std::min((Re/2)*(delx*dely*delx*dely)/((delx*delx) + (dely*dely)), delx/max_eintrag(U,zeilen,spalten));
	else if (max_eintrag(U,zeilen,spalten) != 0 && max_eintrag(V,zeilen,spalten) == 0) return tau * std::min((Re/2)*(delx*dely*delx*dely)/((delx*delx) + (dely*dely)), dely/max_eintrag(V,zeilen,spalten));
	else if (max_eintrag(U,zeilen,spalten) == 0 && max_eintrag(V,zeilen,spalten) == 0) return tau * (Re/2)*(delx*dely*delx*dely)/((delx*delx) + (dely*dely));
	else if (max_eintrag(U,zeilen,spalten) == 0 && max_eintrag(V,zeilen,spalten) == 0) return tau * (Re/2)*(delx*dely*delx*dely)/((delx*delx) + (dely*dely));
	return tau * std::min((Re/2)*(delx*dely*delx*dely)/((delx*delx) + (dely*dely)), std::min(delx/max_eintrag(U,zeilen,spalten), dely/max_eintrag(V,zeilen,spalten)));
}

template <typename T>
void u_randwerte(std::unique_ptr<std::unique_ptr<T[]>[]> &U,int imax,int jmax) {
    for (int j=1;j<jmax+1;j++) {
		U[0][j] = 0.0;
		U[imax][j] = 0.0;
    }
	for (int i=1;i<imax+1;i++) {
		U[i][0] = -U[i][1];
		U[i][jmax+1] = -U[i][jmax];
	}
}

template <typename T>
void v_randwerte(std::unique_ptr<std::unique_ptr<T[]>[]> &V,int imax,int jmax) {
    for (int i=1;i<imax+1;i++) {
		V[i][0] = 0.0;
		V[i][jmax] = 0.0;
    }
	for (int j=1;j<jmax+1;j++) {
		V[0][j] = -V[1][j];
		V[imax+1][j] = -V[imax][j];
	}
}

template <typename T>
void berechne_F(std::unique_ptr<std::unique_ptr<T[]>[]> &F,int imax,int jmax,T delt,T delx,T dely,T Re,std::unique_ptr<std::unique_ptr<T[]>[]> &U,std::unique_ptr<std::unique_ptr<T[]>[]> &V,T GX,T alpha) {
	for(int j=1;j<jmax+1;j++) {      //Randwerte
		F[0][j] = U[0][j];
		F[imax][j] = U[imax][j];
	}
		
	for(int i=1;i<imax;i++) {
		for(int j=1;j<jmax+1;j++) {
			T a = (1/Re)*((U[i+1][j]-2*U[i][j]+U[i-1][j])/ (delx*delx) + (U[i][j+1]-2*U[i][j]+U[i][j-1])/(dely*dely));
			T b = (1/4*delx)*((U[i+1][j]*U[i+1][j]-U[i-1][j]*U[i-1][j]+2*U[i][j]*U[i+1][j]-2*U[i][j]*U[i-1][j]) + alpha*(abs(U[i][j]+U[i+1][j])*(U[i][j]-U[i+1][j])-abs(U[i-1][j]+U[i][j])*(U[i-1][j]-U[i][j])));
			T c = (1/4*dely)*(((V[i][j]+V[i+1][j])*(U[i][j]+U[i][j+1])-(V[i][j-1]+V[i+1][j-1])*(U[i][j-1]+U[i][j]))+alpha*(abs(V[i][j]+V[i+1][j])*(U[i][j]-U[i][j+1])-abs(V[i][j-1]+V[i+1][j-1])*(U[i][j-1]-U[i][j])));
			F[i][j] = U[i][j] + delt*(a-b-c+GX);
		}
	}
}

template <typename T>
void berechne_G(std::unique_ptr<std::unique_ptr<T[]>[]> &G,int imax,int jmax,T delt,T delx,T dely,T Re,std::unique_ptr<std::unique_ptr<T[]>[]> &U,std::unique_ptr<std::unique_ptr<T[]>[]> &V,T GY,T alpha) {
	for(int i=1;i<imax+1;i++) {  //Randwerte
		G[i][0] = V[i][0];
		G[i][jmax] = V[i][jmax];
	}
		
	for(int i=1;i<imax+1;i++) {
		for(int j=1;j<jmax;j++) {
			T a = (1/Re)*((V[i+1][j]-2*V[i][j]+V[i-1][j])/ (delx*delx) + (V[i][j+1]-2*V[i][j]+V[i][j-1])/(dely*dely));
            T b = (1/4*delx)*(((U[i][j]+U[i][j+1])*(V[i][j]+V[i+1][j])-(U[i-1][j]+U[i-1][j+1])*(V[i-1][j]+V[i][j]))+alpha*(abs(U[i][j]+U[i][j+1])*(V[i][j]-V[i+1][j])-abs(U[i-1][j]+U[i-1][j+1])*(V[i-1][j]-V[i][j])));
			T c = (1/4*dely)*((V[i][j+1]*V[i][j+1]-V[i][j-1]*V[i][j-1]+2*V[i][j]*V[i][j+1]-2*V[i][j]*V[i][j-1]) + alpha*(abs(V[i][j]+V[i][j+1])*(V[i][j]-V[i][j+1])-abs(V[i][j-1]+V[i][j])*(V[i][j-1]-V[i][j])));
			G[i][j] = V[i][j] + delt*(a-b-c+GY);
		}
	}
}

template <typename T>
void berechne_rhs(T delt,T delx,T dely,int imax,int jmax,std::unique_ptr<std::unique_ptr<T[]>[]> &F,std::unique_ptr<std::unique_ptr<T[]>[]> &G,std::unique_ptr<std::unique_ptr<T[]>[]> &RHS) {  //Rechte Seite Druckgleichung
	for(int i=1;i<imax+1;i++) {
		for(int j=1;j<jmax+1;j++) {
			RHS[i][j] = (1/delt) * ((F[i][j]-F[i-1][j])/delx + (G[i][j]-G[i][j-1])/dely);
		}
	}
}

template <typename T>
void SOR(int itermax,int imax,int jmax,T omg,T delx,T dely,T eps,std::unique_ptr<std::unique_ptr<T[]>[]> &RHS,std::unique_ptr<std::unique_ptr<T[]>[]> &P, T res) {     
	int it = 0;
	std::unique_ptr<std::unique_ptr<T[]>[]> P_old = std::make_unique<std::unique_ptr<T[]>[]>(imax+2);

	for(int i=0;i<imax+2;i++) {
		P_old[i] = std::make_unique<T[]>(jmax+2);
	}
	
	while (it<itermax && res>=eps) {
		for(int i=0; i<imax+2;i++) {
			for(int j=0;j<jmax+2;j++) {
				P_old[i][j] = P[i][j];
			}
		}
	
		for(int j=1;j<jmax+1;j++) {            //Randwerte
			P[0][j] = P[1][j];
			P[imax+1][j] = P[imax][j];
		}
		for(int i=1;i<imax+1;i++) {
			P[i][0] = P[i][1];
			P[i][jmax+1] = P[i][jmax];
		}	
		for(int i=1;i<imax+1;i++) {
			for(int j=1;j<jmax+1;j++) {
				P[i][j] = (1-omg)*P_old[i][j] + (omg/(2*(1/(delx*delx) + 1/(dely*dely))))*((P_old[i+1][j]+P[i-1][j])/(delx*delx)+(P_old[i][j+1]+P[i][j-1])/(dely*dely)-RHS[i][j]);
			}
		}
		T res = 0.0;
		for(int i=1;i<imax+1;i++) {
			for(int j=1;j<jmax+1;j++) {
				T x = (P[i+1][j]-2*P[i][j]+P[i-1][j])/(delx*delx) + (P[i][j+1]-2*P[i][j]+P[i][j-1])/(dely*dely) - RHS[i][j];
				res += x*x;
			}
		}		
		res = sqrt(res/(imax*jmax));   
        it += 1;
	}
}
	

template <typename T>
void berechne_U(int imax,int jmax,T delt,T delx,std::unique_ptr<std::unique_ptr<T[]>[]> &P,std::unique_ptr<std::unique_ptr<T[]>[]> &F,std::unique_ptr<std::unique_ptr<T[]>[]> &U) {
	for(int i=1;i<imax;i++) { 
		for(int j=1;j<jmax+1;j++) { 
			U[i][j] = F[i][j] - (delt/delx)*(P[i+1][j]-P[i][j]);
		}
	}
}

template <typename T>
void berechne_V(int imax,int jmax,T delt,T dely,std::unique_ptr<std::unique_ptr<T[]>[]> &P,std::unique_ptr<std::unique_ptr<T[]>[]> &G,std::unique_ptr<std::unique_ptr<T[]>[]> &V) {
	for(int i=1;i<imax+1;i++) {
		for(int j=1;j<jmax;j++) {
			V[i][j] = G[i][j] - (delt/dely)*(P[i][j+1]-P[i][j]);
		}
	}
}




template <typename T>
void write_output_into_file(std::string filename, T xlength, T ylength, int imax, int jmax, std::unique_ptr<std::unique_ptr<T[]>[]>& U, std::unique_ptr<std::unique_ptr<T[]>[]>& V, std::unique_ptr<std::unique_ptr<T[]>[]>& P, std::unique_ptr<std::unique_ptr<T[]>[]>& U2, std::unique_ptr<std::unique_ptr<T[]>[]>& V2) {

	for (int i = 1; i < imax + 1; i++) {
		for (int j = 1; j < jmax + 1; j++) {
			U2[i][j] = (U[i][j - 1] + U[i][j]) / 2;
			V2[i][j] = (V[i][j] + V[i + 1][j]) / 2;
		}
	}

	//Daten in File schreiben
	std::ofstream myfile;
	myfile.open(filename);
	if (!myfile.is_open()) {
		std::cout << "ERROR: write_output_into_file failed to open file";
		return;
	}
	myfile << xlength << "\n" << ylength << "\n" << imax << "\n" << jmax << "\n";
	for (int i = 1; i < imax + 1; i++) {
		for (int j = 1; j < jmax + 1; j++) {
			myfile << U2[i][j] << " ";
		}
		myfile << "\n";
	}
	for (int i = 1; i < imax + 1; i++) {
		for (int j = 1; j < jmax + 1; j++) {
			myfile << V2[i][j] << " ";
		}
		myfile << "\n";
	}
	for (int i = 1; i < imax + 1; i++) {
		for (int j = 1; j < jmax + 1; j++) {
			myfile << P[i][j] << " ";
		}
		myfile << "\n";
	}
}

