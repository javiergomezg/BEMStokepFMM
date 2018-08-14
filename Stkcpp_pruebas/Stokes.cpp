/*Prueba del Stresslet en un cilindro con C++*/

#include"cabecera.hpp"

using namespace std;
using namespace Eigen;


/*Codigo para generar un cilindro:
	Rad_cil = radio del cilindro
	h_teta = diferencial de angulo
	n_pan = numero de paneles en que se divide el cilindro
*/

/*
double geom_cil(int, double, double);

double geom_cil(int n, double R, double h){
	double x_p[n+1], y_p[n+1];
	double x_c[n], y_c[n];

	for(int i = 0; i <= n; i++){
		x_p[i] = R*cos(i*h);
		y_p[i] = R*sin(i*h);
		cout<<i<< " x = " << x_p[i] << ", y = " << y_p[i] << endl;
		//cout<< "("<<location[i].x<<","<<location[i].y<<")"<<endl; // visualizacion del contenido de location
	}

	cout << "\n Nodos de colocacion "<<endl;

	for(int i = 0; i < n; i++){
		x_c[i] = 0.5*(x_p[i] + x_p[i+1]);
		y_c[i] = 0.5*(y_p[i] + y_p[i+1]);
		cout<<i<< " x = " << x_c[i] <<endl ;//", y = " << y_col[i] << endl;

	}
	return x_p, y_p, x_c, y_c;
}
*/

class myKernelTxxx: public kernel_Base {
public:
    virtual double kernel_Func(Point r0, Point r1){
        //implement your own kernel here
	double dx = (r0.x - r1.x);
	double dy = (r0.x - r1.x);
        double rSquare	= dx*dx + dy*dy;
	double r = pow(rSquare,0.5);
        return -4.*dx*dx*dx/pow(r,4.);
    }
};

class myKernelTxyx: public kernel_Base { //sirve con Txxy y Tyxx
public:
    virtual double kernel_Func(Point r0, Point r1){
        //implement your own kernel here
	double dx = (r0.x - r1.x);
	double dy = (r0.x - r1.x);
        double rSquare	= dx*dx + dy*dy;
	double r = pow(rSquare,0.5);
        return -4.*dx*dy*dx/pow(r,4.);
    }
};

class myKernelTyxy: public kernel_Base { //sirve con Txyy y Tyyx
public:
    virtual double kernel_Func(Point r0, Point r1){
        //implement your own kernel here
	double dx = (r0.x - r1.x);
	double dy = (r0.x - r1.x);
        double rSquare	= dx*dx + dy*dy;
	double r = pow(rSquare,0.5);
        return -4.*dy*dx*dy/pow(r,4.);
    }
};

class myKernelTyyy: public kernel_Base {
public:
    virtual double kernel_Func(Point r0, Point r1){
        //implement your own kernel here
	double dx = (r0.x - r1.x);
	double dy = (r0.x - r1.x);
        double rSquare	= dx*dx + dy*dy;
	double r = pow(rSquare,0.5);
        return -4.*dy*dy*dy/pow(r,4.);
    }
};


int main(){
	int n_pan, n_cuad;
	double pi = acos(-1.), dpi = 2.*pi, ddpi = 4.*pi;
	double Rad_cil, h_theta ;
	double u_0; // velocidad del flujo libre	

	//cout<<"tamaño location: "<<location.size()<< endl; // con este comando revizo el tamaño del arreglo tipo vector<Point> 

	cout << "Digite numero de paneles: ";
	cin >> n_pan;
	
	cout << "Digite el radio del cilindro: ";
	cin >> Rad_cil;

	cout << "Digite la velocidad del flujo de corriente libre: ";
	cin >> u_0;

	cout << "Digite nodos de cuadratura: ";
	cin >> n_cuad;

	h_theta = dpi/n_pan;

	
	//cout << "\npi = " << endl << pi << endl;
	//cout << "2pi = " << endl << dpi << endl;
	//cout << "h_theta = "<< h_theta << endl;

	/*############################################
	##Inicio definicion del cilindro##############
	#############################################*/

	double x_pan[n_pan + 1], y_pan[n_pan + 1];
	double x_col[n_pan], y_col[n_pan];

	double xc = 0.0, yc = 0.0; // centro del cilindro

	cout << "Digite el centro del cilindro: (xc,yc) "<<endl;
	cin>>xc>>yc;

	for(int i = 0; i <= n_pan; i++){
		x_pan[i] = Rad_cil*cos(i*h_theta) + xc;
		y_pan[i] = Rad_cil*sin(i*h_theta) + yc;
		//cout<<i<< " x = " << x_pan[i] << ", y = " << y_pan[i] << endl;
		//cout<< "("<<location[i].x<<","<<location[i].y<<")"<<endl; // visualizacion del contenido de location
	}

	//cout << "\n Nodos de colocacion "<<endl;

	for(int i = 0; i < n_pan; i++){
		x_col[i] = 0.5*(x_pan[i] + x_pan[i+1]);
		y_col[i] = 0.5*(y_pan[i] + y_pan[i+1]);
		//cout<<i<< " x = " << x_col[i] << ", y = " << y_col[i] << endl;

	}

	double R_pan[n_pan]; // Longitud de los paneles 
	double sr_aux;

	for(int i = 0; i < n_pan; i++){
		sr_aux = (x_pan[i+1] - x_pan[i])*(x_pan[i+1] - x_pan[i]);
		sr_aux += (y_pan[i+1] - y_pan[i])*(y_pan[i+1] - y_pan[i]);
		R_pan[i] = pow(sr_aux,0.5);
		//cout << "Longitud panel "<< i <<" = " << R_pan[i] <<endl;
	}

	double nx[n_pan], ny[n_pan];

	for(int i = 0; i<n_pan; i++){
		nx[i] = -(y_pan[i] - y_pan[i+1])/R_pan[i];
		ny[i] =  (x_pan[i] - x_pan[i+1])/R_pan[i];
		//cout<< "n_"<<i<<" = ("<<nx[i] << "," << ny[i]<< ")" << endl;
	}

	/*############################################
	#### Fin definicion del cilindro##############
	#############################################*/

/*#############################################################################*/

	/*############################################
        #### Inicio Cuadratura Gaussiana##############
 	##############################################*/

	double t_cuad[n_cuad]; // nodos de gauss chebyshev
	
	for(int i = 0 ; i < n_cuad; i++){
		t_cuad[i] = cos(pi*(2.*(i+1)-1.)/(2.*n_cuad));
		//cout << "t_"<< i+1 << "="<< t_cuad[i] << endl;
	}


	/*############################################
        #### fin Cuadratura Gaussiana#################
 	##############################################*/

/*###########################################################################*/

	/*#############################################
	######## inicio source & target ###############
	###############################################*/	
	
	vector<Point> sources, target;

	/*Para rellenar el arreglo de target se requiere la informacion
	del numero de paneles y los puntos de colocacion*/

	/*x_col, y_col, n_pan*/

	for(int i = 0; i < n_pan; i++){
		target.push_back(Point(x_col[i],y_col[i]));
		//cout<<"target_"<<i<<" = ("<<target[i].x<<","<<target[i].y<<")"<<endl;
	}

	/*Para rellenar el arreglo sources se requiere la informacion
	del numero de paneles, posicion de los paneles, numero de nodos 
	de cuadratura y nodos de cuadratura */

	/*x_pan, y_pan, n_pan, t_cuad, n_cuad*/

	int N_sources = n_pan*n_cuad;
	double xg[n_cuad], yg[n_cuad];
	double aux[N_sources][2], w_cuad[N_sources];
	int j = 0;

	for(int p = 0; p < n_pan; p++){
		for(int i = 0; i < n_cuad; i++){
			xg[i] = 0.5*(x_pan[p+1]+x_pan[p]) + 0.5*(x_pan[p+1]-x_pan[p])*t_cuad[i];
			yg[i] = 0.5*(y_pan[p+1]+y_pan[p]) + 0.5*(y_pan[p+1]-y_pan[p])*t_cuad[i];
		}

		for(int i = n_cuad - 1; i >= 0; i--){
			aux[j][0] = xg[i];
			aux[j][1] = yg[i];
			w_cuad[j] = (pi/n_cuad)*pow(1.-pow(t_cuad[i],2.),0.5);
			j++;
		}
	}

	for(int i = 0; i < N_sources;i++){
		sources.push_back(Point(aux[i][0],aux[i][1]));
		//cout << "source_"<<i<<"= ("<< sources[i].x<<","<<sources[i].y<<")"<<endl;
		//cout << "w_cuad_"<<i<<"= "<<w_cuad[i]<<endl;		
	}
		

	/*#############################################
	######## fin source & target ###############
	###############################################*/

	/*#############################################
	########### velocidad #########################
	###############################################*/	

	double ux[n_pan], uy[n_pan];
	double bx[n_pan], by[n_pan];

	for(int i=0; i < n_pan; i++){
		ux[i] = u_0;
		uy[i] = 0.0;

		bx[i] = ux[i];
		by[i] = uy[i];
	}

	/*######### vel ###############################*/

/*##################################################################################*/

	clock_t startBuild	=	clock();

	/*#############################################
	######Inicio calculo streeslet#################
	###############################################*/

	
	
	double Txxx, Txxy, Tyxx, Tyxy;
	double Txyx, Txyy, Tyyx, Tyyy;
	double r_aux, dx, dy;

	for(int q = 0; q < n_pan; q++){
		for(int i=0; i < N_sources; i++){
			dx = sources[i].x - target[q].x;
			dy = sources[i].y - target[q].y;
			r_aux = pow(dx*dx + dy*dy, 0.5);

			Txxx = -4.*(dx*dx*dx)/pow(r_aux,4.);
			Txxy = -4.*(dx*dx*dy)/pow(r_aux,4.);
			
			Tyxx = Txxy;
			Tyxy = -4.*(dy*dx*dy)/pow(r_aux,4.);
			
			Txyx = Txxy;
			Txyy = Tyxy;

			Tyyx = Txyy;
			Tyyy = -4.*(dy*dy*dy)/pow(r_aux,4.);

			if(i/n_cuad != q){
				bx[q] += (-1./dpi)*(R_pan[i/n_cuad]/2.)*w_cuad[i]*(ux[i/n_cuad]*(Txxx*nx[i/n_cuad] 
					+ Txxy*ny[i/n_cuad]) + uy[i/n_cuad]*(Tyxx*nx[i/n_cuad] + Tyxy*ny[i/n_cuad]));

				by[q] += (-1./dpi)*(R_pan[i/n_cuad]/2.)*w_cuad[i]*(ux[i/n_cuad]*(Txyx*nx[i/n_cuad] 
					+ Txyy*ny[i/n_cuad]) + uy[i/n_cuad]*(Tyyx*nx[i/n_cuad] + Tyyy*ny[i/n_cuad]));
			}
		}
	}

	double bj[2*n_pan];

	for(int i=0; i < n_pan; i++){
		bj[i] = bx[i];
		bj[i+n_pan] = by[i];
	}

	clock_t endBuild	=	clock();

	double Timep2p = double(endBuild-startBuild)/double(CLOCKS_PER_SEC);

	cout << endl << "Tiempo p2p: " << Timep2p << endl <<endl;

	/*
	cout<<"b  = "<<endl;

	for(int i=0;i < 2*n_pan; i++){
		cout.precision(8);
		cout<< bj[i] << endl;
	}

	*/

	/*#############################################
	#########Fin calculo streeslet#################
	###############################################*/	

/*################################################################################################################################
##################################################################################################################################
##################################################################################################################################*/

	/*###############BBFMM2D###############################*/

	clock_t starfmm	=	clock();

	unsigned long N = sources.size() + target.size();
	unsigned m = 1;
	vector<Point> location;
	unsigned short nChebNodes = 6;

	for(int i=0; i < N; i++){
		if(i < sources.size()){
			location.push_back(Point(sources[i].x,sources[i].y));
		}
		else{
			location.push_back(Point(target[i-sources.size()].x,target[i-sources.size()].y));
		}
	}

/*
	for(int i=0; i < N; i++){
		cout<<i<<" = ("<<location[i].x<<","<<location[i].y<<")"<<endl;
	}
*/

	double* beta11;
	beta11 = new double[N*m];

	double* beta12;
	beta12 = new double[N*m];

	double* beta21;
	beta21 = new double[N*m];

	double* beta22;
	beta22 = new double[N*m];


	for(int i = 0; i < N_sources; i++){
		beta11[i] = (-1./ddpi)*R_pan[i/n_cuad]*w_cuad[i]*ux[i/n_cuad]*nx[i/n_cuad];
		//cout << "b11_"<<i<<" = "<<beta11[i]<< endl;
	}

	for(int i= N_sources; i < N; i++){
		beta11[i] = 0.0;
		//cout << "b11_"<<i<<" = "<<beta11[i]<< endl;
	}

	
	for(int i = 0; i < N_sources; i++){
		beta12[i] = (-1./ddpi)*R_pan[i/n_cuad]*w_cuad[i]*ux[i/n_cuad]*ny[i/n_cuad];
		//cout << "b12_"<<i<<" = "<<beta12[i]<< endl;
	}

	for(int i= N_sources; i < N; i++){
		beta12[i] = 0.0;
		//cout << "b12_"<<i<<" = "<<beta12[i]<< endl;
	}

	
	for(int i = 0; i < N_sources; i++){
		beta21[i] = (-1./ddpi)*R_pan[i/n_cuad]*w_cuad[i]*uy[i/n_cuad]*nx[i/n_cuad];
		//cout << "b21_"<<i<<" = "<<beta21[i]<< endl;
	}

	for(int i= N_sources; i < N; i++){
		beta21[i] = 0.0;
		//cout << "b21_"<<i<<" = "<<beta21[i]<< endl;
	}

	for(int i = 0; i < N_sources; i++){
		beta22[i] = (-1./ddpi)*R_pan[i/n_cuad]*w_cuad[i]*uy[i/n_cuad]*ny[i/n_cuad];
		//cout << "b22_"<<i<<" = "<<beta22[i]<< endl;
	}

	for(int i= N_sources; i < N; i++){
		beta22[i] = 0.0;
		//cout << "b22_"<<i<<" = "<<beta22[i]<< endl;
	}


	H2_2D_Tree Atree(nChebNodes, beta11, location, N, m);// Build the fmm tree;
	
	double* potentialbxxx;
	potentialbxxx = new double[N*m];

	double* potentialbxyx;	
	potentialbxyx = new double[N*m];

	myKernelTxxx bxxx;
    	bxxx.calculate_Potential(Atree,potentialbxxx);

	myKernelTxyx bxyx;
    	bxyx.calculate_Potential(Atree,potentialbxyx);

	


	H2_2D_Tree Btree(nChebNodes, beta12, location, N, m);// Build the fmm tree;
	
	double* potentialbxxy;
	potentialbxxy = new double[N*m];

	double* potentialbxyy;	
	potentialbxyy = new double[N*m];

	myKernelTxyx bxxy;
    	bxxy.calculate_Potential(Btree,potentialbxxy);

	myKernelTyxy bxyy;
    	bxyy.calculate_Potential(Btree,potentialbxyy);


	H2_2D_Tree Ctree(nChebNodes, beta21, location, N, m);// Build the fmm tree;
	
	double* potentialbyxx;
	potentialbyxx = new double[N*m];

	double* potentialbyyx;	
	potentialbyyx = new double[N*m];

	myKernelTxyx byxx;
    	byxx.calculate_Potential(Ctree,potentialbyxx);

	myKernelTyxy byyx;
    	byyx.calculate_Potential(Ctree,potentialbyyx);


	H2_2D_Tree Dtree(nChebNodes, beta22, location, N, m);// Build the fmm tree;
	
	double* potentialbyxy;
	potentialbyxy = new double[N*m];

	double* potentialbyyy;	
	potentialbyyy = new double[N*m];

	myKernelTyxy byxy;
    	byxy.calculate_Potential(Dtree,potentialbyxy);

	myKernelTyyy byyy;
    	byyy.calculate_Potential(Dtree,potentialbyyy);


/*####################calculo de bx, by#####################################################################################*/


	double bxfmm[n_pan], byfmm[n_pan];
	double bfmm[2*n_pan];

	for(int i = N_sources; i < N; i++){
		bxfmm[i-N_sources] = ux[i] + potentialbxxx[i] + potentialbxxy[i] + potentialbyxx[i] + potentialbyxy[i] ;
		byfmm[i-N_sources] = uy[i] + potentialbxyx[i] + potentialbxyy[i] + potentialbyyx[i] + potentialbyyy[i] ;
		bfmm[i-N_sources] = bxfmm[i-N_sources];
		bfmm[i+ n_pan -N_sources] = byfmm[i-N_sources]; 		
	}

	/*
	cout << "bfmm = "<<endl;
	
	for(int i = 0; i < 2*n_pan; i++){
		cout<< bfmm[i]<<endl;
	}

	*/

	clock_t endfmm = clock();

	double Timefmm = double(endfmm-starfmm)/double(CLOCKS_PER_SEC);

	cout << endl << "Tiempo bbfmm: " << Timefmm << endl <<endl;

	/*
	for(int i = 0; i < N; i++){
		cout << potentialbxxx[i]<< endl;
	}
	*/
	
	return 0;
}


/*
double Cil_make(int n, double R, double h){

}
*/


