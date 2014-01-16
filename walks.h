

#include <fstream>
#include <iostream>
#include <stdio.h>
#include <vector>
#include <string>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_complex.h>	    
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_block.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>




using namespace std; 


// helpful constants to have down
gsl_complex I = gsl_complex_rect(0,1);
gsl_complex I_ = gsl_complex_rect(0,-1);
gsl_complex one = gsl_complex_rect(1,0);
gsl_complex zero = gsl_complex_rect(0,0);
gsl_complex pi_2 = gsl_complex_rect(3.14159/2.0,0);

// nice functions to have
void my_vector_complex_conjugate(gsl_vector_complex * v, int size);
void my_vector_complex_abs(gsl_vector_complex * input, gsl_vector * output, int size);
// print functions
void print_matrix_complex(gsl_matrix_complex * k, int size);
void print_vector_complex(gsl_vector_complex * v, int size);
void print_vector(gsl_vector * v, int size);
void print_complex(gsl_complex c);
void print(float c);

// def the quantum walk class

class q_walk {
	public:
	// length of n by n matrix
	int size;
	// ajacency matrix
	gsl_matrix_complex * aj_matrix;
	// start vector
	gsl_vector_complex * start;
	// spectral decomposition vector // why a vector? I just dont want to deal with arrays
	vector<gsl_matrix_complex*> E_r;
	// eigen values in same order as spec matrices
	gsl_vector * lamda;
	// current walkey
	gsl_vector_complex * current; 
	gsl_vector * current_abs;
	// store max on each edge
	gsl_vector * max_store;
	gsl_vector * max_store_time;
	// what nodes have PGST (with in .995)
	vector<bool> PGST;
	bool UPGST;
	bool any_PGST;
	

	// THIS IS JUST USED FOR THE GRAPHICS STUFF. TAKE OUT IF WANTED
	vector<double> x;
	vector<double> y;
	vector<double> z; 

	

	// constroctors read in the matrix	(uses cin haha sorry) 
	q_walk(string file);
	// constroctor that takes in aj matrix
	q_walk(gsl_matrix_complex * k, int n);
	// calc eigen values and spectral decomposition
	void eigen_stuff();
	// find distrobution for time t
	void current_dis(double time);
	// calc max of each node
	void max(double time);
	// calc max starting at node 0
	void max_special(double time);
	// say if the graph gets close to a node
	void PGST_test(double how_near);
	// spit out the nodes with PGST
	void PGST_print();
	// spit out if they all have PGST
	void UPGST_print();

};


q_walk::q_walk(string file)
{
	// create file stream
	ifstream reading;
	reading.open(file.c_str());
	
	// helpful later
	double real;
	double imag;
	gsl_complex com;
	double x_pos;
	double y_pos;
	double z_pos;

	// read in the first few things. will extend
	int n = 0;
	reading >> n;
	size = n;

	// allocate everything
	aj_matrix = gsl_matrix_complex_alloc(n,n);
	gsl_matrix_complex_set_zero(aj_matrix);
	lamda = gsl_vector_alloc(n);
	gsl_vector_set_zero(lamda);
	start = gsl_vector_complex_alloc(n);
	current = gsl_vector_complex_alloc(n);
	current_abs = gsl_vector_alloc(n);
	gsl_vector_complex_set_zero(start);
	gsl_vector_complex_set_zero(current);
	gsl_vector_set_zero(current_abs);
	max_store = gsl_vector_alloc(n);
	max_store_time = gsl_vector_alloc(n);
	gsl_vector_set_zero(max_store);
	gsl_vector_set_zero(max_store_time);

	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < n; j++)
		{
			reading >> real >> imag;	
			com = gsl_complex_rect(real,imag);
			gsl_matrix_complex_set(aj_matrix,i,j,com);
		}
		// oh what a terrible terrible way to read in a file...
	}
	for(int i = 0; i < n; i++)
	{	
		reading >> x_pos >> y_pos >> z_pos;
		x.push_back(x_pos);
		y.push_back(y_pos);
		z.push_back(z_pos);
	}
	for(int i = 0; i < n; i++)
	{
		reading >> real >> imag;	
		com = gsl_complex_rect(real,imag);
		gsl_vector_complex_set(start,i,com);
		gsl_vector_set(max_store, i, real);
		gsl_vector_set(max_store_time, i, 0);
	}
	reading.close();
	//print_matrix_complex(aj_matrix);
	eigen_stuff();
	current_dis(0);	
}

q_walk::q_walk(gsl_matrix_complex * k, int n)
{
	size = n;
	aj_matrix = gsl_matrix_complex_alloc(n,n);
	gsl_matrix_complex_set_zero(aj_matrix);
	lamda = gsl_vector_alloc(n);
	gsl_vector_set_zero(lamda);
	start = gsl_vector_complex_alloc(n);
	current = gsl_vector_complex_alloc(n);
	current_abs = gsl_vector_alloc(n);
	gsl_vector_complex_set_zero(start);
	gsl_vector_complex_set_zero(current);
	gsl_vector_set_zero(current_abs);
	max_store = gsl_vector_alloc(n);
	max_store_time = gsl_vector_alloc(n);
	gsl_vector_set_zero(max_store);
	gsl_vector_set_zero(max_store_time);
	

	for( int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			gsl_matrix_complex_set(aj_matrix, i, j, gsl_matrix_complex_get(k, i, j));
		}
	}
	
	gsl_vector_complex_set(start, 0, gsl_complex_rect(1,0));

	eigen_stuff();
	current_dis(0);	


}




void q_walk::eigen_stuff()
{
	// stuff I need. remmember to free it at the end
	vector<gsl_vector_complex*> e_i;
	gsl_matrix_complex * MV = gsl_matrix_complex_alloc(size,size);
	gsl_eigen_hermv_workspace * WV = gsl_eigen_hermv_alloc(size);
	vector<gsl_matrix_complex*> E_i_a;	
	vector<gsl_matrix_complex*> E_i_b;	

	for(int i = 0; i<size; i++)
	{	
		E_i_a.push_back(gsl_matrix_complex_alloc(size,size));
		E_i_b.push_back(gsl_matrix_complex_alloc(size,size));
		gsl_matrix_complex_set_zero(E_i_a[i]);
		gsl_matrix_complex_set_zero(E_i_b[i]);

	}

	for(int i = 0; i<size; i++)
	{	
		E_r.push_back(gsl_matrix_complex_alloc(size,size));
		gsl_matrix_complex_set_zero(E_r[i]);
	}

	for(int i = 0; i<size; i++)
	{
		e_i.push_back(gsl_vector_complex_alloc(size));
		gsl_vector_complex_set_zero(e_i[i]);
	}	

	//print_matrix_complex(aj_matrix,size);
	
	gsl_eigen_hermv(aj_matrix,lamda,MV,WV);
	gsl_eigen_hermv_free(WV);

	//print_matrix_complex(MV,size);
	//print_vector(lamda,size);

	for(int i = 0; i<size; i++)
	{
		gsl_matrix_complex_get_col(e_i[i],MV,i);
	}

	for(int i = 0; i<size; i++)
	{
		gsl_matrix_complex_set_col(E_i_a[i],0,e_i[i]);
	}

	for(int i = 0; i<size; i++)
	{
		my_vector_complex_conjugate(e_i[i], size);
	}	

	for(int i = 0; i<size; i++)
	{
		gsl_matrix_complex_set_row(E_i_b[i],0,e_i[i]);
	}

	for(int i = 0; i<size; i++)
	{
		gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, one, E_i_a[i], E_i_b[i], zero, E_r[i]);
	}
	// test the E_r by seeing if E_r=~E_r^2 (property of things)
	for(int i = 0; i<size; i++)
	{
//		print_matrix_complex(E_r[i],size);
	//	gsl_matrix_complex * testery = gsl_matrix_complex_alloc(size,size);
	//	gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, one, E_r[i], E_r[i], zero, testery);
	//	print_matrix_complex(testery,size); 
	}
	// free everything!
	for (int i = 0; i<size; i++)
	{
		gsl_matrix_complex_free(E_i_a[i]);
		gsl_matrix_complex_free(E_i_b[i]);
		gsl_vector_complex_free(e_i[i]);
	}
	gsl_matrix_complex_free(MV);
}

void q_walk::current_dis(double time)
{

	vector<gsl_complex> F_i;
	vector<gsl_complex> V_i;
	gsl_complex R;
	
	vector<gsl_matrix_complex*> work_s;	
	gsl_matrix_complex * All = gsl_matrix_complex_alloc(size,size);

	gsl_matrix_complex_set_zero(All);
	gsl_complex t;
	GSL_SET_COMPLEX(&t, time, 0); 

	for (int i = 0; i < size; i++)
	{
		GSL_SET_COMPLEX(&R, gsl_vector_get(lamda,i),0);
		F_i.push_back(R);
		work_s.push_back(gsl_matrix_complex_alloc(size,size));
		gsl_matrix_complex_memcpy(work_s[i],E_r[i]);
	}

//	print_vector(lamda,size);
//	print_complex(F_i[0]);
//	print_complex(F_i[1]);
	
	for (int i = 0; i < size; i++)
	{
		V_i.push_back(gsl_complex_exp(gsl_complex_mul(gsl_complex_mul(I_,F_i[i]),t)));
	}

//	print_complex(V_i[0]);
//	print_complex(V_i[1]);

	for (int i = 0; i < size; i++)
	{
		gsl_matrix_complex_scale(work_s[i],V_i[i]);	
	}
/*
	for (int i = 0; i < size; i++)
	{
		print_matrix_complex(E_r[i],size);
	}
*/
	//print_matrix_complex(All,size);
	for (int i = 0; i < size; i++)
	{
		gsl_matrix_complex_add(All, work_s[i]);
	}

//	print_vector_complex(start, size);
	//print_matrix_complex(All, size);
	//print_vector_complex(current,size);

	gsl_blas_zgemv (CblasNoTrans, one, All, start, zero, current);





	// free stuff (R, F_i, t, all, V_i)	
	//gsl_complex_free(R);
	//gsl_complex_free(&t);
	gsl_matrix_complex_free(All);
	F_i.clear();
	V_i.clear();
	for (int i = 0; i < size; i++)
	{
		gsl_matrix_complex_free(work_s[i]);
	}
	work_s.clear();
	
	
	my_vector_complex_abs(current, current_abs, size);
	max(time);
	//print_vector(current_abs, size);
	//print_vector(max_store, size);
	//print_vector(max_store_time, size);
}

void q_walk::max(double time)
{

	for (int i = 0; i < size; i++)
	{
		if(gsl_vector_get(current_abs, i) > gsl_vector_get(max_store, i))
		{
			gsl_vector_set(max_store, i, gsl_vector_get(current_abs, i));
			gsl_vector_set(max_store_time, i, time);
		} 
	}

}

void q_walk::PGST_test(double how_near)
{
	any_PGST = false;
	for (int i = 0; i < size; i++)
	{
		if (gsl_vector_get(max_store, i) > (1-how_near))
		{
			PGST.push_back(true); 
		}
		else
		{
			PGST.push_back(false);
		}
	}

	// test
	for (int i = 1; i < size; i++)
	{
		if (gsl_vector_get(max_store, i) > (1-how_near))
		{
			any_PGST = true;
		}
	}	

	
	UPGST = true;
	for (int i = 0; i < size; i++)
	{
		if(!PGST[i])	
		{
			UPGST = false;
		}
	}


}

// To check stuff
void q_walk::PGST_print()
{
	cout << "this graph has PGST at node \n";
	for (int i = 0; i < size; i++)
	{
		if (PGST[i])
		{
			cout << i << "\n";
		}
	}

}

void q_walk::UPGST_print()
{
	if (UPGST)
	{
		cout << "this graph does \n";
	}
	if (!UPGST)
	{
		cout << "nope. :( \n";
	}
	
}
	




// ----------------------------------------------------------------------
// Helpful  functions --------------------------------------------------
// ----------------------------------------------------------------------
void my_vector_complex_conjugate(gsl_vector_complex * v, int size)
{
	for (int i = 0; i < size; i++)
	{
		gsl_vector_complex_set(v,i,gsl_complex_conjugate(gsl_vector_complex_get(v,i)));
	}

}

void print_vector_complex(gsl_vector_complex * v, int size)
{

	for (int i = 0; i < size; i++)  
           	printf ("VC(%d) = %g, %g i\n", i, 
                GSL_REAL(gsl_vector_complex_get (v, i)),
		GSL_IMAG(gsl_vector_complex_get (v, i)));	
	
	printf ("\n");
}

void print_vector(gsl_vector * v, int size)
{

	for (int i = 0; i < size; i++)  
           	printf ("V(%d) = %g, \n", i, 
                gsl_vector_get (v, i));	

	printf ("\n");
}

void print_matrix_complex(gsl_matrix_complex * m, int size)
{
	
	for (int i = 0; i < size; i++) 
	{
        	for (int j = 0; j < size; j++)
		{
           		printf ("CM %g, %g i    ", 
                	GSL_REAL(gsl_matrix_complex_get (m, i, j)),
			GSL_IMAG(gsl_matrix_complex_get (m, i, j)));
		}
	printf ("\n");	
	}

	printf ("\n");
}
void print_complex(gsl_complex c)
{
	printf ("%g, %g i  \n",GSL_REAL(c),GSL_IMAG(c));
}

void print(float c)
{
	printf ("%g \n",c);
}

void my_vector_complex_abs(gsl_vector_complex * input, gsl_vector * output, int size)
{
	for (int i = 0; i < size; i++)
	{
		gsl_vector_set(output, i, (gsl_complex_abs(gsl_vector_complex_get(input, i))*gsl_complex_abs(gsl_vector_complex_get(input, i))));
	}
}


// for UPGST class
void print_matrix_complex_to_file(ofstream & myfile, gsl_matrix_complex * m, int size)
{
	for (int i = 0; i < size; i++) 
	{
        	for (int j = 0; j < size; j++)
		{
                	myfile << GSL_REAL(gsl_matrix_complex_get (m, i, j)) << " " << GSL_IMAG(gsl_matrix_complex_get (m, i, j)) << " ";
		}
	myfile << "\n \n";	
	}

	myfile << "\n \n \n \n ";
	


}

// for printing and then displaying it graphicly
void print_matrix_complex_to_file_graphics(ofstream & myfile, gsl_matrix_complex * m, int size)
{

	gsl_complex circle;
	
	myfile << size << " \n \n \n";	

	for (int i = 0; i < size; i++) 
	{
        	for (int j = 0; j < size; j++)
		{
                	myfile << GSL_REAL(gsl_matrix_complex_get (m, i, j)) << " " << GSL_IMAG(gsl_matrix_complex_get (m, i, j)) << " ";
		}
	myfile << "\n \n";	
	}

	myfile << "\n \n";
	
	for (int i = 0; i < size; i++)
	{
		circle = gsl_complex_polar(1, 2*3.14*i*(1.0/size));
		myfile << GSL_REAL(circle) << " " << GSL_IMAG(circle) << " " << "0" << "\n";
	}

	myfile << "\n \n ";
	
	myfile << "1 0 ";
	for (int i = 0; i < size-1; i++)
	{
		myfile << "0 0 "; 
	}
	


}

// for printing to a dot file and then graphing it
void print_matrix_complex_to_file_dot(ofstream & myfile, gsl_matrix_complex * m, int size)
{

	myfile << "digraph G { \n";
	
	for (int i = 0; i < size; i++) 
	{
        	for (int j = i; j < size; j++)
		{
			if (GSL_IMAG(gsl_matrix_complex_get (m, i, j)) == 1)
			{
				myfile << i << " -> " << j;
			}
			if (GSL_IMAG(gsl_matrix_complex_get (m, i, j)) == -1)
			{
				myfile << j << " -> " << i;
			}
			myfile << "\n";
		}
	}
	
	myfile << "}";



}

// for printing to a latex file and then graphing it
void print_matrix_complex_to_file_latex(ofstream & myfile, gsl_matrix_complex * m, gsl_vector * max, gsl_vector * eigen, double time_interval, double dt, int size)
{

	gsl_complex circle;
	
	char slash = 92;
	myfile << slash << "documentclass[a4paper, 12pt]{article} \n";
	myfile << slash << "usepackage{amsmath} \n";
	myfile << slash << "usepackage{tikz} \n";
	myfile << slash << "usepackage{amsthm} \n";
	myfile << slash << "usetikzlibrary{arrows, snakes,backgrounds} \n";
	myfile << slash << "begin{document} \n";
	myfile << slash << "section{" << size << " node graph} \n" ;
	
	myfile << slash << "tikzstyle{haha} = [circle, draw=blue!50,fill=blue!20,thick] \n";
	myfile << slash << "begin{center} \n";
	myfile << slash << "begin{tikzpicture}[->,>=stealth',shorten >=1pt,auto,node distance=10.0.0cm, semithick] \n";
	
	for (int i = 0; i < size; i++)
	{
	
		circle = gsl_complex_polar(size/1.5, 2*3.14*i*(1.0/size));
		myfile << slash << "node at (" << GSL_REAL(circle) << "," << GSL_IMAG(circle) << ") [haha] (" << i << ") {" << i << "}; \n";
	}

	myfile << slash << "path   ";	

	for (int i = 0; i < size; i++) 
	{
        	for (int j = i; j < size; j++)
		{
			if (GSL_IMAG(gsl_matrix_complex_get (m, i, j)) == 1)
			{
				myfile  << "	(" << i << ") edge node {} (" << j << ") \n";
			}
			if (GSL_IMAG(gsl_matrix_complex_get (m, i, j)) == -1)
			{
				myfile  << "	(" << j << ") edge node {} (" << i << ") \n";
			}
			if (GSL_REAL(gsl_matrix_complex_get (m, i, j)) == 1)
			{
				myfile  << "	(" << j << ") edge[-] node {1} (" << i << ") \n";
			}
			if (GSL_REAL(gsl_matrix_complex_get (m, i, j)) == -1)
			{
				myfile  << "	(" << j << ") edge[-] node {-1} (" << i << ") \n";
			}
			myfile << "\n";
		}
		
	}
	myfile << ";";	

	myfile << slash << "end{tikzpicture} \n";
	myfile << slash << "end{center} \n";

	myfile << "here is the ajacency matrix for the graph \n";
	
	myfile << slash << "begin{equation} \n";
	myfile << slash << "begin{bmatrix} \n";
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			myfile << GSL_REAL(gsl_matrix_complex_get (m, i, j)) << " + " << GSL_IMAG(gsl_matrix_complex_get (m, i, j)) << "i ";
			if (j == size-1)
			{
				myfile << slash << slash;
			}
			else
			{
				myfile << "& ";
			}
		}
		myfile << "\n";
	}
	myfile << slash << "end{bmatrix} \n";
	myfile << slash << "end{equation} \n";
	myfile << "This graph's max at each node during the time interval (0," << time_interval << ") and dt = " << dt << "  is \n";

	
	myfile << slash << "begin{equation} \n";
	myfile << slash << "begin{bmatrix} \n";
	for (int i = 0; i < size; i++)
	{
			myfile << gsl_vector_get (max, i) << " " << slash << slash;
	}
	myfile << slash << "end{bmatrix} \n";
	myfile << slash << "end{equation} \n";

	myfile << "The eigen values of this graph are \n";
	myfile << slash << "begin{equation} \n";
	myfile << slash << "begin{bmatrix} \n";
	for (int i = 0; i < size; i++)
	{
			myfile << gsl_vector_get (eigen, i) << " " << slash << slash;
	}
	myfile << slash << "end{bmatrix} \n";
	myfile << slash << "end{equation} \n";
	

	myfile << slash << "end{document}";	



}
