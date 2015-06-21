/*
===============================================================================
   Implementation of TM-align in C/C++   

   This program is written by Jianyi Yang at
   Yang Zhang lab
   Center for Computational Medicine and Bioinformatics 
   University of Michigan 
   100 Washtenaw Avenue, Ann Arbor, MI 48109-2218 
                                                       
           
   Please report bugs and questions to yangji@umich.edu or zhng@umich.edu
===============================================================================
*/
#define MAXLEN 10000                        //maximum length of filenames
char version[20];                          //version 
 
 
//global variables
double D0_MIN;                             //for d0
double Lnorm;                              //normalization length
double score_d8, d0, d0_search, dcu0;      //for TMscore search
double **score;            			       //Input score table for dynamic programming
bool   **path;                             //for dynamic programming  
double **val;                              //for dynamic programming  
int    xlen, ylen, minlen;                 //length of proteins
double **xa, **ya;                         //for input vectors xa[0...xlen-1][0..2], ya[0...ylen-1][0..2]
                                           //in general, ya is regarded as native structure --> superpose xa onto ya
int    *xresno, *yresno;                   //residue numbers, used in fragment gapless threading 
double **xtm, **ytm;                       //for TMscore search engine
double **xt;                               //for saving the superposed version of r_1 or xtm
char   *seqx, *seqy;                       //for the protein sequence 
int    *secx, *secy;                       //for the secondary structure 
double **r1, **r2;                         //for Kabsch rotation 
double t[3], u[3][3];                      //Kabsch translation vector and rotation matrix



//argument variables
char out_reg[MAXLEN];
double Lnorm_ass, Lnorm_d0, d0_scale, d0A, d0B, d0u, d0a;
bool o_opt, a_opt, u_opt, d_opt, v_opt;
double TM3, TM4, TM5;


#include "basic_fun.h"
#include "NW.h"
#include "Kabsch.h"
#include "TMalign.h"


void print_help(char *arg)
{

	cout <<endl;
	cout << " *****************************************************************************" << endl
		 << " * TM-align (Version "<< version <<"): A protein structural alignment algorithm     *" << endl
		 << " * Reference: Y Zhang and J Skolnick, Nucl Acids Res 33, 2302-9 (2005)       *" << endl
		 << " * Please email your comments and suggestions to Yang Zhang (zhng@umich.edu) *" << endl
		 << " *****************************************************************************" << endl;	
	cout << endl
		 << " Usage: " << arg << " [Options]" << endl << endl
		 << " Options:" << endl
		 << "       -A    input filename of structure A, PDB format, required" << endl << endl
		 << "       -B    input filename of structure B, PDB format, required" << endl
		 << "             structure B will be superimposed onto structure A" <<endl << endl
		 << "       -u    TM-score normalized by user assigned length" << endl
		 << "             warning: it should be >= minimum length of the two structures" << endl
		 << "             otherwise, TM-score may be >1" << endl << endl
		 << "       -a    TM-score normalized by the average length of two structures" << endl 
		 << "             T or F, (default F)" << endl << endl
		 << "       -d    TM-score scaled by an assigned d0, e.g. 5 Angstroms" << endl << endl
		 << "       -o    output filename for superimposed structure of B" << endl
		 << "             To view the superimposed structures of aligned regions by rasmol:" << endl
		 << "             >rasmol -script filename" << endl
		 << "             To view the superimposed structures of all regions by rasmol:" << endl
		 << "             >rasmol -script filename_all" << endl << endl
		 << "       -v    print the version of TM-align" << endl << endl
		 << "       -h    print this help" << endl << endl
		 << "(Options -u, -a, -d -o won't change the final structure alignment)" << endl << endl
		 << " Example usages:" << endl
		 << " "<< arg <<" -A PDB2.pdb -B PDB1.pdb" << endl
		 << " "<< arg <<" -A PDB2.pdb -B PDB1.pdb -u 100 -d 5.0" << endl
		 << " "<< arg <<" -A PDB2.pdb -B PDB1.pdb -a T -o PDB1.sup" << endl << endl;
       
  exit(EXIT_SUCCESS);

}



void parameter_set4search(int xlen, int ylen)
{
	//parameter initilization for searching: D0_MIN, Lnorm, d0, d0_search, score_d8
	D0_MIN=0.5; 
	dcu0=4.25;                       //update 3.85-->4.25
 
	Lnorm=getmin(xlen, ylen);        //normaliz TMscore by this in searching
    if(Lnorm<=19)                    //update 15-->19
    {
        d0=0.168;                   //update 0.5-->0.168
    }
    else
    {
        d0=(1.24*pow((Lnorm*1.0-15), 1.0/3)-1.8);
    }
	D0_MIN=d0+0.8;              //this should be moved to above
    d0=D0_MIN;                  //update: best for search    


	d0_search=d0;	
	if(d0_search>8) d0_search=8;
	if(d0_search<4.5) d0_search=4.5;


    score_d8=1.5*pow(Lnorm*1.0, 0.3)+3.5; //remove pairs with dis>d8 during search & final
}

void parameter_set4final(double len)
{
	D0_MIN=0.5; 
 
	Lnorm=len;            //normaliz TMscore by this in searching
    if(Lnorm<=21)         
    {
        d0=0.5;          
    }
    else
    {
        d0=(1.24*pow((Lnorm*1.0-15), 1.0/3)-1.8);
    }
    if(d0<D0_MIN) d0=D0_MIN;   

	d0_search=d0;	
	if(d0_search>8) d0_search=8;
	if(d0_search<4.5) d0_search=4.5;  

}


void parameter_set4scale(int len, double d_s)
{
 
	d0=d_s;          
	Lnorm=len;            //normaliz TMscore by this in searching

	d0_search=d0;	
	if(d0_search>8) d0_search=8;
	if(d0_search<4.5) d0_search=4.5;  

}




int main(int argc, char *argv[])
{
	strcpy(version, "20120126");

    if (argc < 2) 
    {
        print_help(argv[0]);        
    } 
	
	
	clock_t t1, t2;
	t1 = clock();

    /*********************************************************************************/
	/*                                get argument                                   */ 
    /*********************************************************************************/
    char xname[MAXLEN], yname[MAXLEN],  Lnorm_ave[MAXLEN];
	bool A_opt, B_opt, h_opt=false;
	A_opt = B_opt = o_opt = a_opt = u_opt = d_opt = v_opt = false;

	for(int i = 0; i < argc; i++)
	{
		if ( !strcmp(argv[i],"-B") && i < argc ) { strcpy(xname,    argv[i+1]);      B_opt = true; }
		if ( !strcmp(argv[i],"-A") && i < argc ) { strcpy(yname,    argv[i+1]);      A_opt = true; }
		if ( !strcmp(argv[i],"-o") && i < argc ) { strcpy(out_reg,  argv[i+1]);      o_opt = true; }
		if ( !strcmp(argv[i],"-u") && i < argc ) { Lnorm_ass      = atof(argv[i+1]); u_opt = true; }
		if ( !strcmp(argv[i],"-a") && i < argc ) { strcpy(Lnorm_ave, argv[i+1]);     a_opt = true; }
		if ( !strcmp(argv[i],"-d") && i < argc ) { d0_scale        = atof(argv[i+1]); d_opt = true; }
		if ( !strcmp(argv[i],"-v") ) { v_opt = true; }
		if ( !strcmp(argv[i],"-h") ) { h_opt = true; }
	}


	if(!B_opt || !A_opt)
	{

		if( h_opt )
		{
			print_help(argv[0]);    			
		}
		
		if(v_opt)
		{
			cout <<endl;
			cout << " *****************************************************************************" << endl
				 << " * TM-align (Version "<< version <<"): A protein structural alignment algorithm     *" << endl
				 << " * Reference: Y Zhang and J Skolnick, Nucl Acids Res 33, 2302-9 (2005)       *" << endl
				 << " * Please email your comments and suggestions to Yang Zhang (zhng@umich.edu) *" << endl
				 << " *****************************************************************************" << endl;	
			exit(EXIT_FAILURE);
		}
	}

	if( !B_opt )
	{
		cout << "Please provide structure B" << endl;
		exit(EXIT_FAILURE);
	}		
	if( !A_opt )
	{
		cout << "Please provide structure A" << endl;
		exit(EXIT_FAILURE);
	}


	if( a_opt )
	{
		if(!strcmp(Lnorm_ave, "T"))
		{
		}
		else if(!strcmp(Lnorm_ave, "F"))
		{
			a_opt=false;
		}
		else
		{
			cout << "Wrong value for option -a!  It should be T or F" << endl;
			exit(EXIT_FAILURE);
		}
	}
	if( u_opt )
	{
		if(Lnorm_ass<=0)
		{
			cout << "Wrong value for option -u!  It should be >0" << endl;
			exit(EXIT_FAILURE);
		}
	}
	if( d_opt )
	{
		if(d0_scale<=0)
		{
			cout << "Wrong value for option -d!  It should be >0" << endl;
			exit(EXIT_FAILURE);
		}
	}














    /*********************************************************************************/
	/*                                load data                                      */ 
    /*********************************************************************************/
    load_PDB_allocate_memory(xname, yname);



    

    /*********************************************************************************/
	/*                                parameter set                                  */ 
    /*********************************************************************************/
	parameter_set4search(xlen, ylen);          //please set parameters in the function
    int simplify_step     = 40;               //for similified search engine
    int score_sum_method  = 8;                //for scoring method, whether only sum over pairs with dis<score_d8
        
	int i;
    int *invmap0          = new int[ylen+1]; 
    int *invmap           = new int[ylen+1]; 
    double TM, TMmax=-1;
	for(i=0; i<ylen; i++)
	{
		invmap0[i]=-1;
	}	


	double ddcc=0.4;
	if(Lnorm <= 40) ddcc=0.1;   //Lnorm was setted in parameter_set4search
      

    /*********************************************************************************/
	/*         get initial alignment with gapless threading                          */ 
    /*********************************************************************************/
    get_initial(xa, ya, xlen, ylen, invmap0);
    //find the max TMscore for this initial alignment with the simplified search_engin
    TM=detailed_search(xa, ya, xlen, ylen, invmap0, t, u, simplify_step, score_sum_method);
	if(TM>TMmax)
    {
        TMmax=TM;
    }           
    //run dynamic programing iteratively to find the best alignment
    TM=DP_iter(xa, ya, xlen, ylen, t, u, invmap, 0, 2, 30);
    if(TM>TMmax)
    {        
        TMmax=TM;
        for(int i=0; i<ylen; i++)
        {
            invmap0[i]=invmap[i];
        }
    }
	
    
	


    /*********************************************************************************/
	/*         get initial alignment based on secondary structure                    */ 
    /*********************************************************************************/	
	get_initial_ss(xa, ya, xlen, ylen, invmap);
    TM=detailed_search(xa, ya, xlen, ylen, invmap, t, u, simplify_step, score_sum_method);
    if(TM>TMmax)
    {
        TMmax=TM;
        for(int i=0; i<ylen; i++)
        {
            invmap0[i]=invmap[i];
        }
    } 
    if(TM > TMmax*0.2)
    {
        TM=DP_iter(xa, ya, xlen, ylen, t, u, invmap, 0, 2, 30);
        if(TM>TMmax)
        {
            TMmax=TM;
            for(int i=0; i<ylen; i++)
            {
                invmap0[i]=invmap[i];
            }
        }   
    }
	//output_align(invmap0, ylen);


	
    /*********************************************************************************/
	/*         get initial alignment based on local superposition                    */ 
    /*********************************************************************************/	
	//=initial5 in original TM-align
    if(get_initial_local(xa, ya, xlen, ylen, invmap))
	{		
		TM=detailed_search(xa, ya, xlen, ylen, invmap, t, u, simplify_step, score_sum_method);
		if(TM>TMmax)
		{
			TMmax=TM;
			for(int i=0; i<ylen; i++)
			{
				invmap0[i]=invmap[i];
			}
		} 
		if(TM > TMmax*ddcc)
		{
			TM=DP_iter(xa, ya, xlen, ylen, t, u, invmap, 0, 2, 2);
			if(TM>TMmax)
			{
				TMmax=TM;
				for(int i=0; i<ylen; i++)
				{
					invmap0[i]=invmap[i];
				}
			}   
		}
		//output_align(invmap0, ylen);
	}    
	else
	{
		cout << endl << endl << "Warning: initial alignment from local superposition fail!" << endl << endl <<endl;		
	}
	
	



    /*********************************************************************************/
	/*    get initial alignment based on previous alignment+secondary structure      */ 
    /*********************************************************************************/	
	//=initial3 in original TM-align
    get_initial_ssplus(xa, ya, xlen, ylen, invmap0, invmap);
    TM=detailed_search(xa, ya, xlen, ylen, invmap, t, u, simplify_step, score_sum_method);
    if(TM>TMmax)
    {
        TMmax=TM;
        for(i=0; i<ylen; i++)
        {
            invmap0[i]=invmap[i];
        }
    } 
    if(TM > TMmax*ddcc)
    {
        TM=DP_iter(xa, ya, xlen, ylen, t, u, invmap, 0, 2, 30);
        if(TM>TMmax)
        {
            TMmax=TM;
            for(i=0; i<ylen; i++)
            {
                invmap0[i]=invmap[i];
            }
        }   
    }
	//output_align(invmap0, ylen);   
      
	


	

    /*********************************************************************************/
	/*        get initial alignment based on fragment gapless threading              */ 
    /*********************************************************************************/	    
	//=initial4 in original TM-align
	get_initial_fgt(xa, ya, xlen, ylen, xresno, yresno, invmap);
    TM=detailed_search(xa, ya, xlen, ylen, invmap, t, u, simplify_step, score_sum_method);	
    if(TM>TMmax)
    {
        TMmax=TM;
        for(i=0; i<ylen; i++)
        {
            invmap0[i]=invmap[i];
        }
    } 
    if(TM > TMmax*ddcc)
    {
        TM=DP_iter(xa, ya, xlen, ylen, t, u, invmap, 1, 2, 2);
        if(TM>TMmax)
        {         
            TMmax=TM;
            for(i=0; i<ylen; i++)
            {
                invmap0[i]=invmap[i];
            }
        }   
    } 
	//output_align(invmap0, ylen); 






	
    //*********************************************************************************//
    //     The alignment will not be changed any more in the following                 //
    //*********************************************************************************//
	//check if the initial alignment is generated approately	
	bool flag=false;
	for(i=0; i<ylen; i++)
	{
		if(invmap0[i]>=0)
		{
			flag=true;
			break;			
		}			
	}		
	if(!flag) 
	{
		cout << "There is no alignment between the two proteins!" << endl;
		cout << "Program stop with no result!" << endl;
		return 1;
	}
	//cout << "final alignment" << endl;
	//output_align(invmap0, ylen);










    //*********************************************************************************//
    //       Detailed TMscore search engine  --> prepare for final TMscore             //
    //*********************************************************************************//       
    //run detailed TMscore search engine for the best alignment, and 
	//extract the best rotation matrix (t, u) for the best alginment
    simplify_step=1;
    score_sum_method=8;
    TM=detailed_search(xa, ya, xlen, ylen, invmap0, t, u, simplify_step, score_sum_method);

	//select pairs with dis<d8 for final TMscore computation and output alignment
	int n_ali8, k=0;
	int n_ali=0;
	int *m1, *m2;
	double d;
	m1=new int[xlen]; //alignd index in x
	m2=new int[ylen]; //alignd index in y
	do_rotation(xa, xt, xlen, t, u);
	k=0;
    for(int j=0; j<ylen; j++)
    {
        i=invmap0[j];
        if(i>=0)//aligned
        {
			n_ali++;        
            d=sqrt(dist(&xt[i][0], &ya[j][0]));
			if(d <= score_d8)
			{
				m1[k]=i;
				m2[k]=j;

				xtm[k][0]=xa[i][0];
                xtm[k][1]=xa[i][1];
                xtm[k][2]=xa[i][2];
                    
                ytm[k][0]=ya[j][0];
                ytm[k][1]=ya[j][1];
                ytm[k][2]=ya[j][2];								
				
				k++;
			}
		}
	}
	n_ali8=k;







    //*********************************************************************************//
    //                               Final TMscore                                     //
    //                     Please set parameters for output                            //
    //*********************************************************************************//
    double rmsd, TM1, TM2;
	double d0_out=5.0;  
    simplify_step=1;
    score_sum_method=0;

	double t0[3], u0[3][3];
	double d0_0, TM_0;
	double Lnorm_0=ylen;
	
	
	//normalized by length of structure A
	parameter_set4final(Lnorm_0);
	d0A=d0;
	d0_0=d0A;
	TM1=TMscore8_search(xtm, ytm, n_ali8, t0, u0, simplify_step, score_sum_method, &rmsd);
	TM_0=TM1;

	//normalized by length of structure B
	parameter_set4final(xlen+0.0);
	d0B=d0;
	TM2=TMscore8_search(xtm, ytm, n_ali8, t, u, simplify_step, score_sum_method, &rmsd);




	if(a_opt)
	{
		//normalized by average length of structures A, B
		Lnorm_0=(xlen+ylen)*0.5;
		parameter_set4final(Lnorm_0);
		d0a=d0;
		d0_0=d0a;
		TM3=TMscore8_search(xtm, ytm, n_ali8, t, u, simplify_step, score_sum_method, &rmsd);
		TM_0=TM3;
	}
	if(u_opt)
	{	
		//normalized by user assigned length		
		parameter_set4final(Lnorm_ass);		
		d0u=d0;		
		d0_0=d0u;
		Lnorm_0=Lnorm_ass;
		TM4=TMscore8_search(xtm, ytm, n_ali8, t, u, simplify_step, score_sum_method, &rmsd);	
		TM_0=TM4;
	}
	if(d_opt)
	{	
		//scaled by user assigned d0
		parameter_set4scale(ylen, d0_scale);
		d0_out=d0_scale;
		d0_0=d0_scale;
		//Lnorm_0=ylen;
		Lnorm_d0=Lnorm_0;
		TM5=TMscore8_search(xtm, ytm, n_ali8, t, u, simplify_step, score_sum_method, &rmsd);	
		TM_0=TM5;
	}

   
        
	output_results(xname, yname, xlen, ylen, t0, u0, TM1, TM2, rmsd, d0_out, m1, m2, n_ali8, n_ali, TM_0, Lnorm_0, d0_0);
                    
 









    //*********************************************************************************//
    //                            Done! Free memory                                    //
    //*********************************************************************************//           
    free_memory();
    delete [] invmap0;
    delete [] invmap;
	delete [] m1;
	delete [] m2;


    t2 = clock();    
    float diff = ((float)t2 - (float)t1)/CLOCKS_PER_SEC;
    printf("\nTotal running time is %5.2f seconds\n", diff);        

 	return 0;	
}
