/*
*  Find how much one sweep shifts the ancestry at another locus in space
* This code is to find the statistics of the widths of sweeping and focal allele when they meet.
*
* Modification of Mike McLaren's meta_hap_stepping_stone.c
*
*
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_randist.h>
#include<time.h>
#include<vector>
#include<map>

// Turn off/on logging of genotype counts in alldata.dat file with 0/1
#define LOGGING 0
// Haplotype indices
#define ab 0
#define Ab 1
#define aB 2
#define AB 3
#define L 2000
#define Lsweep 200
// Number of haplotypes (4) and number of genotypes (10)
#define HTYPES 4

// Global variables
const gsl_rng_type * T;
gsl_rng * R;

void next_gen(unsigned int **n, double**x, double**xmig, double w[4], double r,
	double mig, unsigned int N, int wrap, int *nose1, int *nose1b, int *nose2,int *nose2b,int *nose2c,int *nose2d) {
	double D; // coefficient of linkage disequilibrium
	int i, j, k, flag = 1,flag2=1;
	// Recombination
	for (i = 0; i < L; i++) {
		// Calculate haplotype frequencies
		for (j = 0; j < 4; j++) {
			x[i][j] = ((float)n[i][j]) / N;
		}

		// Recombination
		D = x[i][0] * x[i][3] - x[i][1] * x[i][2];
		x[i][0] -= r * D;
		x[i][1] += r * D;
		x[i][2] += r * D;
		x[i][3] -= r * D;
	}

	// Migration
	if (wrap) {
		// Wraparound (no boundaries), so that first and last demes are neighbors
		for (j = 0; j < 4; j++) {
			// leftmost deme:
			xmig[0][j] = (1 - mig) * x[0][j] + 0.5 * mig * (x[L - 1][j] + x[1][j]);
			//rightmost deme:
			xmig[L - 1][j] = (1 - mig) * x[L - 1][j] + 0.5 * mig * (x[L - 2][j] + x[0][j]);
		}
	}
	else {
		// 0 and L are really boundaries
		for (j = 0; j < 4; j++) {
			// leftmost deme:
			xmig[0][j] = (1 - 0.5 * mig) * x[0][j] + 0.5 * mig * x[1][j];
			//rightmost deme:
			xmig[L - 1][j] = (1 - 0.5*mig) * x[L - 1][j] + 0.5 * mig * x[L - 2][j];
		}
	}
	// All interior demes
	for (i = 1; i < L - 1; i++) {
		for (j = 0; j < 4; j++) {
			xmig[i][j] = (1 - mig) * x[i][j] + 0.5 * mig * (x[i - 1][j] + x[i + 1][j]);
		}
	}

	// Selection and sampling within demes
	for (i = 0; i < L; i++) {
		// Selection (within deme)
		for (j = 0; j < 4; j++) {
			if (xmig[i][j] == 1.0) { // a genotype is fixed in the deme, so no need to do sampling
				for (k = 0; k < 4; k++) n[i][k] = 0;
				n[i][j] = N;
				break;
			}
			x[i][j] = xmig[i][j] * w[j];
		}
		// Normalization of x[i] done in multinomial function, so dividing by
		// average fitness is not necessary
		gsl_ran_multinomial(R, 4, N, x[i], n[i]);
		if (n[i][1] == N)
		{
			*nose1 = i;
		}
		if (n[i][1] < N && n[i][1]>0)
		{
			*nose1b = i;
		}
                if(n[i][2]>0)
                   *nose2c=i; 
                if (flag && i<L && n[i][2]>0)
                {
                    *nose2b=i;
                    flag=0;
                }
                if (flag2 && n[i][2]==N)
                {
                    *nose2=i;
                    flag2=0;
                }
                if (n[i][2]==N)
                {
                    *nose2d=i;
                }


	}
}

void loadData(std::vector< std::vector<int> >& dist, char*filename)
{
	FILE *fp;
	int c;
	char ch[BUFSIZ];
	int field;
	int m;
	char *start;
	int i = 0;
	fp = fopen(filename, "r");

	while (!feof(fp))
	{

		while (fgets(ch, sizeof ch, fp) != NULL) {
			start = ch;

			dist.push_back(std::vector <int>());
			while (sscanf(start, "%d%n", &field, &m) == 1)
			{
				dist[i].push_back(field);
				start += m;
			}
			i++;
		}

	}
	fclose(fp);
}

/*function roundint to round a double to an integer*/
int roundint(double x) {
	int result;

	result = floor(x);
	if (x - result >= 0.5) result++;

	return result;
} /*roundint*/



int main(int argc, char *argv[]) {
	int trails;
	long SEED;
	int wrap; //whether the range wraps around
	double r, mig, s1, s2; // recomb rate, migration rate, selection
	double w[4]; // fitnesses of the four genotypes
	unsigned int N, t, tfinal; // deme size, number of demes, deme to paint with neutral allele, time, max number of generations
	unsigned int ntot[4],nfocal; // numbers of the four genotypes in the total population
	FILE *datafile;
	FILE *paramfile;
	char *outfile;
	char *infile;
	char *outfile1 = (char*)malloc(1000);
	int i, j, k;
	int flag, nose1,nose1b,midnose1,nose2, nose2b,nose2c,nose2d,foc;
	int distance;
	int key;
        int tsim=1200,tcount;
        int flagrec;
	std::vector< std::vector<int> > distribution;
	gsl_ran_discrete_t *g = NULL;
	j = 1;
	trails = atof(argv[j++]);
	distance = atof(argv[j++]);
    foc=atof(argv[j++]);
	wrap = roundint(atof(argv[j++]));
	mig = atof(argv[j++]);
	N = (unsigned int)roundint(atof(argv[j++]));
	s1 = 0.05;
	s2 = atof(argv[j++]);
	r = atof(argv[j++]);
	infile = argv[j++];
	outfile = argv[j++];
	w[0] = 1;
	w[1] = 1 + s1;
	w[2] = 1 + s2;
	w[3] = 1 + s1 + s2;

	strcpy(outfile1, infile);
	loadData(distribution, strcat(outfile1, ".txt")); //load the distribution of focal allele occupancies

	double *P = (double*)malloc(distribution.size()*sizeof(double)); //initializing the weigths for the occupancies distribution
	for (i = 0; i < distribution.size(); i++)
		P[i] = 1.;
	g = gsl_ran_discrete_preproc(distribution.size(), P);

	tfinal = 500000;// (unsigned int)roundint(atof(argv[j++]));
	SEED = (-1)*time(NULL);// atof(argv[j++]);
	strcpy(outfile1, outfile);

	datafile = fopen(strcat(outfile, ".txt"), "w");
	strcpy(outfile, outfile1);

	paramfile = fopen(strcat(outfile, "_params.txt"), "w");
	strcpy(outfile, outfile1);
	fprintf(paramfile, " trails = %d\n distance = %d\n foc= %d\n bc = %d\n mig = %f\n N = %d\n s2 = %f\n r = %lf\n", trails, distance,foc,wrap, mig, N, s2, r);
	fclose(paramfile);


	gsl_rng_env_setup();
	T = gsl_rng_default;
	R = gsl_rng_alloc(T);
	gsl_rng_set(R, SEED);
 
	// Initialize population:
	unsigned int **n = (unsigned int**)malloc(L * sizeof(unsigned int*)); // numbers of the four genotypes in each deme
	double **x = (double**)malloc(L * sizeof(double*));
	double **xmig = (double**)malloc(L * sizeof(double*));
	for (i = 0; i < L; i++)
	{
		n[i] = (unsigned int*)malloc(4 * sizeof(unsigned int));
		x[i] = (double*)malloc(4 * sizeof(double));
		xmig[i] = (double*)malloc(4 * sizeof(double));
	}
        
	for (k = 0; k < trails; k++)
	{
                tcount=0;
                flag = 0;
                flagrec=0;
		key = gsl_ran_discrete(R, g);
                nose2=-1;
		nose2b = distribution[key][0];
                nose2c=-1;
                nose2d=-1;
		// Left-most deme fixed for sweeping allele:
		for (i = 0; i < 10; i++)
		{
			n[i][0] = 0;
			n[i][1] = N;
			n[i][2] = 0;
			n[i][3] = 0;
		}
		for (i = 10;i < L; i++)
		{
			n[i][0] = N;
			n[i][1] = 0;
			n[i][2] = 0;
			n[i][3] = 0;
		}
		nose1 = 9;
                nose1b = nose1;
                midnose1 = nose1;
		// Run until one of the alleles fixes (or goes extinct):
		for (t = 0; t < tfinal; t++) {
                        if ((foc-midnose1)<=distance && !flag)
                            flag=1;
                        if (flag==1)
                            tcount++;
			if (flag && tcount==tsim && !flagrec)
			{
                                flag=0;
                                flagrec=1;
				for (i = 1; i < distribution[key].size(); i++)
				{
					n[nose2b + i - 1][2] = distribution[key][i];
					n[nose2b + i - 1][0] = N - n[nose2b + i - 1][2];
				}
			}
			if (flagrec)
			{
				fprintf(datafile,"%d %d %d %d %d %d %d %d %d\n",nose1,nose1b,nose2,nose2b,nose2c,nose2d,ntot[1],ntot[2],nfocal);
			}
                       if(nose1b>=nose2b)
                       {
                         break;
                       }
			// Evolve for one generation
			next_gen(n, x, xmig, w, r, mig, N, wrap, &nose1,&nose1b,&nose2,&nose2b,&nose2c,&nose2d);
                        midnose1=(nose1+nose1b)/2;
			// Add up genotype counts within each deme
			for (j = 0; j < 4; j++) {
				ntot[j] = 0;
			}
                        nfocal=0;
			for (i = 0; i < L; i++) {
                                nfocal+=n[i][2];
				for (j = 0; j < 4; j++) {
					if (n[i][j] == N)
						ntot[j] += 1;
				}
			}

			// Stop if one of the alleles is fixed or extinct:
			//if ((ntot[1] + ntot[3] == L) || (ntot[1] + ntot[3] == 0) /*|| (ntot[2] + ntot[3] == L) || (ntot[2] + ntot[3] == 0)*/)
			if (ntot[0] == L || ntot[1] == L || ntot[2] == L || (ntot[2] + ntot[3]) == L || ntot[3] == L)
                        {
                                fprintf(datafile, "%d\n", (ntot[2] == L || (ntot[2] + ntot[3]) == L || ntot[3] == L));
				break;
                        }
		}

	}
	fclose(datafile);
        for (i = 0; i < L; i++)
	{
            free(n[i]);
            free(x[i]);
            free(xmig[i]);
	}
        free(n);
        free(x);
        free(xmig);
        free(outfile1);
        free(P);
        free(datafile);
        free(paramfile);
        std::vector< std::vector<int> >().swap(distribution);
        printf("done");
	return 0;
}
