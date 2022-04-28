# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <time.h>
# include <math.h>
# include <unistd.h>
# include "RKF78.h"
# include "randombits.h"

typedef struct {
    double x0;
    double phi;
    double beta;
    double lambda;
    double mu;
    double sigma;
    double delta;
} ODE_Parameters;

void ExitError(const char *miss, int errcode)
{
    fprintf (stderr, "\nERROR: %s.\nStopping...\n\n", miss);
    exit(errcode);
}

#define ElliotSigmoidSCALE 1000
#define TwoElliotSigmoidSCALE 2000

double ElliotSigmoid(double x, double sigma, double delta)
{
    x = sigma*(x-delta);
    return x/(ElliotSigmoidSCALE + fabs(x));
}

double Psi(double x, double mu, double sigma, double delta)
{
    //if(fabs(sigma) < ZeRoParsThreshold) return 1.0;
    double ES = ElliotSigmoid(x, sigma, delta);
    sigma *= delta; x /= delta;
    if(x < delta) {
    ES = ES * (x + (mu*(1.0-x)*(sigma + ElliotSigmoidSCALE)) / (sigma + TwoElliotSigmoidSCALE));
    }
    return ((1 - ES)*(sigma + TwoElliotSigmoidSCALE)) / (sigma*(1+mu) + TwoElliotSigmoidSCALE);
}

void MigrationODE(double t, double x, double *der, void *Params){
    ODE_Parameters *par = (ODE_Parameters *) Params; // Pointer cast to save typing and thinking
    *der = par->phi * x - par->beta*x*x - par->lambda*Psi(x, par->mu, par->sigma, par->delta);
}

#define HMAX 1.0
#define HMIN 1.e-6
#define RKTOL 1.e-8

int Generate_EDO_Prediction( double *xt, double x0,unsigned short number_of_years,ODE_Parameters *pars )
{
    register unsigned ty;
    xt[0] = x0; // Storing IC x(0)
    for(ty=1; ty < number_of_years; ty++) xt[ty] = 0.0;
    double t = 0.0, err, h = 1.e-3;
    int status;
    for(ty=1; ty < number_of_years; ty++)
    {
        while(t+h < ty)
        {
            status = RKF78(&t, &x0, &h, &err, HMIN, HMAX, RKTOL, pars, MigrationODE);
            if(status) return status;
        } // Adaptative stepsize h. To assure stopping at t = ty
        h = ty - t;
        status = RKF78(&t, &x0, &h, &err, HMIN, HMAX, RKTOL, pars, MigrationODE);
        if(status) return status;
        xt[ty] = x0;
    }
    return 0;
}

void printBits(size_t const size, void const *const ptr)
{
    unsigned char *b = (unsigned char*) ptr;
    unsigned char byte;
    int i,j;

    for (i=size-1; i>=0; i--){
        for (j=7; j>=0; j--){
            byte = (b[i]>>j)&1;
            printf("%u",byte);
        }
    }
    puts("");
}

typedef struct{
    unsigned long int x0;
    unsigned long int delta;
    unsigned long int lambda;
    unsigned long int mu;
    unsigned long int phi;
    unsigned long int sigma;
    double fitness;
} Individual;

// selects best individual out of t individuals
Individual tournament_selection(Individual *P, unsigned t, unsigned popsize)
{
    Individual best = P[(int)(uniform()*popsize)/1];
    Individual next;
    for(int i=0;i<=t;i++)
    {
        next = P[(int)(uniform()*popsize)/1];
        if(next.fitness < best.fitness)
        {
            best = next;
        }
    }
    return best;
}

// convert to phenotype! 
ODE_Parameters Phenotype(Individual individual)
{
    ODE_Parameters phenotype;
    phenotype.x0 = (double) individual.x0 * 16600/((1UL<<21)-1);
    phenotype.delta = (double) individual.delta * 25000/((1UL<<15)-1);
    phenotype.lambda = (double) individual.lambda * 3000/((1UL<<25)-1);
    phenotype.mu = (double) individual.mu * 20/((1UL<<25)-1);
    phenotype.phi = (double) individual.phi * 100.35/((1UL<<34)-1)-100;
    phenotype.sigma = (double) individual.sigma * 1000/((1UL<<17)-1);
    phenotype.beta = 0.000024382635446;
    return phenotype;
}

double fitness(Individual individual)
{
    double data[] = {15329,14177,13031,9762,11271,8688,7571,6983,4778,2067,1586,793};
    int i;
    double sum = 0;
    unsigned years = 12;
    double x[years],tmp;
    ODE_Parameters pars = Phenotype(individual);
    Generate_EDO_Prediction(x,pars.x0,years,&pars);
    for(i=0;i<years;i++)
    {
        tmp = pow(x[i]-data[i],2);
        sum=sum+tmp;
    }
    return sum;
}

void mutation1(unsigned int *f, float prob){
    if(uniform()<prob)*f = (*f)^(1U<<((unsigned char) uniform()*8*sizeof(*f)));
}

void BitFlipMutation(unsigned long int *f,float prob, const unsigned char len)
{
   register unsigned char i;
   for(i=0; i<len;i++) if(uniform() < prob) *f = (*f)^(1U<<i);
}

void OnePointCrossover(unsigned long int p1, unsigned long int p2, unsigned long int *f1, unsigned long int *f2, const unsigned char bit_cutoff){
    const unsigned char split_position = uniform()*bit_cutoff + 1;
    const unsigned long int mask = 0xFFFFFFFFFFFFFFFFUL <<split_position;
    *f1 = (p1 & mask) | (p2 & ~mask);
    *f2 = (p2 & mask) | (p1 & ~mask);
}

void PrintInd(Individual ind){
    printf("X0=%lu\nDelta=%lu\nLambda=%lu\nMu=%lu\nPhi=%lu\nSigma=%lu\nFitness=%f\n",ind.x0,ind.delta,ind.lambda,ind.mu,ind.phi,ind.sigma,ind.fitness);
    printf("\n-----------------------\n");
}

Individual adjust(Individual individual, int r)
{
    individual.x0 += (int)(r*uniform());
    individual.delta += (int) (r*uniform());
    individual.mu += (int) (r*uniform());
    individual.lambda += (int)(r*uniform());
    individual.sigma += (int)(r*uniform());
    individual.phi += (int)(r*uniform());

    if(individual.x0>0) individual.x0-=(int)(r*uniform());
    if(individual.delta>0) individual.delta -= (int)(r*uniform());
    if(individual.mu>0) individual.mu -= (int)(r*uniform());
    if(individual.lambda>0) individual.lambda -= (int)(r*uniform());
    if(individual.sigma>0) individual.sigma -= (int)(r*uniform());
    if(individual.phi>0) individual.phi -= (int)(r*uniform());

    return individual;
}

//Exploit the solution you have become with the Genetic algorithm by adjusting it a little bit and compare the fitnesses of both. 
Individual HillClimbing(Individual S)
{
    Individual best = S;
    Individual R;
    Individual W;
    int r = 1000;
    int n = 5000, iter = 0, static_iteration = 0;

    do
    {
        R = adjust(S,r);
        R.fitness = fitness(R);
        for(int i =0; i<n-1; i++)
        {
            W = adjust(R,r);
            if((W.fitness = fitness(W))<R.fitness)R=W;
        }
        S=R;
        if(S.fitness<best.fitness)
        {
            best=S;
            static_iteration = 0;
        }
        printf("Iteration= %d Best_fitness= %f\n",iter,best.fitness);
        iter++;
        static_iteration++;
        if(static_iteration > 10 && r>0)r-=100;
    } while (iter<1000 && static_iteration <50);
    return best;
}

int Ind_comparator(const void *v1, const void *v2)
{
    const Individual *p1 = (Individual *)v1;
    const Individual *p2 = (Individual *)v2;
    if (p1->fitness < p2->fitness)
        return -1;
    else if (p1->fitness > p2->fitness)
        return +1;
    else
        return 0;
}

Individual GA(int popsize, int generations)
{
    int generation = 0, static_generations=0;
    int Index;
    int half_popsize = popsize/2;
    float percentage1 = 0.04;
    Individual *P;
    Individual *C;
    Individual Parent1;
    Individual Parent2;
    Individual overall_best;
    Individual best_of_gen;
    if((P=(Individual*)malloc(popsize*sizeof(Individual)))==NULL)
        ExitError("Couldn't allocate memory for individuals",32);
    if((C=(Individual*)malloc(popsize*sizeof(Individual)))==NULL)
        ExitError("Couldn't allocate memory for individuals",32);

    randomize();
    for(int i=0; i<popsize; i++)
    {
        P[i].x0 = ULONGran(21);
        P[i].delta = ULONGran(15);
        P[i].lambda = ULONGran(25);
        P[i].mu = ULONGran(25);
        P[i].phi = ULONGran(34);
        P[i].sigma = ULONGran(17);
        P[i].fitness = fitness(P[i]);
    }
    
    overall_best = P[0];
    overall_best.fitness = fitness(overall_best);
    
    do
    {
        unsigned int tournament_selector_param = (unsigned int) (((0.1*popsize*generation)/generations) +2);
        float mutation_prob = 0.09-(0.085/generations)*generation;
        best_of_gen = P[0];
        best_of_gen.fitness = fitness(best_of_gen);
        for(int i=0; i<popsize;i++)
        {
            if(P[i].fitness<overall_best.fitness)
            {
                overall_best = P[i];
                static_generations = 0;
            }

            if(P[i].fitness<best_of_gen.fitness)
            {
                best_of_gen = P[i];
            }
        }
    printf("Generation= %d Best_fitness= %f\n",generation, best_of_gen.fitness);

        for(int i=0;i<half_popsize;i++)
        {
            Index = half_popsize +i;
            Parent1=tournament_selection(P,tournament_selector_param,popsize);
            Parent2 = tournament_selection(P,tournament_selector_param,popsize);
            OnePointCrossover(Parent1.delta,Parent2.delta,&C[i].delta,&C[Index].delta,15);
            OnePointCrossover(Parent1.x0,Parent2.x0,&C[i].x0,&C[Index].x0,21);
            OnePointCrossover(Parent1.mu,Parent2.mu,&C[i].mu,&C[Index].mu,25);
            OnePointCrossover(Parent1.sigma,Parent2.sigma,&C[i].sigma,&C[Index].sigma,17);
            OnePointCrossover(Parent1.lambda,Parent2.lambda,&C[i].lambda,&C[Index].lambda,25);
            OnePointCrossover(Parent1.phi,Parent2.phi,&C[i].phi,&C[Index].phi,34);
            BitFlipMutation(&C[i].delta,mutation_prob,15);
            BitFlipMutation(&C[i].x0,mutation_prob,21);
            BitFlipMutation(&C[i].mu,mutation_prob,25);
            BitFlipMutation(&C[i].sigma,mutation_prob,17);
            BitFlipMutation(&C[i].lambda,mutation_prob,25);
            BitFlipMutation(&C[i].phi,mutation_prob,34);
            BitFlipMutation(&C[Index].delta,mutation_prob,15);
            BitFlipMutation(&C[Index].x0,mutation_prob,21);
            BitFlipMutation(&C[Index].mu,mutation_prob,25);
            BitFlipMutation(&C[Index].sigma,mutation_prob,17);
            BitFlipMutation(&C[Index].lambda,mutation_prob,25);
            BitFlipMutation(&C[Index].phi,mutation_prob,34);
        }

        for(int i=0; i<popsize; i++)
        {
            C[i].fitness = fitness(C[i]);
        }

        //sort populations here
        qsort(P, popsize, sizeof(Individual), Ind_comparator);
        qsort(C, popsize, sizeof(Individual), Ind_comparator);

        //change population here
        for(int i=0;i<popsize;i++){ 
            if(i<percentage1*popsize){
                P[i]=P[i];
            }
            else{
                P[i]=C[i-(int)(percentage1*popsize)];
            }
            
        }
        generation++;
        static_generations++;
    } while (generation<generations && static_generations<25);
    return overall_best;
}


int main (int argc, char *argv[])
{ 
    randomize();
    unsigned popsize =10000, generations = 60;
    unsigned years = 12;
    
    double x[years];
    int i;
    Individual best = GA(popsize,generations);
    best = HillClimbing(best);
    ODE_Parameters parameters = Phenotype(best);

    Generate_EDO_Prediction(x,parameters.x0,years,&parameters);
    for(i=0; i<years;i++)
    {
        printf("x[%d]=%f\n",i,x[i]);
    }

    printf("X0= %f\nLambda= %f\nPhi= %f\nMu=%f\nSigma= %f\nDelta= %f\n",parameters.x0,parameters.lambda,parameters.phi,parameters.mu,parameters.sigma,parameters.delta);
}