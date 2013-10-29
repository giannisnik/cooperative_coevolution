#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define PI 3.14
#define POPULATION 100
#define SIZE 16
#define RASTRIGIN 0
#define SCHWEFEL 1
#define GRIEWANGK 2
#define ACKLEY 3
#define ROSENBROCK 4

char** initialize_population(int);
double evaluate_individual(char*, int);
char** select_new_population(char**, int, double*);
void mutate_population(char**, int);
char** recombine_population(char**, int);
void GA(int);

char*** initialize_species_populations(int);
char* create_individual(int, int, char*, char*);
char** select_new_species_population(char**, char*,int, int, int);
void mutate_species_population(char**);
char** recombine_species_population(char**);
void CCGA(int);

int main()
{
    CCGA(RASTRIGIN);
    //GA(ROSENBROCK);

    return 0;
}



char** initialize_population(int function){
    int i,j,random,function_variables;
    char** population = malloc(POPULATION*sizeof(char*));

    switch(function){
        case RASTRIGIN:
            function_variables = 20;
            break;
        case SCHWEFEL:
            function_variables = 10;
            break;
        case GRIEWANGK:
            function_variables = 10;
            break;
        case ACKLEY:
            function_variables = 30;
            break;
        case ROSENBROCK:
            function_variables = 2;
            break;
    }

    for(i=0;i<POPULATION;i++){
        population[i] = malloc(SIZE*function_variables*sizeof(char));
    }

    for(i=0;i<POPULATION;i++){
        for(j=0;j<SIZE*function_variables;j++){
            random = rand()%2;
            population[i][j] = (char)(random+48);
        }
    }
    return population;
}



double evaluate_individual(char* individual, int function){
    int i,j,number;
    double fitness;

    if(function==RASTRIGIN){
        double x[20] = {0};

        for(i=0;i<20;i++){
            for(j=0;j<SIZE;j++){
                number = individual[(SIZE*i)+j] - '0';
                x[i] += number*pow(2, (15-j));
            }
            x[i] *= 0.00015625;
            x[i] -= 5.12;
        }
        fitness = 60;
        for(i=0;i<20;i++){
            fitness += pow(x[i], 2) - 3*cos(2*PI*x[i]);
        }
    }
    else if(function==SCHWEFEL){
        double x[10] = {0};

        for(i=0;i<10;i++){
            for(j=0;j<SIZE;j++){
                number = individual[(SIZE*i)+j] - '0';
                x[i] += number*pow(2, (15-j));
            }
            x[i] *= 0.0152587890625;
            x[i] -= 500;
        }
        fitness = 4189.829;
        for(i=0;i<10;i++){
            fitness += x[i]*sin(sqrt(fabs(x[i])));
        }
    }
    else if(function==GRIEWANGK){
        double x[10] = {0};
        double product = 1;

        for(i=0;i<10;i++){
            for(j=0;j<SIZE;j++){
                number = individual[(SIZE*i)+j] - '0';
                x[i] += number*pow(2, (15-j));
            }
            x[i] *= 0.018310546875;
            x[i] -= 600;
        }
        fitness = 1;
        for(i=0;i<10;i++){
            fitness += (pow(x[i], 2)/4000);
            product *= cos(x[i]/sqrt(i+1));
        }
        fitness -= product;
    }
    else if(function==ACKLEY){
        double x[30] = {0};
        double sum1,sum2;
        sum1 = sum2 = 0;

        for(i=0;i<30;i++){
            for(j=0;j<SIZE;j++){
                number = individual[(SIZE*i)+j] - '0';
                x[i] += number*pow(2, (15-j));
            }
            x[i] *= 0.00091552734375;
            x[i] -= 30;
        }
        for(i=0;i<30;i++){
            sum1 += pow(x[i], 2);
            sum2 += cos(2*PI*x[i]);
        }
        fitness = 20 + exp(1) - 20*exp(-0.2*sqrt(sum1/30)) - exp(sum2/30);
    }
    else if(function==ROSENBROCK){
        double x[2] = {0};

        for(i=0;i<2;i++){
            for(j=0;j<SIZE;j++){
                number = individual[(SIZE*i)+j] - '0';
                x[i] += number*pow(2, (15-j));
            }
            x[i] *= 0.0000625;
            x[i] -= 2.048;
        }
        fitness = 100*pow((pow(x[0], 2)-x[1]),2) + pow((1-x[0]),2);

    }
    return fitness;
}



char** select_new_population(char** population, int function, double* overall_best){
    int i,j,k,function_variables,function_evaluations;
    char** new_population;
    double best_individual,random, fitness_sum = 0;
    double fitness[POPULATION];
    double limit[POPULATION];
    double limit_sum = 0;

    switch(function){
        case RASTRIGIN:
            function_variables = 20;
            break;
        case SCHWEFEL:
            function_variables = 10;
            break;
        case GRIEWANGK:
            function_variables = 10;
            break;
        case ACKLEY:
            function_variables = 30;
            break;
        case ROSENBROCK:
            function_variables = 2;
            break;
    }

    function_evaluations = 0;
    best_individual = evaluate_individual(population[0], function);

    for(i=0;i<POPULATION;i++){
        fitness[i] = evaluate_individual(population[i], function);
        function_evaluations++;
        if(fitness[i]<best_individual){
            best_individual = fitness[i];
        }
        fitness_sum += fitness[i];
    }

    for(i=0;i<POPULATION;i++){
        limit_sum += 1/(fitness[i]/fitness_sum);
        limit[i] = limit_sum;
    }

    new_population = malloc(POPULATION*sizeof(char*));
    for(i=0;i<POPULATION;i++){
        j = 0;
        random = limit_sum*((double)rand()/(double)RAND_MAX);
        while(random>limit[j]){
            j++;
        }
        char* new_individual = malloc(SIZE*function_variables*sizeof(char));
        for(k=0;k<SIZE*function_variables;k++){
            new_individual[k] = population[j][k];
        }
        new_population[i] = new_individual;
    }

    for(i=0;i<POPULATION;i++){
        free(population[i]);
    }

    free(population);

    //printf("Best Individual: %f\n",best_individual);
    //printf("Function Evaluations: %d\n",function_evaluations);
    if((best_individual)<(*overall_best))
        (*overall_best) = best_individual;

    return new_population;
}



void mutate_population(char** population, int function){
    int random,i,j,function_variables,length_of_chromosome;

    switch(function){
        case RASTRIGIN:
            function_variables = 20;
            break;
        case SCHWEFEL:
            function_variables = 10;
            break;
        case GRIEWANGK:
            function_variables = 10;
            break;
        case ACKLEY:
            function_variables = 30;
            break;
        case ROSENBROCK:
            function_variables = 2;
            break;
    }

    length_of_chromosome = SIZE*function_variables;

    for(i=0;i<POPULATION;i++){
        for(j=0;j<length_of_chromosome;j++){
            random = rand()%length_of_chromosome;
            if(random==0){
                if(population[i][j]=='0')
                    population[i][j] = '1';
                else
                    population[i][j] = '0';
            }
        }
    }
}



char** recombine_population(char** population, int function){
    int index,random_int,i,j,temp,function_variables,remaining_population,point1,point2;
    double random;
    char temp_char;
    char** new_population;

    switch(function){
        case RASTRIGIN:
            function_variables = 20;
            break;
        case SCHWEFEL:
            function_variables = 10;
            break;
        case GRIEWANGK:
            function_variables = 10;
            break;
        case ACKLEY:
            function_variables = 30;
            break;
        case ROSENBROCK:
            function_variables = 2;
            break;
    }

    new_population = malloc(POPULATION*sizeof(char*));
    remaining_population = POPULATION;

    for(i=0;i<POPULATION/2;i++){
        index = i*2;
        new_population[index] = malloc(SIZE*function_variables*sizeof(char));
        new_population[index+1] = malloc(SIZE*function_variables*sizeof(char));
        random_int = rand()%remaining_population;
        for(j=0;j<SIZE*function_variables;j++){
            new_population[index][j] = population[random_int][j];
        }
        if(random_int!=(remaining_population-1)){
            for(j=0;j<SIZE*function_variables;j++){
                population[random_int][j] = population[remaining_population-1][j];
            }
        }
        remaining_population--;

        random_int = rand()%remaining_population;
        for(j=0;j<SIZE*function_variables;j++){
            new_population[index+1][j] = population[random_int][j];
        }
        if(random_int!=(remaining_population-1)){
            for(j=0;j<SIZE*function_variables;j++){
                population[random_int][j] = population[remaining_population-1][j];
            }
        }
        remaining_population--;

        random = ((double)rand()/(double)RAND_MAX);
        if(random<=0.6){
            point1 = (rand()%((SIZE*function_variables)-1))+1;
            point2 = (rand()%((SIZE*function_variables)-1))+1;
            if(point1==point2){
                while(point1==point2){
                    point2 = (rand()%((SIZE*function_variables)-1))+1;
                }
            }

            if(point1>point2){
                temp = point1;
                point1 = point2;
                point2 = temp;
            }

            for(j=point1;j<point2;j++){
                temp_char = new_population[index][j];
                new_population[index][j] = new_population[index+1][j];
                new_population[index+1][j] = temp_char;
            }
        }
    }

    for(i=0;i<POPULATION;i++){
        free(population[i]);
    }

    free(population);

    return new_population;
}



void GA(int function){
    char** population;
    double best;
    double* best_ptr;
    int i;
    FILE *fp;

    fp = fopen("results.txt", "w");

    population = initialize_population(function);

    best = evaluate_individual(population[0], function);
    best_ptr = &best;

    for(i=0;i<1000;i++){
        population = select_new_population(population, function, best_ptr);
        printf("Best Individual Found: %f\n",*best_ptr);
        fprintf(fp, "%f\n",*best_ptr);
        population = recombine_population(population, function);
        mutate_population(population, function);
    }

    fclose(fp);

}












char*** initialize_species_populations(int function){
    int i,j,k,random,species;
    char*** population;

    switch(function){
        case RASTRIGIN:
            species = 20;
            break;
        case SCHWEFEL:
            species = 10;
            break;
        case GRIEWANGK:
            species = 10;
            break;
        case ACKLEY:
            species = 30;
            break;
        case ROSENBROCK:
            species = 2;
            break;
    }

    population = malloc(species*sizeof(char**));

    for(i=0;i<species;i++){
        population[i] = malloc(POPULATION*sizeof(char*));
    }

    for(i=0;i<species;i++){
        for(j=0;j<POPULATION;j++){
            population[i][j] = malloc(SIZE*sizeof(char));
        }
    }

    for(i=0;i<species;i++){
        for(j=0;j<POPULATION;j++){
            for(k=0;k<SIZE;k++){
                random = rand()%2;
                population[i][j][k] = (char)(random+48);
            }
        }
    }

    return population;
}



char* create_individual(int species, int current_species, char* best_species, char* current_species_individual){
    int j,k,position;
    char* individual;

    individual = malloc(SIZE*species*sizeof(char));
    position = 0;

    for(j=0;j<current_species;j++){
        for(k=0;k<SIZE;k++){
            individual[position] = best_species[position];
            position++;
        }
    }

    for(k=0;k<SIZE;k++){
        individual[position] = current_species_individual[k];
        position++;
    }

    for(j=current_species+1;j<species;j++){
        for(k=0;k<SIZE;k++){
            individual[position] = best_species[position];
            position++;
        }
    }

    return individual;
}



char** select_new_species_population(char** population, char* best_species,int function, int species, int current_species){
    int i,j,k,function_evaluations;
    char** new_population;
    double random, fitness_sum = 0;
    double fitness[POPULATION];
    double limit[POPULATION];
    double limit_sum = 0;
    char* individual;

    function_evaluations = 0;

    for(i=0;i<POPULATION;i++){
        individual = create_individual(species, current_species, best_species, population[i]);
        fitness[i] = evaluate_individual(individual, function);
        free(individual);
        function_evaluations++;
        fitness_sum += fitness[i];
    }

    for(i=0;i<POPULATION;i++){
        limit_sum += 1/(fitness[i]/fitness_sum);
        limit[i] = limit_sum;
    }

    new_population = malloc(POPULATION*sizeof(char*));
    for(i=0;i<POPULATION;i++){
        j = 0;
        random = limit_sum*((double)rand()/(double)RAND_MAX);
        while(random>limit[j]){
            j++;
        }
        char* new_individual = malloc(SIZE*sizeof(char));
        for(k=0;k<SIZE;k++){
            new_individual[k] = population[j][k];
        }
        new_population[i] = new_individual;
    }

    for(i=0;i<POPULATION;i++){
        free(population[i]);
    }

    free(population);

    return new_population;
}


void mutate_species_population(char** population){
    int random,i,j;

    for(i=0;i<POPULATION;i++){
        for(j=0;j<SIZE;j++){
            random = rand()%SIZE;
            if(random==0){
                if(population[i][j]=='0')
                    population[i][j] = '1';
                else
                    population[i][j] = '0';
            }
        }
    }
}



char** recombine_species_population(char** population){
    int index,random_int,i,j,temp,remaining_population,point1,point2;
    double random;
    char temp_char;
    char** new_population;

    new_population = malloc(POPULATION*sizeof(char*));
    remaining_population = POPULATION;

    for(i=0;i<POPULATION/2;i++){
        index = i*2;
        new_population[index] = malloc(SIZE*sizeof(char));
        new_population[index+1] = malloc(SIZE*sizeof(char));
        random_int = rand()%remaining_population;
        for(j=0;j<SIZE;j++){
            new_population[index][j] = population[random_int][j];
        }
        if(random_int!=(remaining_population-1)){
            for(j=0;j<SIZE;j++){
                population[random_int][j] = population[remaining_population-1][j];
            }
        }
        remaining_population--;

        random_int = rand()%remaining_population;
        for(j=0;j<SIZE;j++){
            new_population[index+1][j] = population[random_int][j];
        }
        if(random_int!=(remaining_population-1)){
            for(j=0;j<SIZE;j++){
                population[random_int][j] = population[remaining_population-1][j];
            }
        }
        remaining_population--;

        random = ((double)rand()/(double)RAND_MAX);

        if(random<=0.6){
            point1 = (rand()%(SIZE-1))+1;
            point2 = (rand()%(SIZE-1))+1;
            if(point1==point2){
                while(point1==point2){
                    point2 = (rand()%(SIZE-1))+1;
                }
            }

            if(point1>point2){
                temp = point1;
                point1 = point2;
                point2 = temp;
            }

            for(j=point1;j<point2;j++){
                temp_char = new_population[index][j];
                new_population[index][j] = new_population[index+1][j];
                new_population[index+1][j] = temp_char;
            }
        }
    }

    for(i=0;i<POPULATION;i++){
        free(population[i]);
    }

    free(population);

    return new_population;
}




void CCGA(int function){
    int species,i,j,k,best_individual,random,counter;
    char*** populations;
    char* best_species;
    double best_cost,cost,best_found_cost,current_cost;
    char* individual;

    switch(function){
        case RASTRIGIN:
            species = 20;
            break;
        case SCHWEFEL:
            species = 10;
            break;
        case GRIEWANGK:
            species = 10;
            break;
        case ACKLEY:
            species = 30;
            break;
        case ROSENBROCK:
            species = 2;
            break;
    }

    FILE *fp;

    fp = fopen("results.txt", "w");

    populations = initialize_species_populations(function);

    best_species = malloc(SIZE*species*sizeof(char));
    for(i=0;i<species;i++){
        random = rand()%POPULATION;
        for(j=0;j<SIZE;j++){
            best_species[(i*16)+j] = populations[i][random][j];
        }
    }

    best_found_cost = evaluate_individual(best_species, function);

    counter = 1;

    for(i=0;i<50;i++){
        for(j=0;j<species;j++){
            populations[j] = select_new_species_population(populations[j], best_species, function, species, j);
            populations[j] = recombine_species_population(populations[j]);
            mutate_species_population(populations[j]);

            individual = create_individual(species, j, best_species, populations[j][0]);
            best_cost = evaluate_individual(individual, function);
            best_individual = 0;
            free(individual);
            for(k=1;k<POPULATION;k++){
                individual = create_individual(species, j, best_species, populations[j][k]);
                cost = evaluate_individual(individual, function);
                if(cost<best_cost){
                    best_cost = evaluate_individual(individual, function);
                    best_individual = k;
                }
                free(individual);
            }

            for(k=0;k<SIZE;k++){
                best_species[(j*16)+k] = populations[j][best_individual][k];
            }

            current_cost = evaluate_individual(best_species, function);
            if(current_cost<best_found_cost){
                best_found_cost = current_cost;
            }

            printf("Best Individual Found at %d: %f\n", counter, best_found_cost);
            fprintf(fp, "%f\n",best_found_cost);
            counter++;
        }
    }

    fclose(fp);
}

