#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <random>
//#include <gsl/gsl_math.h>
//#include <gsl/gsl_fit.h>

#define temp_fsize 60
#define class_size 10
#define height 10.0                          // ft
#define width 8.0                            // ft
#define gravity 32.2                         // ft.s^-2
#define discharge_coefficient 0.95            // 
#define room_volume 30000                    // ft^3
#define air_density 0.0752                   //lbs.ft^-3
#define LINEAR_SPLINE 1
#define QUADRATIC_SPLINE 2

struct section{
	struct tm start; //start time
	struct tm end;   //end time
	int num;         //number of students in class section
};


//function prototypes
double calc_heat_escaped(double , double , double );
double calc_new_temperature(double , double );
int calc_expected_people_opening_door(struct tm , struct section [] , int );
int calc_expected_time_door_open(int );
double* run_monte_carlo_simulation(struct tm , double [], double [], int , struct section [], int );
double normal_pdf(double , double , double ); //normal distribution density function
double fTok(double );
double* b_spline(int , int , double [], double [], int , int );
void between(int , int *, int *, double [], int );
int nChoosek( int , int  );


using namespace std;

int main()
{
	time_t rawtime;
    time ( &rawtime );	
	ifstream out_file;
	ifstream in_file;
	ifstream schedule_file;
	ofstream output_file;
	//OPEN INPUT FILES
	out_file.open("out_temp.dat");
	in_file.open("in_temp.dat");
	schedule_file.open("schedule.dat");
	output_file.open("output.csv");
	double temp_in[temp_fsize];
	double temp_out[temp_fsize];
	struct section sections[class_size];
	cout << "PROGRAM STARTING...\n";
	cout << "READING TEMPERATURE FILES\n";
	//READ TEMP DATA FROM FILES
	for(int i = 0; i < temp_fsize; i++)
	{
		out_file >> temp_out[i];
		in_file >> temp_in[i];
	}
	cout << "READING CLASS SCHEDULES FILES\n";
	//READ CLASS SCHEDULES
	for(int i = 0; i < class_size; i++)
	{
		sections[i].start = *(localtime ( &rawtime ));
		sections[i].end = *(localtime ( &rawtime ));
		schedule_file >> sections[i].start.tm_hour >> sections[i].start.tm_min >>sections[i].end.tm_hour >>sections[i].start.tm_min >> sections[i].num;
	}
	struct tm start;
	start = *(localtime ( &rawtime ));
	start.tm_hour = 8;
	start.tm_min = 0;
	double * sim_result;
	//RUN THE SIMULATION
	cout <<"RUNNING MONTE CARLO SIMULATION\n";
	sim_result = run_monte_carlo_simulation(start,temp_out,temp_in,temp_fsize,sections,class_size);
	//DO REGRESSION
	/*
	double c0,c1;
	double cov00,cov01,cov11,sumsq;
	cout << "RUNNING LINEAR REGRESSION\n";
	gsl_fit_linear (sim_result, 1, temp_in, 1,temp_fsize, &c0, &c1, &cov00, &cov01, &cov11, &sumsq); //USE GNU SCIENTIFIC LIBRARY LINEAR REGRESSION FUNCTION
	double reg_result[temp_fsize];
	for(int i = 0; i < temp_fsize; i++)
	{
		reg_result[i] = c0 * sim_result[i] + c1; //computing indoor through linear regression
	}

	*/

	cout << "RUNNING B-SPLINE\n";
	double control[10];
	double x[11];

	x[0] = 0;
	x[1] = 4;
	x[2] = 9;
	x[3] = 18;
	x[4] = 21;
	x[5] = 25;
	x[6] = 35;
	x[7] = 44;
	x[8] = 51;
	x[9] = 54;
	x[10] = 59;
	double * quadratic = b_spline(0,60,x,temp_in,11,QUADRATIC_SPLINE);
	double * linear = b_spline(0,60,x,temp_in,11,LINEAR_SPLINE);



	cout << "SAVING RESULT IN output.csv\n";
	output_file << "indoor temp, model, linear spline run 1, quadratic spline run 1\n";
	for(int i = 0; i < temp_fsize; i++)
	{
		output_file << temp_in[i] << "," <<sim_result[i]<< ","<< linear[i] << "," << quadratic[i] <<"\n";
	}
	out_file.close();
	in_file.close();
	schedule_file.close();
	output_file.close();
	return 0;
}

double calc_heat_escaped(double temp_in, double temp_out, double time_in_seconds)
{
        temp_in  = fTok(temp_in);
        temp_out = fTok(temp_out);
        double air_escaping = ((height * width * discharge_coefficient) * sqrt(gravity * height * (temp_out - temp_in)/ temp_in)) / 3.0; //this computes the air escaping per second
        double heat_escaping = air_escaping * 1.08 * (temp_out - temp_in); // computes heat in BTU per second
        return (heat_escaping * time_in_seconds); // computes total heat in that time frame
}
double calc_new_temperature(double heat_escaped, double temp_in)
{
        double total_heat = room_volume * air_density * temp_in;
        double new_temperature = (total_heat + heat_escaped - 500)/ (room_volume * air_density);
        return new_temperature;
}
int calc_expected_people_opening_door(struct tm now, struct section list[], int size)
{
        int students = 0;
        for(int i = 0; i < size; i++)
        {
                int start_diff = difftime(mktime(&list[i].start), mktime(&now));
                int end_diff = difftime(mktime(&now), mktime(&list[i].end));
                start_diff/=60;
                end_diff /=60;
                students += floor(normal_pdf(start_diff, -2.0,2.5) * list[i].num); //computes expected number of students going to class
                students += floor(normal_pdf(end_diff, -2.0,2.5) * list[i].num);   //computes expected number of students leaving class
        }
        return students;

}
int calc_expected_time_door_open(int people)
{
        int seconds = people * 5;
        if(seconds > 60)
                return 60;
        return seconds;
}

double* run_monte_carlo_simulation(struct tm start, double out[], double in[], int temp_size, struct section list[], int section_size)
{
        double * estimated_in = new double[temp_size];
        estimated_in[0] = in[0]; //the first temperature is not estimated so we copy it
        for(int i = 0;i < temp_size - 1; i++)
        {
                int people = calc_expected_people_opening_door(start,list,section_size);
                int seconds = calc_expected_time_door_open(people);
                double heat_escaped = calc_heat_escaped(estimated_in[i],out[i],seconds);
                double new_temperature = calc_new_temperature(heat_escaped,estimated_in[0]);
                estimated_in[i + 1] = new_temperature;
        }
        return estimated_in;
}

double normal_pdf(double x, double m, double s)
{
    static const double inv_sqrt_2pi = 0.3989422804014327;
    double a = (x - m) / s;

    return inv_sqrt_2pi / s * std::exp(-0.5f * a * a);
}

double fTok(double f)
{
        return ((f - 32)* 5) / 9 + 273.15;
}

double* b_spline(int start, int end, double x[], double y[], int N, int mode)
{
	double * output;
	double p0 = 0;
	double p1 = 0;
	double p2;
	double p3;
	int t2,t3,t4,t5;
	int c0,c1;
	double t = 0;
	output = new double[end - start + 1];
	for(int i = start; i < end +1; i++)
	{
		between(i,&c0,&c1,x,N); // get the control points where i is in between		
		t = ((double) i - (double) x[c0]) / ((double) x[c1] - (double) x[c0]);
		if(mode ==  1) //if linear b-spline
		{
			t2 = (int) x[c0];
			t3 = (int) x[c1];
			p1 = y[t2];
			p2 = y[t3];
			t = i;
			output[i] = (t3 - t)/(t3-t2) * p1 + (t - t2)/(t3 - t2) * p2;//(1 - t) * p0 + t* p1; 
		}
		else if (mode  == 2) // if quadratic b-spline
		{
			
			t = i;
			t2 = (int) x[c0 - 1];
			t3 = (int) x[c0];
			t4 = (int) x[c1];
			t5 = (int) x[c1 + 1];
			p1 = y[t2];
			p2 = y[t4];
			p3 = y[t5];
			double p01 = (double)(t4 - t)/(double)(t4-t2) * p1 + (double) (t - t2)/(double) (t4 - t2) * p2;
			double p11 = (double)(t5 - t)/(double) (t5-t3) * p2 + (double) (t - t3)/(double) (t5 - t3) * p3;
			output[i] = (double) (t4 - t)/(double)(t4 - t3) * p01 + (double)(t - t3)/(double)(t4 - t3) * p11;
			
		}
	}

	return output;
}
void between(int num, int *a, int *b, double x[], int N)
{
	for(int i = 0; i < N; i++)
	{
		if(num < x[i])
		{
			*a = i-1;
			*b = i;
			return;
		}
	}
}
int nChoosek( int n, int k )
{
    if (k > n) return 0;
    if (k * 2 > n) k = n-k;
    if (k == 0) return 1;

    int result = n;
    for( int i = 2; i <= k; ++i ) {
        result *= (n-i+1);
        result /= i;
    }
    return result;
}
