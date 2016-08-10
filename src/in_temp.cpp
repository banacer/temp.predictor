#include <iostream>
#include <fstream>
#include <cmath>
#include <time.h>
#include <random>
//#include <gsl/gsl_math.h>
//#include <gsl/gsl_fit.h>

#define temp_fsize 120
#define class_size 6
#define height 10.0                          // ft
#define width 8.0                            // ft
#define gravity 32.2                         // ft.s^-2
#define discharge_coefficient 0.6            // 
#define room_volume 30000                    // ft^3
#define air_density 0.0752                   //lbs.ft^-3

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
using namespace std;

int main()
{
	time_t rawtime;
	time ( &rawtime );
	std::default_random_engine generator;
  	std::normal_distribution<double> distribution(0.25,0.13);
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
		schedule_file >> sections[i].start.tm_hour >> sections[i].start.tm_min >>sections[i].end.tm_hour >>sections[i].end.tm_min >> sections[i].num;
	} 
	struct tm start;
	start = *(localtime ( &rawtime ));
	start.tm_hour = 8;
	start.tm_min = 0;
	double * sim_result;
	//RUN THE SIMULATION
	cout <<"RUNNING MONTE CARLO SIMULATION\n";
	sim_result = run_monte_carlo_simulation(start,temp_out,temp_in,temp_fsize,sections,class_size);
	/*
	//DO REGRESSION
	double c0,c1;
	double cov00,cov01,cov11,sumsq;
	cout << "RUNNING LINEAR REGRESSION\n";
	gsl_fit_linear (sim_result, 1, temp_in, 1,temp_fsize, &c0, &c1, &cov00, &cov01, &cov11, &sumsq); //USE GNU SCIENTIFIC LIBRARY LINEAR REGRESSION FUNCTION
	double reg_result[temp_fsize];
	for(int i = 0; i < temp_fsize; i++)
	{
		reg_result[i] = c0 * sim_result[i] + c1; //computing indoor through linear regression
	}

	cout << "SAVING RESULT IN output.csv\n";
	output_file << "indoor temp, model temp, regression temp\n";
	for(int i = 0; i < temp_fsize; i++)
	{
		output_file << temp_in[i] << "," << sim_result[i] << "," << reg_result[i] <<"\n";
	}
	*/
	for(int i = 0; i < temp_fsize; i++)
	{
		double number = distribution(generator);
		cout << sim_result[i] + number << endl;
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
		start_diff/=60;
		int end_diff = difftime(mktime(&now), mktime(&list[i].end));
		end_diff/=60;
		students += floor(fabs(normal_pdf(start_diff, -2.0,2.5)) * list[i].num); //computes expected number of students going to class
		students += floor(fabs(normal_pdf(end_diff, -2.0,2.5)) * list[i].num);   //computes expected number of students leaving class		
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
	estimated_in[0] = 70.21; //the first temperature is not estimated so we copy it
	for(int i = 0;i < temp_size - 1 ; i++)
	{
		int people = calc_expected_people_opening_door(start,list,section_size);
		int seconds = calc_expected_time_door_open(people);
		double heat_escaped = calc_heat_escaped(estimated_in[i],out[i],seconds);
		double new_temperature = calc_new_temperature(heat_escaped,estimated_in[i]);
		estimated_in[i + 1] = new_temperature;
		start.tm_min+=1;
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

