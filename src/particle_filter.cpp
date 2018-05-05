/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

/**
 * init Initializes particle filter by initializing particles to Gaussian
 *   distribution around first position and all the weights to 1.
 * @param x Initial x position [m] (simulated estimate from GPS)
 * @param y Initial y position [m]
 * @param theta Initial orientation [rad]
 * @param std[] Array of dimension 3 [standard deviation of x [m], standard deviation of y [m]
 *   standard deviation of yaw [rad]]
 */
void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	cout << "init" << endl;

	num_particles = 10;

	default_random_engine gen;

	m_std_x = std[0];
	m_std_y = std[1];
	m_std_theta = std[2];


	// create a normal (Gaussian) distribution
	normal_distribution<double> dist_x(x, m_std_x);
	normal_distribution<double> dist_y(y, m_std_y);
	normal_distribution<double> dist_theta(theta, m_std_theta);


	for (int i = 0; i < num_particles; ++i) {

		Particle p;
		p.id = i;

		p.x = dist_x(gen);
		p.y = dist_y(gen);
		p.theta = dist_theta(gen);
		p.weight = 1;

		particles.push_back(p);

		// Print your samples to the terminal.
		cout << "p " << i  << " " << particles[i].x << " " << particles[i].y << " " << particles[i].theta << " " << particles[i].weight << endl;
	}

	is_initialized = true;

}
/**
 * prediction Predicts the state for the next time step
 *   using the process model.
 * @param delta_t Time between time step t and t+1 in measurements [s]
 * @param std_pos[] Array of dimension 3 [standard deviation of x [m], standard deviation of y [m]
 *   standard deviation of yaw [rad]]
 * @param velocity Velocity of car from t to t+1 [m/s]
 * @param yaw_rate Yaw rate of car from t to t+1 [rad/s]
 */
void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	cout << "prediction" << endl;

	for (int i = 0; i < num_particles; i++){
		Particle& p = particles[i];
		//avoid division by zero
		if (fabs(yaw_rate) > 0.001) {
			p.x = p.x + velocity/yaw_rate * ( sin (p.theta + yaw_rate*delta_t) - sin(p.theta));
			p.y = p.y + velocity/yaw_rate * ( cos(p.theta) - cos(p.theta + yaw_rate*delta_t) );
		}
		else {
			p.x = p.x + velocity * delta_t * cos(p.theta);
			p.y = p.y + velocity * delta_t * sin(p.theta);
		}

		// Print your samples to the terminal.
		cout << "p " << i  << " " << particles[i].x << " " << particles[i].y << " " << particles[i].theta << " " << particles[i].weight << endl;
	}

}

/**
 * dataAssociation Finds which observations correspond to which landmarks (likely by using
 *   a nearest-neighbors data association).
 * @param predicted Vector of predicted landmark observations
 * @param observations Vector of landmark observations
 */
void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	for (int i = 0; i < observations.size(); i++) {

		// Initialize with big number.
		double minDistance = 1.0e99;

		int mapId = -1;

		for (int j = 0; j < predicted.size(); j++) {

			double xDistance = observations[i].x - predicted[j].x;
			double yDistance = observations[i].y - predicted[j].y;

			double distance = xDistance * xDistance + yDistance * yDistance;

			if (distance < minDistance) {
				minDistance = distance;
				mapId = predicted[j].id;
			}
		}

		observations[i].id = mapId;
	}

}

/**
 * updateWeights Updates the weights for each particle based on the likelihood of the
 *   observed measurements.
 * @param sensor_range Range [m] of sensor
 * @param std_landmark[] Array of dimension 2 [Landmark measurement uncertainty [x [m], y [m]]]
 * @param observations Vector of landmark observations
 * @param map Map class containing map landmarks
 */
void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html



	for (int p = 0; p < num_particles; p++){

		// find predicted landmarks
		vector<LandmarkObs> predictedLandmarks;
		for(int lm = 0; lm < map_landmarks.landmark_list.size(); lm++) {

			LandmarkObs predictedLandmark;
			predictedLandmark.id 	= map_landmarks.landmark_list[lm].id_i;
			predictedLandmark.x		= map_landmarks.landmark_list[lm].x_f;
			predictedLandmark.y 	= map_landmarks.landmark_list[lm].y_f;

			double dx = predictedLandmark.x - particles[p].x;
			double dy = predictedLandmark.y - particles[p].y;

			if ( dx*dx + dy*dy <= sensor_range*sensor_range) {
				predictedLandmarks.push_back(predictedLandmark);
			}
		}

		vector<LandmarkObs> mappedObservations = map(particles[p], observations);
		dataAssociation(predictedLandmarks, mappedObservations);

		double weight = 1.0;

		for(int i =0; i < mappedObservations.size(); i++){
			weight *= calculateWeight(mappedObservations[i], map_landmarks);
		}

		particles[p].weight = weight;
	}

	normalizeParticlesWeights();


}

double ParticleFilter::calculateWeight(LandmarkObs observation, const Map &map_landmarks){

	double weight = 0;

	if(observation.id < map_landmarks.landmark_list.size()){

		double sig_x= m_std_x;
		double sig_y= m_std_y;
		double x_obs= observation.x;
		double y_obs= observation.y;
		double mu_x= map_landmarks.landmark_list[observation.id].x_f;
		double mu_y= map_landmarks.landmark_list[observation.id].y_f;

		// calculate normalization term
		double gauss_norm= (1.0/(2.0 * M_PI * sig_x * sig_y));

		// calculate exponent
		double exponent= ((x_obs - mu_x)*(x_obs - mu_x))/(2 * sig_x*sig_x) + ((y_obs - mu_y)*(y_obs - mu_y))/(2 * sig_y*sig_y);

		// calculate weight using normalization terms and exponent
		weight = gauss_norm * exp(-exponent);
	}

	if(weight == 0){
		weight = .00001;
	}
	else{
		cout << "weight calculated\n";
	}

	return weight;

}

vector<LandmarkObs> ParticleFilter::map(Particle particle, const std::vector<LandmarkObs> observations){

	vector<LandmarkObs> mappedObservations;

	for (int i = 0; i < observations.size(); i++){

		const LandmarkObs& observation = observations[i];
		// transform to map x coordinate
		double x_map = particle.x + (cos(particle.theta) * observation.x) - (sin(particle.theta) * observation.y);
		// transform to map y coordinate
		double y_map = particle.y  + (sin(particle.theta) * observation.x) + (cos(particle.theta) * observation.y);

		LandmarkObs mappedObservation;
		mappedObservation.id = observation.id;
		mappedObservation.x = x_map;
		mappedObservation.y = y_map;

		mappedObservations.push_back(mappedObservation);
	}

	return mappedObservations;

}

void ParticleFilter::normalizeParticlesWeights(void){

	double weightSum = 0;

	for (int i = 0; i < particles.size(); i++){
		weightSum = particles[i].weight;
	}

	for (int i = 0; i < particles.size(); i++){
		particles[i].weight /=  weightSum;
	}
}

/**
 * resample Resamples from the updated set of particles to form
 *   the new set of particles.
 */
void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	static default_random_engine gen;
	gen.seed(123);
	discrete_distribution<> dist_particles(weights.begin(), weights.end());
	vector<Particle> new_particles;
	new_particles.resize(num_particles);
	for (int i = 0; i < num_particles; i++) {
		new_particles[i] = particles[dist_particles(gen)];
	}
	particles = new_particles;

}

/*
 * Set a particles list of associations, along with the associations calculated world x,y coordinates
 * This can be a very useful debugging tool to make sure transformations are correct and assocations correctly connected
 */
Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
		const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	particle.associations= associations;
	particle.sense_x = sense_x;
	particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
	copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
	string s = ss.str();
	s = s.substr(0, s.length()-1);  // get rid of the trailing space
	return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
	copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
	string s = ss.str();
	s = s.substr(0, s.length()-1);  // get rid of the trailing space
	return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
	copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
	string s = ss.str();
	s = s.substr(0, s.length()-1);  // get rid of the trailing space
	return s;
}
